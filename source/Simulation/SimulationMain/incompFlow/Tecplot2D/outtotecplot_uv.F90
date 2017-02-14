!!****if* source/Simulation/SimulationMain/incompFlow/Tecplot2D/outtotecplot_uv
!!
!! NAME
!!
!!  outtotecplot_uv
!!
!! SYNOPSIS
!!
!!  call outtotecplot_uv(integer(in) :: mype,
!!                       real(in) :: time,
!!                       real(in) :: dt,
!!                       integer(in) :: istep,
!!                       integer(in) :: count,
!!                       real(in) :: timer,
!!                       integer(in) :: blocklist,
!!                       integer(in) :: blockcount,
!!                       integer(in) :: firstfileflag)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   mype : my Processor Number
!!
!!   time : 
!!
!!   dt : 
!!
!!   istep : 
!!
!!   count : 
!!
!!   timer : 
!!
!!   blocklist : 
!!
!!   blockcount : 
!!
!!   firstfileflag : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

! Subroutine outtotecplot
!
! Subroutine to write out to Tecplot data in binary form.
!
! ---------------------------------------------------------------------------
#include "constants.h"
#include "Flash.h"

subroutine outtotecplot_uv(mype,time,dt,istep,count, &
           timer,blockList,blockCount,firstfileflag)

   use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
       Grid_getBlkBoundBox, Grid_getBlkCenterCoords

#ifdef FLASH_GRID_UG
#else
   use physicaldata, ONLY : interp_mask_unk,interp_mask_unk_res
#endif

   use IncompNS_data, ONLY : ins_invRe

   implicit none
   include "Flash_mpi.h"

   integer, intent(in) :: mype,istep,count,firstfileflag
   integer, intent(in) :: blockCount
   integer, intent(in) :: blockList(MAXBLOCKS)
   real, intent(in) :: time,dt,timer

   ! Local Variables:
   integer :: numblocks,var,i,j,k,lb,nxc,nyc,nzc
   character(25) :: filename
   character(6) :: index_lb,index_mype

   real xedge(NXB+1),xcell(NXB+1)
   real yedge(NYB+1),ycell(NYB+1)
   real intsx(NXB+1), intsy(NYB+1)

   real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData

   real facevarxx(NXB+2*NGUARD+1,NYB+2*NGUARD),   &
        facevarxxan(NXB+2*NGUARD+1,NYB+2*NGUARD), &
        difffacevarx(NXB+2*NGUARD+1,NYB+2*NGUARD)

   real facevaryy(NXB+2*NGUARD,NYB+2*NGUARD+1),   &
        facevaryyan(NXB+2*NGUARD,NYB+2*NGUARD+1), &
        difffacevary(NXB+2*NGUARD,NYB+2*NGUARD+1)

   real presvar(NXB+2*NGUARD,NYB+2*NGUARD),       &
        presvar_an(NXB+2*NGUARD,NYB+2*NGUARD),    &
        diffpresvar(NXB+2*NGUARD,NYB+2*NGUARD)

   real xi,yj

   real*4 arraylb(NXB+1,NYB+2,1),        &
          arraylb2(NXB+2,NYB+1,1),       &
          arraylb3(NXB,NYB,1)

   integer blockID
   real del(MDIM),dx,dy
   real, dimension(MDIM)  :: coord,bsize
   real ::  boundBox(2,MDIM)

   integer*4 TecIni,TecDat,TecZne,TecNod,TecFil,TecEnd
   integer*4 VIsdouble
   integer*4 Debug,ijk,ijk2,ijk3,Npts,NElm
   character*1 NULLCHR

   real :: pi, tn, mvisc

   mvisc = ins_invRe

   pi = 4.*atan(1.)

   tn = time

#ifdef WRITE_TECPLOT
!-----------------------------------------------------------------------
!                                                         TecPlot set-up
!-----------------------------------------------------------------------
   Debug     = 0
   VIsdouble = 0
   NULLCHR   = CHAR(0)
   ijk       = (NXB+1)*(NYB+2)
   ijk2      = (NXB+2)*(NYB+1)
   ijk3      = NXB*NYB
!-----------------------------------------------------------------------


! -- filetime.XX --

   write(filename, '("./IOData/dauv_time.", i2.2)') mype

   ! create/clear filetime.XX if time = 0
   if(firstfileflag .eq. 0) then
      open(unit=33, file=filename, status='replace')
#ifdef FLASH_GRID_UG
      write(33,*) 'NXB, NYB, NZB'
      write(33,'(3i4.1)') NXB, NYB, NZB
#else
      write(33,*) 'NXB, NYB, NZB, interp. order (prolong, restrict)'
      write(33,'(5i4.1)') NXB, NYB, NZB, interp_mask_unk(1), &
              interp_mask_unk_res(1)
#endif
      write(33,'(a23,a43,a49,a12)') 'file number, time, dt, ',      &
               'step number, ',                                     &
               'total number of blocks, number of blocks output, ', &
               'elapsed time'
      close(33)
   endif

   ! write timestep data to filetime.XX on each processor
   open(unit=33, file=filename, status='old', position='append')
   write(33,66) count, tn, dt, istep,blockcount,timer
   close(33)
66    format(i4.4,g23.15,g23.15,i8.1,i5.1,g23.15)

! -- data.XXXX.XX --

   nxc = NXB + NGUARD + 1
   nyc = NYB + NGUARD + 1

   ! write solution data to UVEL.XXXX.XX
   write(filename,'("./IOData/UVEL.",i4.4,".",i2.2,".plt")') &
         count, mype

   i = TecIni('AMR2D'//NULLCHR,              &
              'x y unum'//NULLCHR,           &
              filename//NULLCHR,             &
              './IOData/'//NULLCHR,          &
              Debug,VIsdouble)
   call int2char(mype,index_mype)

   do lb = 1,blockcount

      blockID =  blockList(lb)

      ! Get blocks dx, dy ,dz:
      call Grid_getDeltas(blockID,del)
      dx = del(IAXIS)
      dy = del(JAXIS)
      ! Get Coord and Bsize for the block:
      ! Bounding box:
      call Grid_getBlkBoundBox(blockId,boundBox)
      bsize(:) = boundBox(2,:) - boundBox(1,:)

      call Grid_getBlkCenterCoords(blockId,coord)

      ! Point to blocks center and face-x vars:
      call Grid_getBlkPtr(blockID,facexData,FACEX)

      do j=NGUARD,nyc
         do i=NGUARD+1,nxc
            if (j .eq. NGUARD) then

               facevarxx(i,j) =  0.5*(facexData(VELC_FACE_VAR,i,j,1)+    &
                                      facexData(VELC_FACE_VAR,i,j+1,1))

            elseif (j .eq. nyc) then

               facevarxx(i,j) =  0.5*(facexData(VELC_FACE_VAR,i,j,1)+    &
                                      facexData(VELC_FACE_VAR,i,j-1,1))
            else

               facevarxx(i,j) = facexData(VELC_FACE_VAR,i,j,1)

            endif
         enddo
      enddo

      ! Write Block Results into data file:
      call int2char(blockID,index_lb)

      i = TecZne(                                                  &
          'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,    &
           NXB+1,NYB+2,1,                                          &
          'BLOCK'//NULLCHR,                                        &
           CHAR(0))

      ! Write x:
      do j=1,NYB+2
         do i=1,NXB+1
            xi =  coord(IAXIS) - bsize(IAXIS)/2.0 + real(i-1)*dx
            arraylb(i,j,1) = sngl(xi)
         enddo
      enddo
      i = TecDat(ijk,arraylb,0)

      ! Write y:
      do j=1,NYB+2

         if (j .eq. 1) then
            yj = coord(JAXIS) - bsize(JAXIS)/2.0
         elseif (j .eq. NYB+2) then
            yj = coord(JAXIS) + bsize(JAXIS)/2.0
         else
            yj = coord(JAXIS) - bsize(JAXIS)/2.0 + real(j-2)*dy + dy/2.0
         endif

         do i=1,NXB+1
            arraylb(i,j,1) = sngl(yj)
         enddo
      enddo
      i = TecDat(ijk,arraylb,0)

      ! Write unum:
      arraylb(:,:,1) =    &
      sngl(facevarxx(NGUARD+1:nxc,NGUARD:nyc))
      i = TecDat(ijk,arraylb,0)

   enddo
   i = TecEnd()

   ! write solution data to VVEL.XXXX.XX
   write(filename,'("./IOData/VVEL.",i4.4,".",i2.2,".plt")') &
         count, mype

   i = TecIni('AMR2D'//NULLCHR,        &
       'x y vnum'//NULLCHR,            &
       filename//NULLCHR,              &
       './IOData/'//NULLCHR,           &
       Debug,VIsdouble)

   call int2char(mype,index_mype)

   do lb = 1,blockcount

      blockID =  blockList(lb)

      ! Get blocks dx, dy ,dz:
      call Grid_getDeltas(blockID,del)
      dx = del(IAXIS)
      dy = del(JAXIS)
      ! Get Coord and Bsize for the block:
      ! Bounding box:
      call Grid_getBlkBoundBox(blockId,boundBox)
      bsize(:) = boundBox(2,:) - boundBox(1,:)

      call Grid_getBlkCenterCoords(blockId,coord)

      ! Point to blocks center and face-y vars:
      call Grid_getBlkPtr(blockID,faceyData,FACEY)

      do j=NGUARD+1,nyc
         do i=NGUARD,nxc
            if (i .eq. NGUARD) then

               facevaryy(i,j) =  0.5*(faceyData(VELC_FACE_VAR,i,j,1)+    &
                                      faceyData(VELC_FACE_VAR,i+1,j,1))

            elseif (i .eq. nxc) then

               facevaryy(i,j) =  0.5*(faceyData(VELC_FACE_VAR,i,j,1)+   &
                                      faceyData(VELC_FACE_VAR,i-1,j,1))

            else

               facevaryy(i,j) = faceyData(VELC_FACE_VAR,i,j,1)

            endif
         enddo
      enddo

      ! Write Block Results into data file:
      call int2char(blockID,index_lb)

      i = TecZne(                                                 &
          'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,   &
          NXB+2,NYB+1,1,                                          &
          'BLOCK'//NULLCHR,                                       &
          CHAR(0))

      ! Write x:
      do j=1,NYB+1
         do i=1,NXB+2 !NGUARD
            if (i .eq. 1) then
               xi = coord(IAXIS) - bsize(IAXIS)/2.0
            elseif (i .eq. NXB+2) then !NGUARD
               xi = coord(IAXIS) + bsize(IAXIS)/2.0
            else
               xi = coord(IAXIS) - bsize(IAXIS)/2.0 +  &
                    real(i-2)*dx + dx/2.0
            endif
            arraylb2(i,j,1) = sngl(xi)
         enddo
      enddo
      i = TecDat(ijk2,arraylb2,0)

      ! Write y:
      do j=1,NYB+1
         yj =  coord(JAXIS) - bsize(JAXIS)/2.0 +  &
               real(j-1)*dy
         do i=1,NXB+2 !NGUARD
            arraylb2(i,j,1) = sngl(yj)
         enddo
      enddo
      i = TecDat(ijk2,arraylb2,0)

      ! Write vnum:
      arraylb2(:,:,1) = &
      sngl(facevaryy(NGUARD:nxc,NGUARD+1:nyc))
      i = TecDat(ijk2,arraylb2,0)

   enddo
   i = TecEnd()

   ! write solution data to PRES.XXXX.XX
   write(filename,'("./IOData/PRES.",i4.4,".",i2.2,".plt")') &
         count, mype

   i = TecIni('AMR2D'//NULLCHR,        &
       'x y pnum'//NULLCHR,   &
       filename//NULLCHR,              &
       './IOData/'//NULLCHR,           &
       Debug,VIsdouble)

   call int2char(mype,index_mype)

   do lb = 1,blockcount

      blockID =  blockList(lb)

      ! Get blocks dx, dy ,dz:
      call Grid_getDeltas(blockID,del)
      dx = del(IAXIS)
      dy = del(JAXIS)
      ! Get Coord and Bsize for the block:
      ! Bounding box:
      call Grid_getBlkBoundBox(blockId,boundBox)
      bsize(:) = boundBox(2,:) - boundBox(1,:)

      call Grid_getBlkCenterCoords(blockId,coord)

      ! Point to blocks center and center vars:
      call Grid_getBlkPtr(blockID,solnData,CENTER)

      do j=NGUARD+1,nyc-1
         do i=NGUARD+1,nxc-1
             presvar(i,j) = solnData(PRES_VAR,i,j,1)
          enddo
       enddo

       ! Write Block Results into data file:
       call int2char(blockID,index_lb)

       i = TecZne(                                                  &
           'ZONE T=BLKPROC'//index_lb//'.'//index_mype//NULLCHR,    &
           NXB,NYB,1,                                               &
           'BLOCK'//NULLCHR,                                        &
           CHAR(0))

       ! Write x:
       do j=1,NYB
          do i=1,NXB
             xi = coord(IAXIS) - bsize(IAXIS)/2.0 + &
                  real(i-1)*dx + dx/2.0
             arraylb3(i,j,1) = sngl(xi)
          enddo
       enddo
       i = TecDat(ijk3,arraylb3,0)

       ! Write y:
       do j=1,NYB
          yj =  coord(JAXIS) - bsize(JAXIS)/2.0 +   &
                real(j-1)*dy + dy/2.0
          do i=1,NXB
             arraylb3(i,j,1) = sngl(yj)
          enddo
       enddo
       i = TecDat(ijk3,arraylb3,0)

       ! Write pnum:
       arraylb3(:,:,1) =  &
       sngl(presvar(NGUARD+1:nxc-1,NGUARD+1:nyc-1))
       i = TecDat(ijk3,arraylb3,0)
    enddo
    i = TecEnd()

#endif

   return

End subroutine outtotecplot_uv
