subroutine amr_cp_loc(idest,isource,iface,iopt,idiag,nlayers)


! $RCSfile: amr_cp_loc.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:08 $

!------------------------------------------------------------------------
!
! This routine copies guard cell information to face iface in block
! idest, from the appropriate face of the neighboring block, assuming
! that both blocks are on the same processor.
! It can be easily edited to alter the data pattern required for schemes
! of different order.
!
! Arguments:
!      idest id of the destination block which requires guard 
!       cell info
!      isource id of the block from which guard cell info is to
!       be obtained
!      iface the id of the face of the block which requires
!       guard cell info
!       iopt            a switch to control which data source is to be used
!                       iopt=1 will use 'unk'
!                       iopt=2 will use 'work'
!       nlayers         the number of guard cell layers at each boundary
!
!
! Written :     Peter MacNeice          February 1997
!------------------------------------------------------------------------
  use physicaldata
  use workspace
  implicit none

  
      integer idest,isource,iface,iopt,idiag,nlayers

      integer istart,jstart,kstart
      integer iend,jend,kend
      integer i_dest,i_source,j_dest,j_source,k_dest,k_source
      integer id,is,ivar
      integer i,j,k


#ifdef TIMINGS
#include "timer.fh"
#endif

!------------------------------------------------------------------------

#ifdef TIMINGS
      itimer1 = irtc()
#endif

      kstart = 1
      kend   = ku_bnd
      jstart = 1
      jend   = ju_bnd
      istart = 1
      iend   = iu_bnd

      if (iface.eq.1.or.iface.eq.2) then
        kstart = nzlo
        kend = nzhi
        jstart = nylo
        jend = nyhi
        istart = 0
        iend = nlayers - 1
      else if (iface.eq.3.or.iface.eq.4) then
        istart = 1
        iend = iu_bnd
        kstart = nzlo
        kend = nzhi
        jstart = 0
        jend = nlayers - 1
      else if(iface.eq.5.or.iface.eq.6) then
        istart = 1
        iend = iu_bnd
        jstart = 1
        jend = ju_bnd
        kstart = 0
        kend = nlayers - 1
      end if

!------------------------------------
      if(iopt.eq.1) then
!------------------------------------

      if(iface.eq.2) then

       i_dest   = nxb+1+nguard
       i_source = 1+nguard+gc_off_x

       do k = kstart,kend
          do j = jstart,jend
             do i=istart,iend
                id = i_dest+i
                is = i_source+i
                do ivar=1,nvar
                   unk(ivar,id,j,k,idest) = & 
     &                  unk(ivar,is,j,k,isource)
                enddo
             enddo
          enddo
       enddo

      elseif(iface.eq.1) then

        i_dest   = nguard
        i_source = nxb+nguard-gc_off_x
        
       do k = kstart,kend
          do j = jstart,jend
             do i=istart,iend
                id = i_dest-i
                is = i_source-i
                do ivar=1,nvar
                   unk(ivar,id,j,k,idest) = & 
     &                  unk(ivar,is,j,k,isource)
                enddo
             enddo
          enddo
       enddo

      elseif(iface.eq.4) then

       j_dest   = nyb+1+nguard
       j_source = 1+nguard+gc_off_y

       do k = kstart,kend
          do j=jstart,jend
             id = j_dest+j
             is = j_source+j
             do i = istart,iend
                do ivar=1,nvar
                   unk(ivar,i,id,k,idest) = & 
     &                  unk(ivar,i,is,k,isource)
                enddo
             end do
          end do
       end do

      elseif(iface.eq.3) then

       j_dest   = nguard
       j_source = nyb+nguard-gc_off_y

       do k = kstart,kend
          do j=jstart,jend
             id = j_dest-j
             is = j_source-j
             do i = istart,iend
                do ivar=1,nvar
                   unk(ivar,i,id,k,idest) = & 
     &                  unk(ivar,i,is,k,isource)
                enddo
             enddo
          end do
       end do

      elseif(iface.eq.6) then

       k_dest   = nzb+1+nguard
       k_source = 1+nguard+gc_off_z

       do k=kstart,kend
         id = k_dest+k
         is = k_source+k
         do j = jstart,jend
            do i = istart,iend
               do ivar=1,nvar
                  unk(ivar,i,j,id,idest) = & 
     &                 unk(ivar,i,j,is,isource)
               enddo
            enddo
         enddo
      enddo
        
      elseif(iface.eq.5) then

        k_dest   = nguard
        k_source = nzb+nguard-gc_off_z
        
        do k=kstart,kend
           id = k_dest-k
           is = k_source-k
           do j = jstart,jend
              do i = istart,iend
                 do ivar=1,nvar
                    unk(ivar,i,j,id,idest) = & 
     &                   unk(ivar,i,j,is,isource)
                 enddo
              enddo
           enddo
        enddo
        
      endif


!------------------------------------
      elseif(iopt.eq.2) then
!------------------------------------
 

      if(iface.eq.2) then

       i_dest   = nxb+1+nguard_work
       i_source = 1+nguard_work+gc_off_x
       
       do i=0,nlayers-1
            work(i_dest+i,:,:,idest,1) = & 
     &       work(i_source+i,:,:,isource,1)
       enddo

      elseif(iface.eq.1) then

       i_dest   = nguard_work
       i_source = nxb+nguard_work-gc_off_x

       do i=0,nlayers-1
            work(i_dest-i,:,:,idest,1) = & 
     &       work(i_source-i,:,:,isource,1)
       enddo

      elseif(iface.eq.4) then

       j_dest   = nyb+1+nguard_work
       j_source = 1+nguard_work+gc_off_y

       do j=0,nlayers-1
            work(:,j_dest+j,:,idest,1) = & 
     &       work(:,j_source+j,:,isource,1)
       enddo

      elseif(iface.eq.3) then

       j_dest   = nguard_work
       j_source = nyb+nguard_work-gc_off_y

       do j=0,nlayers-1
            work(:,j_dest-j,:,idest,1) = & 
     &       work(:,j_source-j,:,isource,1)
       enddo

      elseif(iface.eq.6) then

       k_dest   = nzb+1+nguard_work
       k_source = 1+nguard_work+gc_off_z

       do k=0,nlayers-1
            work(:,:,k_dest+k,idest,1) = & 
     &       work(:,:,k_source+k,isource,1)
       enddo

      elseif(iface.eq.5) then

       k_dest   = nguard_work
       k_source = nzb+nguard_work-gc_off_z

       do k=0,nlayers-1
            work(:,:,k_dest-k,idest,1) = & 
     &       work(:,:,k_source-k,isource,1)
       enddo

      endif

!------------------------------------
      endif
!------------------------------------

#ifdef TIMINGS
      itimer2 = irtc()
      irtc_cploc = itimer2-itimer1+irtc_cploc
#endif


      return
      end
