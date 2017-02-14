!!****if* source/physics/IncompNS/IncompNSExtras/IncompNS_velomg2center
!!
!! NAME
!!
!!  IncompNS_velomg2center
!!
!! SYNOPSIS
!!
!!  call IncompNS_velomg2center(integer(in) :: blockList,
!!                              integer(in) :: blockCount)
!!
!! DESCRIPTION
!!
!!   Compute cell-centered (or volume-averaged) versions of some
!!   quantities, by interpolation or averaging from face-centered
!!   versions.
!!
!! ARGUMENTS
!!
!!   blockList : list of block IDs
!!
!!   blockCount : number of block IDs in the list
!!
!! AUTOGENROBODOC
!!
!! NOTES
!!
!!  To write cell centered velocities and vorticity velx,vely,velz,omgx,omgy,omgz,
!!  add REQUIRES physics/IncompNS/IncompNSExtras to Config.
!!
!!***

subroutine IncompNS_velomg2center(blockList,blockCount)


  use Grid_interface, only : Grid_getDeltas,         &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr


  implicit none
#include "constants.h"
#include "Flash.h"

  integer,intent(in) :: blockCount
  integer,intent(in) :: blockList(MAXBLOCKS)

  integer :: blkID

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezdata

  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: omgvertz
#if NDIM == MDIM
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: omgvertx
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: omgverty   
#endif

  real :: del(MDIM)

  integer :: lb,i,j,k

  do lb = 1,blockCount

      blkID = blockList(lb)
  
      ! Get Deltas for blkID
      call Grid_getDeltas(blkID,del)

      ! Get Data in cell center, facex and facey
      call Grid_getBlkPtr(blkID,solnData,CENTER)
      call Grid_getBlkPtr(blkID,facexData,FACEX)
      call Grid_getBlkPtr(blkID,faceyData,FACEY)
#if NDIM == MDIM
      call Grid_getBlkPtr(blkID,facezData,FACEZ)
#endif

      ! Linear interpolation of face velocities to cell centered velx,vely, velz
      ! feed into solnData(VELX_VAR,i,j,1), solnData(VELY_VAR,i,j,1)
      solnData(VELX_VAR,GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC) = &
                                 0.5*(facexData(VELC_FACE_VAR,GRID_ILO_GC  :GRID_IHI_GC,     &  
                                                              GRID_JLO_GC  :GRID_JHI_GC,     &
                                                              GRID_KLO_GC  :GRID_KHI_GC)   + & 
                                      facexData(VELC_FACE_VAR,GRID_ILO_GC+1:GRID_IHI_GC+1,   &           
                                                              GRID_JLO_GC  :GRID_JHI_GC,     &
                                                              GRID_KLO_GC  :GRID_KHI_GC)) 

      solnData(VELY_VAR,GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC) = &
                                 0.5*(faceyData(VELC_FACE_VAR,GRID_ILO_GC  :GRID_IHI_GC,     &  
                                                              GRID_JLO_GC  :GRID_JHI_GC,     &
                                                              GRID_KLO_GC  :GRID_KHI_GC)   + &
                                      faceyData(VELC_FACE_VAR,GRID_ILO_GC  :GRID_IHI_GC,     &  
                                                              GRID_JLO_GC+1:GRID_JHI_GC+1,   &
                                                              GRID_KLO_GC  :GRID_KHI_GC))

#if NDIM == MDIM
      solnData(VELZ_VAR,GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC) = &
                                 0.5*(facezData(VELC_FACE_VAR,GRID_ILO_GC  :GRID_IHI_GC,     &  
                                                              GRID_JLO_GC  :GRID_JHI_GC,     &
                                                              GRID_KLO_GC  :GRID_KHI_GC)   + &
                                      facezData(VELC_FACE_VAR,GRID_ILO_GC  :GRID_IHI_GC,     &  
                                                              GRID_JLO_GC  :GRID_JHI_GC,     &
                                                              GRID_KLO_GC+1:GRID_KHI_GC+1))
#endif 

 
      ! Compute omgvert on z direction, minimum number of NGUARD == 2.
      omgvertz(GRID_ILO_GC+1:GRID_IHI_GC,GRID_JLO_GC+1:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC) = &
              (faceyData(VELC_FACE_VAR,GRID_ILO_GC+1:GRID_IHI_GC,             &
                                       GRID_JLO_GC+1:GRID_JHI_GC,             &
                                       GRID_KLO_GC  :GRID_KHI_GC)   -         &
               faceyData(VELC_FACE_VAR,GRID_ILO_GC  :GRID_IHI_GC-1,           &
                                       GRID_JLO_GC+1:GRID_JHI_GC,             &
                                       GRID_KLO_GC  :GRID_KHI_GC))/del(IAXIS) &
             -(facexData(VELC_FACE_VAR,GRID_ILO_GC+1:GRID_IHI_GC,             &
                                       GRID_JLO_GC+1:GRID_JHI_GC,             &
                                       GRID_KLO_GC  :GRID_KHI_GC)   -         &
               facexData(VELC_FACE_VAR,GRID_ILO_GC+1:GRID_IHI_GC,             &
                                       GRID_JLO_GC  :GRID_JHI_GC-1,           &
                                       GRID_KLO_GC  :GRID_KHI_GC))/del(JAXIS)

      ! Bilinear interpolation to cell centered omgz -> feed into
      ! solndata(OMGZ_VAR,i,j,1)
      solnData(OMGZ_VAR,GRID_ILO_GC+1:GRID_IHI_GC-1,GRID_JLO_GC+1:GRID_JHI_GC-1,GRID_KLO_GC:GRID_KHI_GC)= &
         0.25*(omgvertz(GRID_ILO_GC+1:GRID_IHI_GC-1,GRID_JLO_GC+1:GRID_JHI_GC-1,GRID_KLO_GC:GRID_KHI_GC)+ &
               omgvertz(GRID_ILO_GC+2:GRID_IHI_GC  ,GRID_JLO_GC+1:GRID_JHI_GC-1,GRID_KLO_GC:GRID_KHI_GC)+ &
               omgvertz(GRID_ILO_GC+1:GRID_IHI_GC-1,GRID_JLO_GC+2:GRID_JHI_GC  ,GRID_KLO_GC:GRID_KHI_GC)+ &
               omgvertz(GRID_ILO_GC+2:GRID_IHI_GC  ,GRID_JLO_GC+2:GRID_JHI_GC  ,GRID_KLO_GC:GRID_KHI_GC))

#if NDIM == MDIM

      ! Compute omgvert on x direction, minimum number of NGUARD == 2.
      omgvertx(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC+1:GRID_JHI_GC,GRID_KLO_GC+1:GRID_KHI_GC) = &
              (facezData(VELC_FACE_VAR,GRID_ILO_GC  :GRID_IHI_GC,               &
                                       GRID_JLO_GC+1:GRID_JHI_GC,               &
                                       GRID_KLO_GC+1:GRID_KHI_GC)   -           &
               facezData(VELC_FACE_VAR,GRID_ILO_GC  :GRID_IHI_GC,               &
                                       GRID_JLO_GC  :GRID_JHI_GC-1,             &
                                       GRID_KLO_GC+1:GRID_KHI_GC))/del(JAXIS)   &
             -(faceyData(VELC_FACE_VAR,GRID_ILO_GC  :GRID_IHI_GC,               &
                                       GRID_JLO_GC+1:GRID_JHI_GC,               &
                                       GRID_KLO_GC+1:GRID_KHI_GC)   -           &
               faceyData(VELC_FACE_VAR,GRID_ILO_GC  :GRID_IHI_GC,               &
                                       GRID_JLO_GC+1:GRID_JHI_GC,               &
                                       GRID_KLO_GC  :GRID_KHI_GC-1))/del(KAXIS)      

      ! Bilinear interpolation to cell centered omgx -> feed into
      ! solndata(OMGX_VAR,i,j,1)
      solnData(OMGX_VAR,GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC+1:GRID_JHI_GC-1,GRID_KLO_GC+1:GRID_KHI_GC-1)= &
         0.25*(omgvertx(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC+1:GRID_JHI_GC-1,GRID_KLO_GC+1:GRID_KHI_GC-1)+ &
               omgvertx(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC+1:GRID_JHI_GC-1,GRID_KLO_GC+2:GRID_KHI_GC)  + &
               omgvertx(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC+2:GRID_JHI_GC  ,GRID_KLO_GC+1:GRID_KHI_GC-1)+ &
               omgvertx(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC+2:GRID_JHI_GC  ,GRID_KLO_GC+2:GRID_KHI_GC))


      ! Compute omgvert on y direction, minimum number of NGUARD == 2.
      omgverty(GRID_ILO_GC+1:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC+1:GRID_KHI_GC) = &
              (facexData(VELC_FACE_VAR,GRID_ILO_GC+1:GRID_IHI_GC,             &
                                       GRID_JLO_GC  :GRID_JHI_GC,               &
                                       GRID_KLO_GC+1:GRID_KHI_GC)   -         &
               facexData(VELC_FACE_VAR,GRID_ILO_GC+1:GRID_IHI_GC,             &
                                       GRID_JLO_GC  :GRID_JHI_GC,               &
                                       GRID_KLO_GC  :GRID_KHI_GC-1))/del(KAXIS) &
             -(facezData(VELC_FACE_VAR,GRID_ILO_GC+1:GRID_IHI_GC,             &
                                       GRID_JLO_GC  :GRID_JHI_GC,               &
                                       GRID_KLO_GC+1:GRID_KHI_GC)   -         &
               facezData(VELC_FACE_VAR,GRID_ILO_GC  :GRID_IHI_GC-1,             &
                                       GRID_JLO_GC  :GRID_JHI_GC,               &
                                       GRID_KLO_GC+1:GRID_KHI_GC))/del(IAXIS) 


      ! Bilinear interpolation to cell centered omgy -> feed into
      ! solndata(OMGY_VAR,i,j,1)
      solnData(OMGY_VAR,GRID_ILO_GC+1:GRID_IHI_GC-1,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC+1:GRID_KHI_GC-1)= &
         0.25*(omgvertx(GRID_ILO_GC+1:GRID_IHI_GC-1,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC+1:GRID_KHI_GC-1)+ &
               omgvertx(GRID_ILO_GC+1:GRID_IHI_GC-1,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC+2:GRID_KHI_GC)  + &
               omgvertx(GRID_ILO_GC+2:GRID_IHI_GC  ,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC+1:GRID_KHI_GC-1)+ &
               omgvertx(GRID_ILO_GC+2:GRID_IHI_GC  ,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC+2:GRID_KHI_GC))

#endif

      ! Release data pointers for cell center, facex, facey
      call Grid_releaseBlkPtr(blkID,solnData,CENTER)
      call Grid_releaseBlkPtr(blkID,facexData,FACEX)
      call Grid_releaseBlkPtr(blkID,faceyData,FACEY)
#if NDIM == MDIM
      call Grid_releaseBlkPtr(blkID,facezData,FACEZ)
#endif

  enddo



  ! Compute cell centerd vorticity magnitude:
  do lb = 1,blockCount
     blkID = blockList(lb)
     ! Get Blocks internal limits indexes:
     call Grid_getBlkPtr(blkID,solnData,CENTER)
     solnData(OMGM_VAR,:,:,:) = 0.
#if NDIM == MDIM
     do k=GRID_KLO_GC+1,GRID_KHI_GC-1 
        do j=GRID_JLO_GC+1,GRID_JHI_GC-1
           do i=GRID_ILO_GC+1,GRID_IHI_GC-1
             solnData(OMGM_VAR,i,j,k)=sqrt(solnData(OMGX_VAR,i,j,k)**2. + &
                                           solnData(OMGY_VAR,i,j,k)**2. + &
                                           solnData(OMGZ_VAR,i,j,k)**2.)
           enddo
        enddo
     enddo
#else
     k = CONSTANT_ONE
     do j=GRID_JLO_GC+1,GRID_JHI_GC-1
        do i=GRID_ILO_GC+1,GRID_IHI_GC-1
           solnData(OMGM_VAR,i,j,k)=sqrt(solnData(OMGZ_VAR,i,j,k)**2.)
        enddo
     enddo
#endif
     call Grid_releaseBlkPtr(blkID,solnData,CENTER)   
  enddo


  return

end subroutine IncompNS_velomg2center
