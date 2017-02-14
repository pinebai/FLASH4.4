!!****if* source/physics/IncompNS/IncompNSStats/ins_statsVelpTimeAvg
!!
!! NAME
!!
!!  ins_statsVelpTimeAvg
!!
!! SYNOPSIS
!!
!!  call ins_statsVelpTimeAvg(integer(in) :: n)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   n : 
!!
!! AUTOGENROBODOC
!!
!!
!!***



subroutine ins_statsVelpTimeAvg(n)

  use Grid_interface, ONLY : Grid_getListOfBlocks,Grid_getBlkPtr,Grid_releaseBlkPtr

  implicit none
#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: n

  ! Local Variables
  real :: re_n,re_inp1

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,blockID,i,j,k


  re_n    = real(n)
  re_inp1 = 1./(re_n+1.)

  ! Get List of leaf Blocks:
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! Loop on Blocks
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)


     ! Loop on U velocities:
     do k=GRID_KLO,GRID_KHI
       do j=GRID_JLO,GRID_JHI
         do i=GRID_ILO,GRID_IHI+1
           facexData(VAVG_FACE_VAR,i,j,k) = re_inp1*(re_n*facexData(VAVG_FACE_VAR,i,j,k) + &
                                                       1.*facexData(VELC_FACE_VAR,i,j,k))
         enddo
       enddo
     enddo

     ! Loop on V velocities:
     do k=GRID_KLO,GRID_KHI
       do j=GRID_JLO,GRID_JHI+1
         do i=GRID_ILO,GRID_IHI
           faceyData(VAVG_FACE_VAR,i,j,k) = re_inp1*(re_n*faceyData(VAVG_FACE_VAR,i,j,k) + &
                                                       1.*faceyData(VELC_FACE_VAR,i,j,k))
         enddo
       enddo
     enddo

     ! Loop on W velocities:
     do k=GRID_KLO,GRID_KHI+1
       do j=GRID_JLO,GRID_JHI
         do i=GRID_ILO,GRID_IHI
           facezData(VAVG_FACE_VAR,i,j,k) = re_inp1*(re_n*facezData(VAVG_FACE_VAR,i,j,k) + &
                                                       1.*facezData(VELC_FACE_VAR,i,j,k))
         enddo
       enddo
     enddo
     ! Use DUST var to store average W vel, so to compute domain average later:
     solnData(DUST_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI) =              &
     0.5*(facezData(VAVG_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI) +   &
          facezData(VAVG_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO+1:GRID_KHI+1))

     ! Loop on pressure:
     do k=GRID_KLO,GRID_KHI
       do j=GRID_JLO,GRID_JHI
         do i=GRID_ILO,GRID_IHI
           solnData(PAVG_VAR,i,j,k) = re_inp1*(re_n*solnData(PAVG_VAR,i,j,k) + &
                                                 1.*solnData(PRES_VAR,i,j,k))
         enddo
       enddo
     enddo

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo

  return

end subroutine


