!!****if* source/physics/IncompNS/IncompNSStats/ins_statsRestressesTimeavg
!!
!! NAME
!!
!!  ins_statsRestressesTimeavg
!!
!! SYNOPSIS
!!
!!  call ins_statsRestressesTimeavg(integer(in) :: n)
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



subroutine ins_statsRestressesTimeavg(n)

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

  real :: uc,vc,wc

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

     ! Loop on cell centers:
     do k=GRID_KLO,GRID_KHI
       do j=GRID_JLO,GRID_JHI
         do i=GRID_ILO,GRID_IHI

           ! Cell center velocities
           uc = 0.5*(facexData(VELC_FACE_VAR,i,j,k)+facexData(VELC_FACE_VAR,i+1,j,k))
           vc = 0.5*(faceyData(VELC_FACE_VAR,i,j,k)+faceyData(VELC_FACE_VAR,i,j+1,k))
           wc = 0.5*(facezData(VELC_FACE_VAR,i,j,k)+facezData(VELC_FACE_VAR,i,j,k+1))

           ! Reynolds Stresses:
           solnData(UUAV_VAR,i,j,k) = re_inp1*(re_n*solnData(UUAV_VAR,i,j,k) + uc*uc)
           solnData(VVAV_VAR,i,j,k) = re_inp1*(re_n*solnData(VVAV_VAR,i,j,k) + vc*vc)
           solnData(WWAV_VAR,i,j,k) = re_inp1*(re_n*solnData(WWAV_VAR,i,j,k) + wc*wc)
           solnData(UVAV_VAR,i,j,k) = re_inp1*(re_n*solnData(UVAV_VAR,i,j,k) + uc*vc)
           solnData(UWAV_VAR,i,j,k) = re_inp1*(re_n*solnData(UWAV_VAR,i,j,k) + uc*wc)
           solnData(VWAV_VAR,i,j,k) = re_inp1*(re_n*solnData(VWAV_VAR,i,j,k) + vc*wc)

           ! Pressure:
           solnData(PPAV_VAR,i,j,k) = re_inp1*(re_n*solnData(PPAV_VAR,i,j,k) + &
                                      solnData(PRES_VAR,i,j,k)*solnData(PRES_VAR,i,j,k))
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


