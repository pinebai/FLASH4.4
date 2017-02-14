!!****if* source/Grid/GridSolvers/IsoBndMultipole/gr_isoFindMassCenter
!!
!! NAME
!!
!!  gr_isoFindMassCenter
!!
!! SYNOPSIS
!!
!!  gr_isoFindMassCenter(integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes the centroid of the specified variable and returns
!!  its location in the mpole_common variables Xcm, Ycm, and Zcm.
!!  Also computes the total value of the quantity and leaves it
!!  in the variable Mtot.  If Mtot=0, the routine aborts.
!!
!! ARGUMENTS
!!
!!  idensvar -- the index of the density variable
!!
!!***
#ifdef DEBUG_ALL
#define DEBUG_GRAVITY
#endif


subroutine gr_isoFindMassCenter (idensvar)

!======================================================================

  use Driver_interface, ONLY : Driver_abortFlash
  use gr_isoMpoleData, ONLY : Xcm, Ycm,Zcm,mpole_geometry,Mtot
  use Grid_interface, ONLY : Grid_getListOfBlocks
  use Grid_data, ONLY : gr_meshComm

  implicit none
  
#include "constants.h"
#include "Flash.h"

  include "Flash_mpi.h"
  
  integer, intent(IN)  :: idensvar
  integer              :: blkCount
  integer, dimension(MAXBLOCKS) :: blkList
  
  integer, parameter :: nsum = 4    ! Number of summed quantities;
  real    :: sum(nsum), lsum(nsum)
  
  integer :: i,error
  !=================================================================
  call Grid_getListOfBlocks(LEAF, blklist, blkCount)
  !               Sum quantities over all locally held leaf blocks.
  if(.not.((mpole_geometry==SPHERICAL).and.(NDIM==2))) then
     sum(:)  = 0.
     lsum(:) = 0.
     do i = 1, blkCount
        call gr_isoSumLocal (lsum, nsum, blkList(i), idensvar)
     enddo
     
     !============================================================
     
     !               Give all processors a copy of the global sums.
     
     call mpi_allreduce (lsum, sum, nsum, FLASH_REAL, & 
          MPI_SUM, gr_meshComm, error)
     
     !               Now normalize the center-of-mass coordinates.
     
     Mtot = sum(1)
!     if (abs(Mtot) < 1.E-99) &  
!CD (11/2/07): Rejected by Absoft compiler as too small.  Strange 
!because the N117 compiler flag should promote reals to doubles.
     if (abs(Mtot) < tiny(1.)) &
          call Driver_abortFlash ("FATAL:  find_center_of_mass:  Mtot = 0")
     Xcm = sum(2) / Mtot
     Ycm = sum(3) / Mtot
     Zcm = sum(4) / Mtot
!  write (*,*) Xcm, Ycm, Zcm
  end if
  !=================================================================

  return
end subroutine gr_isoFindMassCenter

