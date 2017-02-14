!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhCheckAccuracy
!!
!! NAME
!!
!!  Gravity_bhCheckAccuracy
!!
!!
!! SYNOPSIS
!!
!!  logical res = gr_bhCheckAccuracy(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Evaluates whether gravitational acceleration/potential within 
!!  a given block calculated by tree solver is accurate enough and 
!!  does not need to be recalculated.
!!
!! ARGUMENTS
!!
!!  blockno     : number of block into which the target cell belongs
!!  point       : indeces of the target cell in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!! RESULT
!!
!!  Returns TRUE if the gravitational potential/aceleration within 
!!  a given block calculated by tree solver is accurate enough and 
!!  does not need to be recalculated.
!!
!!***

logical function Gravity_bhCheckAccuracy(blockno, point, blkLimits, solnData)
  use Grid_interface, ONLY : Grid_getBlkPhysicalSize
  use Gravity_data, ONLY : grv_bhUseRelAccErr, grv_bhAccErr
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  integer :: ii, jj, kk
  real :: delx_inv, dely_inv, delz_inv
  real :: blockSize(MDIM), grav(MDIM), grol(MDIM)
  real :: err_lim2, gerr2

#ifdef GRAV_TREE_ACC
  grav = solnData(GACX_VAR:GACZ_VAR, point(IAXIS), point(JAXIS), point(KAXIS))  
  grol = solnData(GAOX_VAR:GAOZ_VAR, point(IAXIS), point(JAXIS), point(KAXIS))  
#else
  call Grid_getBlkPhysicalSize(blockno, blockSize)
  delx_inv = (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) / blockSize(IAXIS)
  dely_inv = (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) / blockSize(JAXIS)
  delz_inv = (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1) / blockSize(KAXIS)

  ii = point(IAXIS)
  jj = point(JAXIS)
  kk = point(KAXIS)

  ! only ii,jj,kk point was updated, therefore derivatives are calculated using
  ! the old potential at +1 points
  grav(IAXIS) = (solnData(GPOT_VAR,ii,jj,kk) - solnData(GPOL_VAR,ii+1,jj,kk))*delx_inv
  grav(JAXIS) = (solnData(GPOT_VAR,ii,jj,kk) - solnData(GPOL_VAR,ii,jj+1,kk))*dely_inv
  grav(KAXIS) = (solnData(GPOT_VAR,ii,jj,kk) - solnData(GPOL_VAR,ii,jj,kk+1))*delz_inv

  grol(IAXIS) = (solnData(GPOL_VAR,ii,jj,kk) - solnData(GPOL_VAR,ii+1,jj,kk))*delx_inv
  grol(JAXIS) = (solnData(GPOL_VAR,ii,jj,kk) - solnData(GPOL_VAR,ii,jj+1,kk))*dely_inv
  grol(KAXIS) = (solnData(GPOL_VAR,ii,jj,kk) - solnData(GPOL_VAR,ii,jj,kk+1))*delz_inv
#endif

  gerr2 = (grav(IAXIS) - grol(IAXIS))**2 &
  &     + (grav(JAXIS) - grol(JAXIS))**2 &
  &     + (grav(KAXIS) - grol(KAXIS))**2

  if (grv_bhUseRelAccErr) then
    err_lim2 = grv_bhAccErr**2 * (grav(IAXIS)**2+grav(JAXIS)**2+grav(KAXIS)**2)
  else
    err_lim2 = grv_bhAccErr**2
  endif

  !print *, "ChA: ", err_lim2, gerr2, grav, grol, grv_bhAccErr
    
  if (gerr2 .gt. err_lim2) then
    Gravity_bhCheckAccuracy = .false.            
  else
    Gravity_bhCheckAccuracy = .true.
  endif

  !print *, "ChA return: ", Gravity_bhCheckAccuracy 

  return
end function Gravity_bhCheckAccuracy


