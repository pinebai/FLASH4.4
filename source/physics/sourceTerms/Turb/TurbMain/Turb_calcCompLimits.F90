!!****if* source/physics/sourceTerms/Turb/TurbMain/Turb_calcCompLimits
!!
!! NAME
!!
!!  Turb_calcCompLimits
!!
!! SYNOPSIS
!!
!!  call Turb_calcCompLimits(integer, intent(IN), dimension(LOW:HIGH,MDIM)  :: blklimits,
!!                           integer, intent(OUT), dimension(LOW:HIGH,MDIM)  :: complimits,
!!                           integer, intent(IN)  :: local_stepsize)
!!
!! DESCRIPTION
!!
!! Aaron Jackson 2010
!! Calculate the indices over which the first operator should be 
!! performed. If turb_stepSize is small enough compared to the number
!! of guard cells and the type of operator being performed, then an
!! intermediate guard cell fill is not necessary and we should
!! calculate the first operator into the guard cells
!!
!! ARGUMENTS
!!
!!   blklimits : array holding upper and lower index limits of interior block cells (no GC)
!!
!!   complimits : array holding upper and lower index limits of block cells where the computation takes place
!!
!!   local_stepsize : local step size
!!
!!
!!
!!***

#include "constants.h"
#include "Flash.h"
subroutine Turb_calcCompLimits(blkLimits, compLimits, local_stepSize)

  use Turb_data, only : turb_stepSize, turb_fillGC

  implicit none
  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: blkLimits
  integer, intent(OUT), dimension(LOW:HIGH,MDIM) :: compLimits
  integer, intent(IN) :: local_stepSize

  ! check if we are not at flame refinement level
  ! or if we are performing guard cell fill between operators
  if ( turb_fillGC .or. (local_stepSize .ne. turb_stepSize) ) then

     ! we only need to compute for interior cells
     compLimits(:,:) = blkLimits(:,:)

  else

     ! we both operators can be performed without guard cell filling
     ! need to compute first operator into guard cells for second
     ! operator
     compLimits(LOW,IAXIS) = blkLimits(LOW,IAXIS)-2*turb_stepSize
     compLimits(HIGH,IAXIS) = blkLimits(HIGH,IAXIS)+2*turb_stepSize
     compLimits(LOW,JAXIS) = blkLimits(LOW,JAXIS)-K2D*2*turb_stepSize
     compLimits(HIGH,JAXIS) = blkLimits(HIGH,JAXIS)+K2D*2*turb_stepSize
     compLimits(LOW,KAXIS) = blkLimits(LOW,KAXIS)-K3D*2*turb_stepSize
     compLimits(HIGH,KAXIS) = blkLimits(HIGH,KAXIS)+K3D*2*turb_stepSize

  endif

  return
end subroutine Turb_calcCompLimits
