!!****if* source/physics/materialProperties/Opacity/localAPI/op_setPEarrayJmax
!!
!! NAME
!!
!!  op_setPEarrayJmax
!!
!! SYNOPSIS
!!
!!  call op_setPEarrayJmax ()
!!
!! DESCRIPTION
!!
!!  This routine sets the Jmax array for determining the photoelectric
!!  cross sections according to the F.Biggs and R.Lighthill report:
!!
!!       Analytical Approximations for X-Ray Cross Sections II
!!       Frank Biggs and Ruth Lighthill
!!       Weapons Effects Research Department
!!       Sandia Laboratories, December 1971
!!
!!  The data is NOT the updated data from the 1988 update of the report!
!!  This routine can only be called after all the necessary arrays have
!!  been allocated.
!!
!! ARGUMENTS
!!
!!***
subroutine op_setPEarrayJmax ()

  implicit none

  return
end subroutine op_setPEarrayJmax
