!!****f* source/Grid/Grid_setInterpValsGcell
!!
!!
!! NAME
!!
!!  Grid_setInterpValsGcell
!!
!!
!! SYNOPSIS
!!
!!  call Grid_setInterpValsGcell(logical(IN) :: setval)
!!
!!
!! DESCRIPTION
!!
!! Sets interpolation values for guardcell filling within INS fractional
!! step method.
!!
!!
!! ARGUMENTS
!!
!! setval = .true.  set values 
!!          .false. restore all interpolation-restriction values to quadratic
!!
!! NOTES
!!
!!  The current functionality is specifically for use by the IncompNS unit.
!!
!!  Only implemented for PARAMESH 4.
!!
!! DEV: This should be generalized, to be useful for other code than IncompNS!
!!***

subroutine Grid_setInterpValsGcell(setval)

  implicit none

  logical, intent(IN) :: setval

  return

end subroutine Grid_setInterpValsGcell

