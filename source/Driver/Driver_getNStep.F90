!!****f* source/Driver/Driver_getNStep
!!
!! NAME
!!  Driver_getNStep
!!
!! SYNOPSIS
!!  
!!
!!  Driver_getNStep(real(out) :: nstep)
!!  
!! DESCRIPTION 
!!
!!  Accessor funtion that returns the current step number in the
!! simulation
!!
!! ARGUMENTS
!!  nstep - returned value, current step number 
!!
!!
!!***


!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in FORTRAN
!! module Driver_data (in file Driver_data.F90). The other variables
!! are local to the specific routine and do not have the prefix "dr_"

subroutine Driver_getNStep(nstep)

implicit none
  integer, intent(out) :: nstep

  !dummy value for stub
  nstep = 0

end subroutine Driver_getNStep

