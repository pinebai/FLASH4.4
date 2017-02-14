!!****if* source/Driver/DriverMain/Driver_getNStep
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
!!  Accessor function that returns the current step number in the
!! simulation from the Driver. 
!!
!! ARGUMENTS
!!  nstep - returned value, current step number 
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***

subroutine Driver_getNStep(nstep)

  use Driver_data, ONLY : dr_nStep

implicit none
  integer, intent(out) :: nstep
  
  nstep = dr_nStep

end subroutine Driver_getNStep

