!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_laserIOCheckWrite
!!
!! NAME
!!
!!  ed_laserIOCheckWrite
!!
!! SYNOPSIS
!!
!!  call ed_laserIOCheckWrite (integer, optional, intent (in) :: passSplitDriver)
!!
!! DESCRIPTION
!!
!!  This routine is called each cycle to check to see if it is time to write laser ray
!!  information to the plot file. If so, the variable 'ed_laserIOWrite' is set to true.
!!
!!  It is time to write ray data if:
!!
!!        1. The run-time parameter ed_useLaserIO is true
!!        2. This is NOT the second half of the time step when split driver is in use
!!        3. This is the cycle immediately following a plot
!!
!! ARGUMENTS
!!
!!  passSplitDriver : indicates first/second half of time step for split driver
!!
!!***

subroutine ed_laserIOCheckWrite (passSplitDriver)

  implicit none

  integer, optional, intent (in) :: passSplitDriver

  return
end subroutine ed_laserIOCheckWrite
