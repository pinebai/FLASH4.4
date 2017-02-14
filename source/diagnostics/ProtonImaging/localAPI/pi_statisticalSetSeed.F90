!!****if* source/diagnostics/ProtonImaging/localAPI/pi_statisticalSetSeed
!!
!! NAME
!!
!!  pi_statisticalSetSeed
!!
!! SYNOPSIS
!!
!!  call pi_statisticalSetSeed (integer, intent (in) :: seed)
!!
!! DESCRIPTION
!!
!!  Initializes the random number generator with the supplied seed. After calling this
!!  routine, the random number generator is ready for use.
!!
!! ARGUMENTS
!!
!!  seed : the seed value for the random number generator
!!
!!***

subroutine pi_statisticalSetSeed (seed)

  implicit none

  integer, intent (in) :: seed

  return
end subroutine pi_statisticalSetSeed
