!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonUtilities/pi_statisticalSetSeed
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

  use ProtonImaging_data, ONLY : pi_randomNumberSeedArray
  use ut_randomInterface, ONLY : ut_randomSeed

  implicit none

  integer, intent (in) :: seed
!
!
!     ...Set the seed for the random number generator.
!
!
  pi_randomNumberSeedArray (:) = seed                       ! set the seed array
  call ut_randomSeed   (ut_put = pi_randomNumberSeedArray)  ! commit seed array to random number generator
!
!
!     ...Ready!
!
!
  return
end subroutine pi_statisticalSetSeed
