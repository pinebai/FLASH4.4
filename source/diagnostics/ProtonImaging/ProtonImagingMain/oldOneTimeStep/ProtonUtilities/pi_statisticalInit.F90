!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonUtilities/pi_statisticalInit
!!
!! NAME
!!
!!  pi_statisticalInit
!!
!! SYNOPSIS
!!
!!  call pi_statisticalInit ()
!!
!! DESCRIPTION
!!
!!  Initializes the statistical environment to be used by the proton imaging code. It currently
!!  consists in allocation of the seed array needed for the random number generator.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine pi_statisticalInit ()

  use ProtonImaging_data, ONLY : pi_randomNumberSeedArray

  use ut_randomInterface, ONLY : ut_randomSeed

  implicit none

  integer :: seedArraySize
!
!
!     ...Allocate the needed seed array.
!
!
  call ut_randomSeed         (ut_size = seedArraySize)      ! get compiler dependent seed array size
  allocate (pi_randomNumberSeedArray (1:seedArraySize))     ! allocate properly sized seed array
!
!
!     ...Ready!
!
!
  return
end subroutine pi_statisticalInit
