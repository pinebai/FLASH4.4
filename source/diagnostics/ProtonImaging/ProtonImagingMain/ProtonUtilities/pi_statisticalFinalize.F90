!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonUtilities/pi_statisticalFinalize
!!
!! NAME
!!
!!  pi_statisticalFinalize
!!
!! SYNOPSIS
!!
!!  call pi_statisticalFinalize ()
!!
!! DESCRIPTION
!!
!!  Finalizes the statistical environment that was used by the proton imaging code. It currently
!!  consists in deallocation of the seed array that was needed for the random number generator.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine pi_statisticalFinalize ()

  use ProtonImaging_data, ONLY : pi_randomNumberSeedArray

  implicit none
!
!
!     ...Deallocate the seed array.
!
!
  if (allocated (pi_randomNumberSeedArray)) deallocate (pi_randomNumberSeedArray)
!
!
!     ...Ready!
!
!
  return
end subroutine pi_statisticalFinalize
