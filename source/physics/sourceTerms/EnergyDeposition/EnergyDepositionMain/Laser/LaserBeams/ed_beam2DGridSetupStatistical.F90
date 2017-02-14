!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam2DGridSetupStatistical
!!
!! NAME
!!
!!  ed_beam2DGridSetupStatistical
!!
!! SYNOPSIS
!!
!!  call ed_beam2DGridSetupStatistical (logical, intent (in)    :: seedInitialize,
!!                                      logical, intent (in)    :: seedIncrement,
!!                                      integer, intent (inout) :: seedMaximum,
!!                                      integer, intent (inout) :: seedStepping,
!!                                      integer, intent (out)   :: seed)
!!
!! DESCRIPTION
!!
!!  Sets up information needed for establishing a statistical linear cross sectional grid for
!!  a particular planar 2D beam. Currently this consists only in setting a seed value for
!!  a random uniform number generator supplied by Fortran 90 and the values for the maximum
!!  seed and the seed stepping number.
!!
!! ARGUMENTS
!!
!!  seedInitialize : indicates, whether the initial seed values need to be set
!!  seedIncrement  : indicates, whether the seed needs to be incremented
!!  seedMaximum    : the maximum seed value allowed
!!  seedStepping   : the seed stepping value
!!  seed           : the seed value
!!
!! NOTES
!!
!!  Only one seed value is needed to generate a definite (reproducible) sequence of random
!!  numbers using the canned f90 random number generator.
!!
!!***

subroutine ed_beam2DGridSetupStatistical (seedInitialize,               &
                                          seedIncrement,                &
                                                          seedMaximum,  &
                                                          seedStepping, &
                                                          seed          )

  use EnergyDeposition_data,  ONLY : ed_largestPositiveInteger

  implicit none

  logical, intent (in)    :: seedInitialize
  logical, intent (in)    :: seedIncrement
  integer, intent (inout) :: seedMaximum
  integer, intent (inout) :: seedStepping
  integer, intent (out)   :: seed
!
!
!     ...If requested, set the initial seed, the seed stepping value and the maximum
!        seed value.
!
!
  if (seedInitialize) then
      seed         = 0
      seedStepping = 1
      seedMaximum  = ed_largestPositiveInteger
  end if
!
!
!     ...If requested, increment the seed by the provided stepping value. If the incremented
!        seed value exceeds the maximum seed value, reset the seed to zero.
!
!
  if (seedIncrement) then
      if (seed > seedMaximum - seedStepping) then
          seed = 0
      else
          seed = seed + seedStepping
      end if
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beam2DGridSetupStatistical
