!!****f* source/Particles/Particles_sinkAdvanceParticles
!!
!! NAME
!!
!!  Particles_sinkAdvanceParticles
!!
!! SYNOPSIS
!!
!!  call Particles_sinkAdvanceParticles(real, INTENT(in) :: dr_dt)
!!
!! DESCRIPTION
!!
!!  Updates sink particle postions and velocities based on gravitational accelerations,
!!  using Leapfrog or Euler, depending on the user's choice. There is also a special
!!  implementation of Leapfrog for cosmological simulations.
!!  This routine performs subcycling on sink-sink interactions, in case of very close
!!  encounters and highly eccentric orbits.
!!
!! ARGUMENTS
!!
!!   dr_dt - the current time step
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2008-2015
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   removed call to AccelGasOnSinks (now called in Gravity_PotentialListOfBlocks)
!!   (CTSS, CF 2013)
!!
!!***

subroutine Particles_sinkAdvanceParticles(dr_dt)
  implicit none
  real, INTENT(in) :: dr_dt
end subroutine Particles_sinkAdvanceParticles
