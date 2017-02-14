!!****if* source/Simulation/SimulationMain/Orbit/pt_initPositions
!!
!! NAME
!!
!!  pt_initPositions
!!
!! SYNOPSIS
!!
!!  pt_initPositions(integer, INTENT(in)  :: blockid)
!!
!! DESCRIPTION
!!    An override for the Orbit simulation.  Sets up only two particles; they should orbit!
!!
!!
!! ARGUMENTS
!!
!!   blockid : ID of block in current processor
!!
!!
!!***

subroutine pt_initPositions (blockID, success)

  use Simulation_data, ONLY: sim_xmin, sim_xmax, sim_ymin, sim_ymax, sim_zmin, sim_zmax, &
       sim_extfield, sim_Newton, sim_ptMass, sim_Sep, sim_nPtot
  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, &
    Grid_getBlkCenterCoords
  use Particles_data, ONLY: pt_numLocal, particles, pt_maxPerProc,pt_meshMe

#include "Flash.h"
#include "constants.h"

  implicit none

  integer, INTENT(in) :: blockID
  logical, intent(out) :: success
  integer       :: i, p
  logical       :: IsInBlock
  real          :: xpos, ypos, zpos, bxLower, byLower, bzLower, bxUpper, byUpper, bzUpper
  real          :: xvel, yvel, zvel, blockSize(3), blockCenter(3)

!-------------------------------------------------------------------------------

! Particle slot number (incremented and saved between calls)

  p = pt_numLocal

!-------------------------------------------------------------------------------


! Get locations of block faces.

  call Grid_getBlkPhysicalSize(blockID, blockSize)
  call Grid_getBlkCenterCoords(blockID, blockCenter)
  bxLower    = blockCenter(1) - 0.5*blockSize(1)
  bxUpper    = blockCenter(1) + 0.5*blockSize(1)
  if (NDIM >= 2) then
     byLower = blockCenter(2) - 0.5*blockSize(2)
     byUpper = blockCenter(2) + 0.5*blockSize(2)
  endif
  if (NDIM == 3) then
     bzLower = blockCenter(3) - 0.5*blockSize(3)
     bzUpper = blockCenter(3) + 0.5*blockSize(3)
  endif

! Loop over both particles and compute their positions.

  ypos = 0.5*(sim_yMin+sim_yMax) * K2D
  zpos = 0.5*(sim_zMin+sim_zMax) * K3D

  xvel = 0.
  if (sim_nPtot == 1) then                 ! only one particle
     yvel = 0.
  else if (sim_extField) then              ! fixed central potential
     yvel = sqrt(sim_Newton*sim_ptMass/(0.5*sim_Sep))
  else                                     ! self-gravitating particles
     yvel = 0.5*sqrt(sim_Newton/(0.5*sim_Sep))
  endif
  zvel = 0.

  do i = 1, sim_nPtot       ! 1 or 2 particles in orbit
     xpos = 0.5*(sim_xMin+sim_xMax) + 0.5*sim_Sep*(2*i-3)

! Check to see if the particle lies within this block.

     IsInBlock = (xpos >= bxLower) .and. (xpos < bxUpper)
     if (NDIM >= 2) &
          IsInBlock = IsInBlock .and. ((ypos >= byLower) .and. (ypos < byUpper))
     if (NDIM == 3) &
          IsInBlock = IsInBlock .and. ((zpos >= bzLower) .and. (zpos < bzUpper))

! If yes, and adequate particle buffer space is available, initialize it.

     if (IsInBlock) then
        p = p + 1
        if (p > pt_maxPerProc) then
           call Driver_abortFlash &
                ("InitParticlePositions:  Exceeded max # of particles/processor!")
        endif

! Particle current block number.
        particles(BLK_PART_PROP,p) = real(blockID)
        particles(PROC_PART_PROP,p) = real(pt_meshMe)
! Particle mass.
#ifdef MASS_PART_PROP
        if (MASS_PART_PROP > 0) particles(MASS_PART_PROP,p) = 1.
#endif 
! Particle position and velocity.
        particles(POSX_PART_PROP,p) = xpos
        particles(VELX_PART_PROP,p) = xvel
        particles(POSY_PART_PROP,p) = ypos
        particles(VELY_PART_PROP,p) = (2*i-3)*yvel
        particles(POSZ_PART_PROP,p) = zpos
        particles(VELZ_PART_PROP,p) = zvel

     endif
  enddo
    
! Set the particle database local number of particles.

  pt_numLocal = p

  success = .true.
  return

!-------------------------------------------------------------------------------
  
end subroutine pt_initPositions
