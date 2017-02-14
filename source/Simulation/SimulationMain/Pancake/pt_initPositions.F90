!!****if* source/Simulation/SimulationMain/Pancake/pt_initPositions
!!
!! NAME
!!
!!   pt_initPositions
!!
!! SYNOPSIS
!!
!!   pt_initPositions(integer, INTENT(in) :: blockId, 
!!                    logical, INTENT(out) :: success)
!!
!! DESCRIPTION
!!   An override for the Zel'dovich pancake problem.
!!
!!
!! ARGUMENTS
!!
!!   blockId : Id of block in current processor
!!   success: logical argument indicating whether particles initalised 
!!            correctly on this block.
!!
!!
!!***

subroutine pt_initPositions(blockId, success)

use Simulation_data, ONLY : sim_nxp, sim_nyp, sim_nzp, sim_xmin, sim_ymin, sim_zmin, &
     sim_xcos, sim_ycos, sim_zcos, sim_xmax, sim_ymax, sim_zmax, sim_hubble, &
     sim_zcaustic, sim_zinitial, sim_kmag, sim_temp2, sim_mr
use Particles_data, ONLY : pt_numLocal, pt_meshMe
use pt_interface, ONLY : pt_createTag
use Grid_interface, ONLY : Grid_getBlkPhysicalSize, &
                           Grid_getBlkCenterCoords

#include "Flash.h"
#include "constants.h"

implicit none

integer, intent(in) :: blockId
logical, intent(out) :: success

integer       :: i, j, k, tr
real          :: xr, yr, zr, vxr, vyr, vzr
real          :: xe, xl, xl_unpert, yl_unpert, zl_unpert, v
integer, save :: p

logical       :: IsInBlock
real          :: bxLower, bxUpper, byLower, byUpper, bzLower, bzUpper
real, dimension(3) :: blockSize, blockCenter

! Particle slot number (incremented and saved between calls)
p = pt_numLocal

! Location of block faces
call Grid_getBlkPhysicalSize(blockId, blockSize)
call Grid_getBlkCenterCoords(blockId, blockCenter)
bxLower = blockCenter(1) - 0.5*blockSize(1)
bxUpper = blockCenter(1) + 0.5*blockSize(1)
if (NDIM >= 2) then
   byLower = blockCenter(2) - 0.5*blockSize(2)
   byUpper = blockCenter(2) + 0.5*blockSize(2)
endif
if (NDIM == 3) then 
   bzLower = blockCenter(3) - 0.5*blockSize(3)
   bzUpper = blockCenter(3) + 0.5*blockSize(3)
endif

! Compute perturbed particle positions

if (pt_meshMe == MASTER_PE) &
  print *, 'Particles_initPositions:  computing perturbed particle positions...'

do k = 1, sim_nzp
  zl_unpert = sim_zmin + (k-0.5)*(sim_zmax-sim_zmin)/sim_nzp
  do j = 1, sim_nyp
    yl_unpert = sim_ymin + (j-0.5)*(sim_ymax-sim_ymin)/sim_nyp
    do i = 1, sim_nxp
      xl_unpert = sim_xmin + (i-0.5)*(sim_xmax-sim_xmin)/sim_nxp

! Unperturbed position is "Lagrangian" coordinate; compute perturbed "Eulerian"
! coordinate

      xl = (xl_unpert - 0.5*(sim_xmin+sim_xmax)) * sim_xcos + &
           (yl_unpert - 0.5*(sim_ymin+sim_ymax)) * sim_ycos + &
           (zl_unpert - 0.5*(sim_zmin+sim_zmax)) * sim_zcos
      xe = xl - sim_temp2 * sin(sim_kmag*xl)

      xr = xl_unpert + (xe-xl)*sim_xcos
      yr = (yl_unpert + (xe-xl)*sim_ycos) * K2D
      zr = (zl_unpert + (xe-xl)*sim_zcos) * K3D

! Enforce periodic boundary conditions

      if (xr < sim_xmin) xr = sim_xmax + (xr-sim_xmin)
      if (xr > sim_xmax) xr = sim_xmin + (xr-sim_xmax)
      if (yr < sim_ymin) yr = sim_ymax + (yr-sim_ymin)
      if (yr > sim_ymax) yr = sim_ymin + (yr-sim_ymax)
      if (zr < sim_zmin) zr = sim_zmax + (zr-sim_zmin)
      if (zr > sim_zmax) zr = sim_zmin + (zr-sim_zmax)

! Compute velocity

      v = -sim_hubble * (1.0+sim_zcaustic) * sqrt(1.0+sim_zinitial) * & 
                        sin(sim_kmag*xl) / sim_kmag
      vxr = v * sim_xcos
      vyr = v * sim_ycos
      vzr = v * sim_zcos

! Check if particle is in this block 
      IsInBlock = (xr >= bxLower) .and. (xr < bxUpper)
      if (NDIM >= 2) &
         IsInBlock = IsInBlock .and. ((yr >= byLower) .and. (yr < byUpper))
      if (NDIM == 3) & 
         IsInBlock = IsInBlock .and. ((zr >= bzLower) .and. (zr < bzUpper))

! If it is, keep it; otherwise discard it.
      if (IsInBlock) then
        p = p + 1
        call InitSingleParticle(p, xr, yr, zr, vxr, vyr, vzr, sim_mr, blockId)
      endif

    enddo
  enddo
enddo

! Set the particle database local number of particles.
pt_numLocal = p
success = .true.

return
end subroutine pt_initPositions


!**********************************************************************
!  Routine:     InitSingleParticle

!  Description: Initialize a single particle with given characteristics.


subroutine InitSingleParticle (p, xpos, ypos, zpos, xvel, yvel, zvel, &
                               mass, block)

  use Particles_data, ONLY : pt_maxPerProc, particles,pt_meshMe
  use Driver_interface, ONLY : Driver_abortFlash

#include "Flash.h"
#include "constants.h"

  implicit none

  real, intent(in)    :: xpos, ypos, zpos, xvel, yvel, zvel, mass
  integer, intent(in) :: p, block

  if (p > pt_maxPerProc) &
    call Driver_abortFlash("InitSingleParticle:  Exceeded max # of particles!")

  particles(BLK_PART_PROP,p)  = real(block)
  particles(PROC_PART_PROP,p) = real(pt_meshMe)
  particles(MASS_PART_PROP,p) = mass
  particles(POSX_PART_PROP,p) = xpos
  particles(VELX_PART_PROP,p) = xvel
  particles(POSY_PART_PROP,p) = ypos
  particles(VELY_PART_PROP,p) = yvel
  particles(POSZ_PART_PROP,p) = zpos
  particles(VELZ_PART_PROP,p) = zvel

return
end subroutine InitSingleParticle
