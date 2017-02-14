!!****if* source/Particles/ParticlesMain/passive/Euler/pt_advanceEuler_passive
!!
!! NAME
!!
!!  pt_advanceEuler_passive
!!
!! SYNOPSIS
!!
!!  call pt_advanceEuler_passive(real(in)   :: dtOld,
!!                         real(in)   :: dtNew,
!!                         real(inout):: particles(:,p_count),
!!                         integer(in):: p_count,
!!                         integer(in):: ind)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  
!!  This version implements 1st order Euler integration in time.
!!  SCHEMATICALLY:
!!      x(t+1) = x(t) + dtNew * v(t)
!!
!! In detail:
!!      x(t+dtNew)  = x(t)      + dtNew * v(x(t),      t)
!!      v(t+dtNew)  =                     v(x(t+dtNew),t+dtNew)
!!
!!
!! ARGUMENTS
!!
!!   dtOld -- not used in this first-order scheme
!!   dtNew -- current time increment
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind   -- index for type into pt_typeInfo data structure
!!  
!!
!! SIDE EFFECTS
!!
!!  Updates the POS{X,Y,Z} and VEL{X,Y,Z} properties of particles in the particles structure.
!!  Sorts particles in the particles structure by calling Grid_sortParticles (indirectly).
!!
!! NOTES
!!
!!  No special handling is done for the first call - it is assumed that particle
!!  initialization fills in initial velocity components properly.
!!***

!===============================================================================

subroutine pt_advanceEuler_passive (dtOld,dtNew,particles,p_count, ind)
    
  use Particles_data, ONLY: useParticles, pt_typeInfo,&
       pt_posAttrib, pt_velNumAttrib,pt_velAttrib



  use Grid_interface, ONLY : Grid_mapMeshToParticles
  
  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Particles.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles
  integer       :: i,nstep
  real          :: jumpx,jumpy,jumpz
  integer :: part_props = NPART_PROPS

!!------------------------------------------------------------------------------
  
  integer :: mapType 
!!------------------------------------------------------------------------------
  
  ! Don't do anything if runtime parameter isn't set
  if (.not.useParticles ) return

  mapType=pt_typeInfo(PART_MAPMETHOD,ind)


     
  ! Update the particle positions.
  do i = 1, p_count
 
     jumpx = dtNew * particles(VELX_PART_PROP,i)
     particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + jumpx

     jumpy = dtNew * particles(VELY_PART_PROP,i)
     particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + jumpy

     jumpz = dtNew * particles(VELZ_PART_PROP,i)
     particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + jumpz

  enddo
     
  ! Map the updated gas velocity field onto the current particle positions to
  ! obtain the updated particle velocities.

  call Grid_mapMeshToParticles(particles,&
       part_props, BLK_PART_PROP,p_count,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)
  
  
  return
!!------------------------------------------------------------------------------
  
end subroutine pt_advanceEuler_passive


