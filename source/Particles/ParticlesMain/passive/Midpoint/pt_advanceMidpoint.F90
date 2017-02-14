!!****if* source/Particles/ParticlesMain/passive/Midpoint/pt_advanceMidpoint
!!
!! NAME
!!
!!  pt_advanceMidpoint
!!
!! SYNOPSIS
!!
!!  pt_advanceMidpoint(real(in) :: dtOld,
!!                    real(in) :: dtNew,
!!                    real(inout):: particles(:,p_count),
!!                    integer(in):: p_count,
!!                    integer(in):: ind)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  
!!  This version implements 2nd-order midpoint advancement for
!!  passive particles. SCHEMATICALLY:
!!      x(t) = x(t-2) + 2*dt*v(t-1)
!!  or, for cases where the time step can change
!!      x(t) = x(t-2) + (dt(t-1)+dt(t-2))*v(t-1)
!!
!!
!! In detail:
!! (a) Normal case:
!!     x(t+dtNew) = x(t-dtOld) + (dtOld+dtNew)* v(x(t), t)
!!
!! (b) Initial step (first call): "Euler step"
!!     x(t+dtNew) = x(t)       +        dtNew * v(x(t), t)
!!
!! The velocity components are set as
!!      v(t+dtNew)  =                     v(x(t+dtNew),t+dtNew)
!!
!! ARGUMENTS
!!
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind -- index of this particle type in pt_typeInfo
!!
!! SIDE EFFECTS
!!
!!  Updates the POS{X,Y,Z}, POS2PREV{X,Y,Z}, and VEL{X,Y,Z} properties of particles in the particles structure.
!!  Sorts particles in the particles structure by calling Grid_sortParticles (indirectly).
!!
!! NOTES
!!
!!  At the first call, it is assumed that particle initialization
!!  has filled in initial velocity components properly.
!!  
!!***

!===============================================================================

subroutine pt_advanceMidpoint (dtOld,dtNew,particles,p_count, ind)
  
  use Particles_data, ONLY: useParticles, pt_typeInfo,&
       pt_velNumAttrib,pt_posAttrib, pt_velAttrib
       
  use Grid_interface, ONLY : Grid_mapMeshToParticles

  
  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles
  
  integer :: mapType

  integer       :: i,nstep
  real          :: posXminus1, posYminus1, posZminus1
  logical, save :: first_call = .true.
  integer       :: part_props=NPART_PROPS
!!--------------------------------------------------------------------------------------
  
  ! Don't do anything if runtime parameter isn't set
  if (.not.useParticles ) return
  mapType=pt_typeInfo(PART_MAPMETHOD,ind)
  
  !! If this is the first sweep through, calculate the new positions by an Euler step.
  !! Then save them for the next jump step
  if (first_call) then
     
     ! Map the updated (xyz sweep) velocity field onto the current particle positions to
     ! obtain the particle velocities at the midpoints
     
     
     ! Save the initial positions for the second time step and update the particle positions
     do i = 1, p_count
        particles(POS2PREVX_PART_PROP,i) = particles(POSX_PART_PROP,i)
        particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + &
             dtNew * particles(VELX_PART_PROP,i)
        if (NDIM .GT. 1) then
           particles(POS2PREVY_PART_PROP,i) = particles(POSY_PART_PROP,i)
           particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + &
                dtNew * particles(VELY_PART_PROP,i)
        endif
        if (NDIM .EQ. 3) then
           particles(POS2PREVZ_PART_PROP,i) = particles(POSZ_PART_PROP,i)
           particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + &
                dtNew * particles(VELZ_PART_PROP,i)
        endif
     enddo  !! end of loop over all particles
     
     
     first_call = .false.
     
     !! end of initial call and Euler step
     
  else  ! (.not. first_call)
     
     !! Subsequent sweeps, calculate the final value from the midpoint velocity and two previous location
     !! Also save positions for the next step
     
     
     !! Update positions using the locations from two steps ago and the velocity
     !!    from one step ago.
     !! Complicated by needing to save the current particle position for future use
     !!    in a temporary variable
     do i = 1, p_count
        posXminus1 = particles(POSX_PART_PROP,i)
        particles(POSX_PART_PROP,i) = particles(POS2PREVX_PART_PROP,i) + & 
             (dtOld + dtNew) * particles(VELX_PART_PROP,i)
        particles(POS2PREVX_PART_PROP,i) = posXminus1
        
        if (NDIM .GT. 1) then
           posYminus1 = particles(POSY_PART_PROP,i)
           particles(POSY_PART_PROP,i) = particles(POS2PREVY_PART_PROP,i) + &
                (dtOld + dtNew) * particles(VELY_PART_PROP,i)
           particles(POS2PREVY_PART_PROP,i) = posYminus1
           
        endif
        if (NDIM .EQ. 3) then
           posZminus1 = particles(POSZ_PART_PROP,i)
           particles(POSZ_PART_PROP,i) = particles(POS2PREVZ_PART_PROP,i) + &
                (dtOld + dtNew) * particles(VELZ_PART_PROP,i)
           particles(POS2PREVZ_PART_PROP,i) = posZminus1
        endif
     enddo
     
     
  endif  ! end of (.not. first_call)
  
  
  call Grid_mapMeshToParticles(particles,&
       part_props, BLK_PART_PROP,p_count,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)
  
  !!--------------------------------------------------------------------------------------
  
  return
  
end subroutine pt_advanceMidpoint

