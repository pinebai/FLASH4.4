!!****if* source/Particles/ParticlesMain/passive/EstiMidpoint2/pt_prepareEsti
!!
!! NAME
!!
!!  pt_prepareEsti
!!
!! SYNOPSIS
!!
!!  call pt_prepareEsti(real(in)   :: dtOld,
!!                         real(in)   :: dtNew,
!!                         real(inout):: particles(:,p_count),
!!                         integer(in):: p_count,
!!                         integer(in):: ind)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the passive particle module.
!!  This portion cleans up and prepares for the next time step.
!!  
!!  This version implements an "estimated midpoint" advancement for
!!  passive particles. The scheme is formally second-order if the timestep
!!  remains unchanged.
!!  SCHEMATICALLY:
!!      x(t+1) = x(t) + dt * (1/2)*( vp(t+1) + vp(t))
!!   where vp is velocity evaluated at an estimated midpoint location
!!   reached between timesteps.
!!
!!  If time step has changed significantly (e.g. dtOld != dtNew, more
!!  or less, see pt_dtChangeTolerance under PARAMETERS below), then
!!  the previously estimated midpoint location is regarded as unusable
!!  and an Euler step is taken.
!!
!!  A new estimated midpoint is always calculated at the end.
!!
!! Actually in detail, normally:
!!      x(t+dtNew)      = x(t)       + (1/2)dtNew* ( vp(t+1/2*dtNew) + v(xp(t+1/2*dtNew),t+dtNew) )
!!      v(t+dtNew)      =                                              v(x(t+dtNew),t+dtNew)
!!  or, when an Euler step is taken instead:
!!        x(t+dtNew)    = x(t)      + dtNew *  v(x(t),      t)
!!      v(t+dtNew)      =                      v(x(t+dtNew),t+dtNew)
!! Followed in either case by:
!!      xp(t+1.5*dtNew) = x(t+dtNew) + (1/2)dtNew * v(t+dtNew) 
!!      vp(t+1.5*dtNew) =                           v(xp(t+1.5*dtNew),t+dtNew)
!!
!!
!! Problem:
!!  This is only a second-order method when dtNew == dtOld (exactly).
!!
!! ARGUMENTS
!!
!!   dtOld -- previous time interval
!!   dtNew -- current time interval
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind  -- index of this type into pt_infoType
!!  
!! PARAMETERS
!!
!!    pt_dtChangeTolerance REAL [0.4] Do Euler step if change in time step is greater than this
!!                                    percentage.  Set to 0 to always do Euler, set to a huge
!!                                    number to always use estimated midpoint velocities.
!!
!! SIDE EFFECTS
!!
!!  Updates the POS{X,Y,Z}, VEL{X,Y,Z}, POSPRED{X,Y,Z}, and VELPRED{X,Y,Z} properties of
!!  particles in the particles structure. On return, they will contain
!!   *  current positions at t+dtNew,
!!   *  current velocities evaluated (at t+dtNew) for current positions,
!!   *  positions of next estimated midpoints ("predicted" to be reached at t+1.5*dtNew),
!!      and
!!   *  velocities (evaluated at t+dtNew) for these next estimated midpoints,
!!  respectively.
!!
!!
!! NOTES
!!
!!   At the first call, it is assumed that particle initialization
!!   has filled in initial velocity components properly.
!!
!!   This is a slightly corrected version of the scheme that used to be called
!!   "PredictorCorrector" before the FLASH3.0 release.
!!
!!***

!===============================================================================
#ifdef DEBUG_ALL
#define DEBUG_PARTICLES
#endif

subroutine pt_prepareEsti (dtOld,dtNew,particles,p_count, ind)
  
  use Particles_data, ONLY: useParticles, pt_typeInfo,&
  pt_velNumAttrib,pt_posPredAttrib, pt_velPredAttrib
       
       use Grid_interface, ONLY : Grid_mapMeshToParticles
  
  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles
  
  integer :: mapType,i
  
  mapType=pt_typeInfo(PART_MAPMETHOD,ind)
  
  ! Compute the predicted midpoint positions.  These are staggered forward 1/2 step.
  !print *, 'about to compute estimated positions'
  do i = 1, p_count
     particles(POSPREDX_PART_PROP,i) =  particles(POSX_PART_PROP,i) + &
          0.5*dtNew * particles(VELX_PART_PROP,i)
     particles(POSPREDY_PART_PROP,i) =  particles(POSY_PART_PROP,i) + &
          0.5*dtNew * particles(VELY_PART_PROP,i)
     particles(POSPREDZ_PART_PROP,i) =  particles(POSZ_PART_PROP,i) + &
          0.5*dtNew * particles(VELZ_PART_PROP,i)
  enddo
  


  ! Map the current gas velocity field onto the estimated positions to get the
  ! predicted velocities.
  
  call Grid_mapMeshToParticles(particles,NPART_PROPS, BLK_PART_PROP,p_count,&
       pt_posPredAttrib,pt_velNumAttrib,pt_velPredAttrib,mapType)
  
  
  return
  
end subroutine pt_prepareEsti

