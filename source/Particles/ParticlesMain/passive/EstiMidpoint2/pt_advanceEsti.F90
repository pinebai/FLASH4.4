!!****if* source/Particles/ParticlesMain/passive/EstiMidpoint2/pt_advanceEsti
!!
!! NAME
!!
!!  pt_advanceEsti
!!
!! SYNOPSIS
!!
!!  pt_advanceEsti(real(in) :: dtOld,
!!                    real(in) :: dtNew,
!!                    real(inout):: particles(:,p_count),
!!                    integer(in):: p_count,
!!                    integer(in):: ind)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the passive particle module.
!!  
!!  This version implements a "corrected estimated midpoint" advancement
!!  for passive particles. The correction should keep the scheme second-order
!!  even when the timestep changes.
!!
!!  If time step has changed significantly (e.g. dtOld very different
!!  from dtNew, see pt_dtChangeTolerance under PARAMETERS below), then
!!  the previously estimated midpoint location is ignored
!!  and an Euler step is taken.
!!
!!  A new estimated midpoint is always calculated at the end.
!!
!! Actually in detail, normally:
!!      x(t+dtNew)      = x(t)       + dtNew* vcombined ,
!!  where
!!      vcombined =  c1 * vp(t+1/2*dtNew)
!!                 + c2 * v(xp(t+1/2*dtNew),t+dtNew)
!!                 + c3 * v(t)
!!                 + c4 * v(x(t) ,          t+dtNew)
!!               ==  c1 * v(xp(t+1/2*dtNew),t)
!!                 + c2 * v(xp(t+1/2*dtNew),t+dtNew)
!!                 + c3 * v(x(t) ,          t)
!!                 + c4 * v(x(t) ,          t+dtNew)
!!               ==  c1 * v(x(t)+(1/2)dtOld*v(t),t)
!!                 + c2 * v(x(t)+(1/2)dtOld*v(t),t+dtNew)
!!                 + c3 * v(x(t) ,               t)
!!                 + c4 * v(x(t) ,               t+dtNew)
!!  -- with weights c1,c2,c3,c4 chosen appropriately, see below --
!!  or, when an Euler step is taken instead:
!!        x(t+dtNew)    = x(t)      + dtNew *  v(x(t),      t)
!! Followed in either case by:
!!      v(t+dtNew)      =                           v(x(t+dtNew),t+dtNew)
!!      xp(t+1.5*dtNew) = x(t+dtNew) + (1/2)dtNew * v(t+dtNew) 
!!      vp(t+1.5*dtNew) =                           v(xp(t+1.5*dtNew),t+dtNew)
!!
!!  Note that with c1=c2=0.5 and c3=c4=0, this corresponds to the (uncorrected)
!!  "estimated midpoint" scheme. By choosing appropriate values for these
!!  weights, the scheme is corrected so that it doesn't lose its second-orderness
!!  when dtNew!=dtOld.
!!
!!  It turns out that the constants have to satisfy three conditions:
!!   c1+c2+c3+c4 == 1  ,  (otherwise the scheme will not even be of first order)
!!      c2   +c4 == 0.5, 
!!   c1+c2       == dtNew / dtOld .
!!
!!  So there remains one parameter that can be chose freely.  The implementation
!!  code here assumes c4=0, but this could be changed if necessary.
!!
!! ARGUMENTS
!!
!!   dtOld -- previous time interval
!!   dtNew -- current time interval
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind -- index of this particle type into pt_typeInfo
!!  
!! PARAMETERS
!!
!!    pt_dtChangeToleranceUp REAL [5.0] Do Euler step if time step increases and relative change is
!!                                      greater than this.  Set to 0 to always do Euler, set to a huge
!!                                      number to always use estimated midpoint velocities.
!!    pt_dtChangeToleranceDown REAL [0.8] Do Euler step if time step decreases and relative change is
!!                                      greater than this.  Set to 0 to always do Euler, set to a large
!!                                      number to always use estimated midpoint velocities.
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
!!   This is a generalization of the scheme that used to be called
!!   "PredictorCorrector" before the FLASH3.0 release.
!!
!!   The pt_dtChangeTolerance runtime parameter can be set to a
!!   larger value here than for the original "EstiMidpoint" scheme,
!!   since this scheme corrects for changes in the time step.
!!***

!===============================================================================
#ifdef DEBUG_ALL
#define DEBUG_PARTICLES
#endif

subroutine pt_advanceEsti (dtOld,dtNew,particles,p_count, ind)
  
  use Particles_data, ONLY: useParticles, pt_typeInfo,&
       pt_velNumAttrib,pt_posPredAttrib, pt_velPredAttrib,&
       pt_posAttrib, pt_velAttrib,pt_dtChangeTolerance, pt_restart,&
       pt_meshMe
       
  use Grid_interface, ONLY : Grid_mapMeshToParticles
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  
  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles
  
  integer :: mapType
  

  integer       :: i
  logical, save :: firstCall = .true.  ! initialized on compilation
  logical       :: doEstiMidPoint, tolerableStepChange
  integer, parameter :: part_props=NPART_PROPS

  real, save ::  pt_dtChangeToleranceUp, pt_dtChangeToleranceDown
  real :: c1,c2,c3,c4, timeStepStretchFactor
!! ----------------------------------------------------------------------------

  ! Don't do anything if runtime parameter isn't set
  if (.not.useParticles ) return
  mapType=pt_typeInfo(PART_MAPMETHOD,ind)
  
  
  ! Initializations
  if (firstCall) then
     call RuntimeParameters_get ("pt_dtChangeToleranceUp", pt_dtChangeToleranceUp)
     call RuntimeParameters_get ("pt_dtChangeToleranceDown", pt_dtChangeToleranceDown)
  end if



  timeStepStretchFactor = 0.5 * dtOld / dtNew
  c4 = 0                        ! May be changed here
  c2 = 0.5                      ! May be changed here; keep c2 + c4 == 0.5 !
  c1 = 0.5/timeStepStretchFactor - c2
  c3 = 1 - c1 - c2 - c4


  tolerableStepChange = ( ((dtNew - dtOld) <= dtOld*pt_dtChangeToleranceUp) &
                    .AND. ((dtOld - dtNew) <= dtOld*pt_dtChangeToleranceDown) )
  doEstiMidPoint = tolerableStepChange .AND. ((.NOT.firstCall) .OR. pt_restart)
   
  ! For the first timestep, or perhaps if timesteps have changed very much, we do an Euler step.
  if (.NOT. doEstiMidPoint) then
     if (pt_meshMe .eq. MASTER_PE) print*,'pt_advanceEsti: Euler only at this time step!'
    ! Update the particle positions.
     do i = 1, p_count
        particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + &
             dtNew * particles(VELX_PART_PROP,i)
        particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + &
             dtNew * particles(VELY_PART_PROP,i)
        particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + &
             dtNew * particles(VELZ_PART_PROP,i)
     enddo
     

         ! end of know-nothing Euler step  
  else   ! time steps are similar enough, and previous values exist
     
     !-------------------------------------------------------------------------------
     ! For subsequent timesteps, we have predicted positions and velocities, so do
     ! a "corrected estimated midpoint" step.
#ifdef DEBUG_PARTICLES
99   format ('c[1:4]=(',3(1PG9.3,','),(1PG9.3),')')
  print 99,c1,c2,c3,c4
#endif
     do i = 1, p_count
        particles(VELPREDX_PART_PROP,i) = &
             (c3*particles(VELX_PART_PROP,i) + &
             c1*particles(VELPREDX_PART_PROP,i))
        particles(VELPREDY_PART_PROP,i) = &
             (c3*particles(VELY_PART_PROP,i) + &
             c1*particles(VELPREDY_PART_PROP,i))
        particles(VELPREDZ_PART_PROP,i) = &
             (c3*particles(VELZ_PART_PROP,i) + &
             c1*particles(VELPREDZ_PART_PROP,i))
     enddo
     ! Map the updated gas velocity field onto the predicted positions to obtain the
     !  predicted particle velocities.
     call Grid_mapMeshToParticles(particles,&
          part_props, BLK_PART_PROP,p_count,&
          pt_posPredAttrib,pt_velNumAttrib,pt_velAttrib,mapType)
     
     ! Update the particle positions using the predicted particle velocities ???
     !  (obtained from predicted positions) ???
     
     do i = 1, p_count
        particles(POSPREDX_PART_PROP,i) =  particles(POSX_PART_PROP,i) + &
             dtNew * (c2*particles(VELX_PART_PROP,i) + &
             particles(VELPREDX_PART_PROP,i))
        particles(POSPREDY_PART_PROP,i) = particles(POSY_PART_PROP,i) + &
             dtNew * (c2*particles(VELY_PART_PROP,i) + &
             particles(VELPREDY_PART_PROP,i))
        particles(POSPREDZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + &
             dtNew * (c2*particles(VELZ_PART_PROP,i) + &
             particles(VELPREDZ_PART_PROP,i))
     enddo
     call Grid_mapMeshToParticles(particles,&
          part_props, BLK_PART_PROP,p_count,&
          pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

     ! Update the particle positions using the predicted particle velocities
     !  (obtained from predicted positions)
     
     do i = 1, p_count
        particles(POSX_PART_PROP,i) = particles(POSPREDX_PART_PROP,i) + &
             dtNew * c4*particles(VELX_PART_PROP,i)
        particles(POSY_PART_PROP,i) = particles(POSPREDY_PART_PROP,i) + &
             dtNew * c4*particles(VELY_PART_PROP,i)
        particles(POSZ_PART_PROP,i) = particles(POSPREDZ_PART_PROP,i) + &
             dtNew * c4*particles(VELZ_PART_PROP,i)
     enddo
     ! end of correction section

  endif  ! end of split between doing correction and doing Euler timestep

    
  ! Map the updated gas velocity field onto the current particle positions to
  ! obtain the updated particle velocities.
  ! Note that this code is shared with the Euler step case.
  
  call Grid_mapMeshToParticles(particles,NPART_PROPS, BLK_PART_PROP,p_count,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)
     
  
  ! update initialization
  if (firstCall) then
     firstCall = .false.
  endif
  

  return
  
end subroutine pt_advanceEsti


