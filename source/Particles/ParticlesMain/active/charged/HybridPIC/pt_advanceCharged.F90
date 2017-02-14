!!****if* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_advanceCharged
!!
!! NAME
!!
!!  pt_advanceCharged
!!
!! SYNOPSIS
!!
!!  call pt_advanceCharged(real(in)    :: dtold,
!!                        real(in)    :: dtnew,
!!                        real(inout) :: particlesunused(NPART_PROPS,p_countUnused),
!!                        integer(in) :: p_countunused)
!!
!! DESCRIPTION
!!
!!   Advances particles in time
!!
!! ARGUMENTS
!!
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particlesunused -- particles on which to operate
!!   p_countunused - the number of particles in the list to advance
!!
!!
!!
!!***

subroutine pt_advanceCharged (dtOld,dtNew,particlesUnused,p_countUnused)

  use Driver_interface, only : Driver_getDt, Driver_getSimTime
  use Particles_data, only: useParticles, pt_restart, particles, &
       pt_numLocal, pt_maxPerProc, pt_indexCount, pt_indexList
  use Grid_interface, only: Grid_moveParticles,Grid_fillGuardCells, &
       Grid_mapMeshToParticles, Grid_mapParticlesToMesh
  use Timers_interface, only : Timers_start, Timers_stop

  use pt_picInterface, ONLY : pt_picUpdatePositions, &
       pt_picSetCurrents, pt_picAdvanceB,&
       pt_picEfield, pt_picAdvanceVel, pt_picApplyBoundary
  implicit none

#include "Flash.h"
#include "constants.h"  
#include "Particles.h"

  real, intent(in)   :: dtOld, dtNew
  ! Cannot use these arguments since we want to add/delete particles
  integer,intent(in) :: p_countUnused
  real,dimension(NPART_PROPS,p_countUnused),intent(inout) :: particlesUnused

  integer, parameter :: add2grid = 0 ! =0 for zero grid before map to grid
  logical, save      :: first_call = .true.
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk
  real    :: dt, t
  logical :: globalMove, halfDt

  if (.not. useParticles) return
  call Timers_start("particles")

  call Driver_getDt(dt)
  call Driver_getSimTime(t)
  t  = t-dt
  
  if (first_call .and. .not. pt_restart) then
     globalMove = .true.
     call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,pt_numLocal, &
       pt_indexList, pt_indexCount, globalMove) 
     call pt_picApplyBoundary()
     halfDt=.false.
     call pt_picUpdatePositions(-dtNew, halfDt)   ! Moves the particles -dt/2
     globalMove = .true.
     call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,pt_numLocal, &
       pt_indexList, pt_indexCount, globalMove) 
     call pt_picApplyBoundary()

     first_call = .false.
  endif

  call Timers_start("move particles")
  halfDt=.false.
  call pt_picUpdatePositions(dt, halfDt)   ! move to time level n
  ! From here until the next call TMP?_PART_PROP stores position at n+1/2
  globalMove = .true.
  call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,pt_numLocal, &
       pt_indexList, pt_indexCount, globalMove) 
  call pt_picApplyBoundary()

  call pt_picSetCurrents ! prepare particles for MapToMesh

#define MAPPARTICLESTOMESH(propPart,vG,m) Grid_mapParticlesToMesh(particles,NPART_PROPS,pt_numLocal,pt_maxPerProc,propPart,vG,m)
  ! Deposit charge densities and ionic currents on the grid at time n
  call MAPPARTICLESTOMESH(CHARGE_PART_PROP, CDEN_VAR, add2grid)
  call MAPPARTICLESTOMESH(JX_PART_PROP, GJIX_VAR, add2grid)
  call MAPPARTICLESTOMESH(JY_PART_PROP, GJIY_VAR, add2grid)
  call MAPPARTICLESTOMESH(JZ_PART_PROP, GJIZ_VAR, add2grid)
  ! not needed for pic
  call MAPPARTICLESTOMESH(MASS_PART_PROP, PDEN_VAR, add2grid)
  call Grid_fillGuardCells( CENTER, ALLDIR)

  halfdt=.true.
  call pt_picUpdatePositions(dt, halfdt)   ! move to time level n+1/2

  call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,pt_numLocal, &
       pt_indexList, pt_indexCount, globalMove) 
  call pt_picApplyBoundary()

  call Timers_start("particle forces")
  call pt_picAdvanceB(dt)     ! Update B to n+1/2 according to CAM-CL
  call Timers_stop("particle forces")
  
  ! Deposit charge at n+1/2
  call MAPPARTICLESTOMESH(CHARGE_PART_PROP, CDEN_VAR, add2grid)

  call Grid_fillGuardCells( CENTER, ALLDIR)

  call pt_picEfield(GRBX_VAR, GRBY_VAR, GRBZ_VAR)
  ! From here until AdvanceVel, TMP?_PART_PROP stores velocity at n
  halfDt=.true.
  call pt_picAdvanceVel(dt, halfDt)
  call pt_picSetCurrents()
  ! Deposit ionic current at n+1/2
  call MAPPARTICLESTOMESH(JX_PART_PROP, GJIX_VAR, add2grid)
  call MAPPARTICLESTOMESH(JY_PART_PROP, GJIY_VAR, add2grid)
  call MAPPARTICLESTOMESH(JZ_PART_PROP, GJIZ_VAR, add2grid)
  
  call Grid_fillGuardCells( CENTER, ALLDIR)
  
  call pt_picEfield(GRBX_VAR, GRBY_VAR, GRBZ_VAR)
  halfDt=.false.
  call pt_picAdvanceVel(dt,halfDt)

  call Timers_stop("move particles")
  call Timers_stop ("particles")
  
end subroutine pt_advanceCharged

