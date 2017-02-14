!!****if* source/Particles/ParticlesMain/passive/unitTest/Particles_unitTest
!!
!! NAME
!!
!!  Particles_unitTest
!!
!! SYNOPSIS
!!
!!  call Particles_unitTest(
!!                          integer(in)     :: fileUnit,
!!                          logical(inout)  :: perfect  )
!!
!! DESCRIPTION
!!
!!  This subroutine advances particles by a step (possibly a half time
!!  step) and also advances a corresponding analytical solution for
!!  comparison.  To be used within tests that examine the accuracy
!!  of passive particle time advance.
!!
!! ARGUMENTS
!!
!!  
!!  fileUnit:  integer   number of file for output;
!!                       ignored in this implementation.
!!  perfect:   logical   should be set to .FALSE. to indicate
!!                       failure of a test, should otherwise
!!                       be left unchanged; ignored in this
!!                       implementation.
!!
!!
!!***

subroutine Particles_unitTest(fileUnit,perfect)
  
#include "Flash.h"
#include "constants.h"
  use Driver_interface, ONLY : Driver_getDt, Driver_getSimTime
  use Particles_interface, ONLY : Particles_advance
  use pt_utData, ONLY : pt_utAnalyticParticlePositions
  use pt_interface, ONLY : pt_utFakeParticlesAdvance,pt_utUpdateAnaPosns,&
       pt_utComputeError
  implicit none
  
  
  integer, intent(IN)           :: fileUnit ! Output to file
  logical, intent(INOUT)        :: perfect  ! Flag to indicate errors
  
  real,save :: dt, simTime, dtOld
  logical,save :: firstCall = .true.

  
  call Driver_getDt(dt)
  if(firstCall) then
     dtOld=dt
     firstCall = .false.
  end if
  call Driver_getSimTime(simTime)

  if (pt_utAnalyticParticlePositions) then
     call pt_utFakeParticlesAdvance(dtOld, dt, simTime,1)
  else
     call Particles_advance(dtOld, dt)
  end if
  
  call pt_utUpdateAnaPosns(dtOld, dt, simTime)
  call pt_utComputeError(dtOld, dt, simTime)
  
#ifdef DEBUG_DRIVER
  print*, 'return from Particles_advance '
#endif
  dtOld=dt

  return
  
end subroutine Particles_unitTest
