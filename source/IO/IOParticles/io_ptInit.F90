!!****if* source/IO/IOParticles/io_ptInit
!!
!! NAME
!!  io_ptInit
!!
!! SYNOPSIS
!!
!!  io_ptInit()
!!
!! DESCRIPTION
!!
!!  Perform IOParticles initialization
!!
!!
!!  The IOParticles unit uses a number of runtime parameters to determine
!!  if, when, and how various particle output files need to be written.
!!
!!  To determine exactly which runtime parameters control these
!!  files, please check the Config file in IO/IOParticles or the
!!  setup_params file in the object directory.
!!
!!
!! 
!! ARGUMENTS
!!
!!
!!
!!***

subroutine io_ptInit()

  use IOParticles_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  use IO_data, ONLY : io_restart, io_summaryOutputOnly
  use IO_interface, ONLY : IO_getPrevScalar
  implicit none

#include "Flash.h"
#include "constants.h"



  integer :: step
  real :: simTime
  integer :: error
  real :: zinitial

  call RuntimeParameters_get('particleFileIntervalTime', io_particleFileIntervalTime)
  call RuntimeParameters_get('particleFileIntervalStep', io_particleFileIntervalStep)
  call RuntimeParameters_get('particleFileNumber', io_particleFileNumber)
  call RuntimeParameters_get('useParticles', useParticles)
  call RuntimeParameters_get('particleFileIntervalZ', io_particleFileIntervalZ)
  call RuntimeParameters_get('writeParticleSubset', io_writeParticleSubset)
  call RuntimeParameters_get('writeParticleAll', io_writeParticleAll)

  io_dumpParticleFileExist = .false.
  
 
  if(io_restart) then
     call IO_getPrevScalar("nstep", step, error)
     if (error /= NORMAL) call Driver_abortFlash("ERROR: Can't find nstep.")
     call IO_getPrevScalar("time", simTime, error)
     if (error /= NORMAL) call Driver_abortFlash("ERROR: Can't find time.")
     call IO_getPrevScalar('nextParticleFileZ', io_nextParticleFileZ, error)
     if(error /= NORMAL) then
        if(error == NOTFOUND) then
           io_nextParticleFileZ = HUGE(1.)
        else
           call Driver_abortFlash("ERROR: Error in looking up nextParticleFileZ scalar!")
        end if
     end if
     call IO_getPrevScalar('nextParticleFileTime', io_nextParticleFileTime, error)
     if(error /= NORMAL) then
        if(error == NOTFOUND) then
           ! backward compatibility for checkpoints which don't have this scalar
           if ( io_particleFileIntervalTime > 0.e0 ) then
              ! assume we were trying for an even-tempered sequence from the beginning
              io_nextParticleFileTime = ceiling(simTime/io_particleFileIntervalTime)*io_particleFileIntervalTime
           else
              io_nextParticleFileTime = 0.e0
           endif
        else
           call Driver_abortFlash("ERROR: Error in looking up nextParticleFileTime scalar!")
        end if
     end if
  else
     call RuntimeParameters_get("nbegin", step)
     call RuntimeParameters_get("tinitial",simTime)
     call RuntimeParameters_get("zInitial", zinitial)
     io_nextParticleFileTime = 0.0 !any number to make sure logical expression below is valid, overwritten below
  end if
 


  
  
  !initialize next step at which to write a particle file
  
  if(io_particleFileIntervalStep > 0) then
     io_nextParticleFileStep = step + io_particleFileIntervalStep
  else
     io_nextParticleFileStep = 0
  end if

  !calculate next simulation time in which to 
  !write a Particle file
  if (.not. io_restart .or. (io_nextParticleFileTime <= simTime) .or. &
       (io_nextParticleFileTime > simTime + io_particleFileIntervalTime)) then
     if ( io_particleFileIntervalTime > 0.e0 ) then
        io_nextParticleFileTime = simTime + io_particleFileIntervalTime

     else
        io_nextParticleFileTime = 0.e0
     end if
  end if

 !If we are using redshift controls, set the next redshift we should output.
  if(.not. io_restart)then
     if( io_particleFileIntervalZ < HUGE(1.)) then
        io_nextParticleFileZ = zinitial - io_particleFileIntervalZ
     else
        io_nextParticleFileZ = HUGE(1.)
     end if
  end if


  if (io_summaryOutputOnly) then
     !Ensure that the variables nextParticleFileZ and nextParticleFileTime
     !retain their original values so that they will contain expected values
     !if we write an emergency checkpoint file, e.g. with .dump_restart.

     !IO_writeParticles.F90 particle file conditions are shown in parentheses:
     !------------------------------------------------------------------------
     !(io_nextParticleFileStep == nstep)
     io_nextParticleFileStep = -1 !Normally set from io_particleFileIntervalStep
     !(io_particleFileIntervalTime > 0.e0 .AND. lsimTime >= io_nextParticleFileTime)
     io_particleFileIntervalTime = -1.0 !No need to modify io_nextParticleFileTime
     !(io_particleFileIntervalZ < HUGE(1.) .AND. currentRedshift <= io_nextParticleFileZ)
     io_particleFileIntervalZ = HUGE(1.)
     !(io_justCheckpointed)
     !Note that io_justCheckpointed can never be set to .true. when
     !io_summaryOutputOnly is .true., see the values assigned in IO_init.F90
     !when io_summaryOutputOnly is .true. and the logic in IO_output.F90.
  end if

end subroutine io_ptInit
