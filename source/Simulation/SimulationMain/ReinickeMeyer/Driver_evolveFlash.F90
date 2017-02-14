!!****if* source/Simulation/SimulationMain/ReinickeMeyer/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!! This routine implements a directionally unsplit integrator for time 
!! advancement. A single step in this driver accounts for all 
!! multidimensional fluxes and updates conserved variables for the next
!! time step. The routine is a default driver for the Unsplit Staggered
!! Mesh (USM) MHD and Unsplit Hydro solvers.
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***

!#define DEBUG_DRIVER
#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif


subroutine Driver_evolveFlash()

  use Driver_data,         ONLY : dr_globalMe, dr_globalNumProcs, dr_nbegin, &
                                  dr_nend, dr_dt, dr_wallClockTimeLimit, &
                                  dr_tmax, dr_simTime, dr_redshift,      &
                                  dr_nstep, dr_dtOld, dr_dtNew,          &
                                  dr_simGeneration,                      &
                                  dr_restart, dr_elapsedWCTime,          &
                                  dr_redshiftOld, dr_useRedshift,        &
                                  dr_redshiftfinal,                      &
                                  dr_useSTS, dr_nuSTS, dr_nstepTotalSTS, &
                                  dr_dtSTS,dr_dt_subSTS,                 &
                                  dr_dtAdvect, dr_dtDiffuse,             &
                                  dr_useSTSforDiffusion,                 &
                                  dr_tstepChangeFactor,                  &
                                  dr_allowDtSTSDominate
  use Driver_interface,    ONLY : Driver_sourceTerms, Driver_computeDt, &
                                  Driver_getElapsedWCTime,&
                                  Driver_superTimeStep
  use Logfile_interface,   ONLY : Logfile_stamp, Logfile_close
  use Timers_interface,    ONLY : Timers_start, Timers_stop, &
                                  Timers_getSummary
  use Diffuse_interface,   ONLY : Diffuse
  use Particles_interface, ONLY : Particles_advance, Particles_dump
  use Grid_interface,      ONLY : Grid_getLocalNumBlks, &
                                  Grid_getListOfBlocks, &
                                  Grid_updateRefinement
  use Grid_interface,      ONLY : Grid_fillGuardCells
  use Hydro_interface,     ONLY : Hydro
  use Gravity_interface,   ONLY : Gravity_potentialListOfBlocks
  use IO_interface,        ONLY : IO_output,IO_outputFinal
  use Cosmology_interface, ONLY : Cosmology_redshiftHydro, &
                                  Cosmology_solveFriedmannEqn, &
                                  Cosmology_getRedshift
  use RadTrans_interface,  ONLY : RadTrans
  use Eos_interface,       ONLY : Eos_logDiagnostics
  use Simulation_interface, ONLY: Simulation_adjustEvolution
  use Simulation_interface, ONLY: Simulation_computeAnalytical
  use Profiler_interface, ONLY : Profiler_start, Profiler_stop

  implicit none

#include "constants.h"
#include "Flash.h"

  integer   :: localNumBlocks

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: sweepDummy

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(4,2) :: strBuff
  character(len=15) :: numToStr

  logical :: gridChanged
  logical :: endRun !Should we end our run on this iteration?
  logical :: shortenedDt !Is the last timestep being shortened to reach dr_tmax?

  ! for super-time-stepping
  integer :: nstepSTS
  real    :: dt_diffuse_temp
  logical :: useSTS_local
  integer :: nstepTotalSTS_local 

  integer, parameter :: driftUnk_flags = DRIFT_NO_PARENTS
  integer :: lb

  endRun = .false.

  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Profiler_start("FLASH_evolution")
  call Timers_start("evolution")

  do dr_nstep = dr_nBegin, dr_nend

     ! Initialize the local STS switch to user's setup.
     ! Later, useSTS_local can be turned-off when dr_dtDiffusion > dr_dtAdvect
     useSTS_local = dr_useSTS


     call dr_shortenLastDt(dr_dt, dr_simTime, dr_tmax, shortenedDt, 1)
     !!Step forward in time. See bottom of loop for time step calculation.
     call Grid_getLocalNumBlks(localNumBlocks)
     call Grid_getListOfBlocks(LEAF,blockList,blockCount)



     if (dr_globalMe == MASTER_PE) then

        write (numToStr(1:), '(I10)') dr_nstep
        write (strBuff(1,1), "(A)") "n"
        write (strBuff(1,2), "(A)") trim(adjustl(numToStr))
        
        write (numToStr(1:), "(1PE12.6)") dr_simTime
        write (strBuff(2,1), "(A)") "t"
        write (strBuff(2,2), "(A)") trim(adjustl(numToStr))

        if (.not. dr_useRedshift) then
           if (.not. dr_useSTS) then
              write (numToStr(1:), "(1PE12.6)") dr_dt
           else
              write (numToStr(1:), "(1PE12.6)") max(dr_dt,dr_dtSTS)
           endif
           write (strBuff(3,1), "(A)") "dt"
           write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))

           call Logfile_stamp( strBuff(1:3,:), 3, 2, "step")

        else

           write (numToStr(1:), "(F8.3)") dr_redshift
           write (strBuff(3,1), "(A)") "z"
           write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))

           if (.not. dr_useSTS) then
              write (numToStr(1:), "(1PE12.6)") dr_dt
           else
              write (numToStr(1:), "(1PE12.6)") max(dr_dt,dr_dtSTS)
           endif
           write (strBuff(4,1), "(A)") "dt"
           write (strBuff(4,2), "(A)") trim(adjustl(NumToStr))
           
           call Logfile_stamp( strBuff, 4, 2, "step")

        endif

     end if

     call Simulation_adjustEvolution(blockCount, blockList, dr_nstep, dr_dt, dr_simTime)

     if (dr_useSTSforDiffusion) then
        !! Note: This setup will use the STS for overcoming small diffusion time steps
        !!       (assuming dr_dtDiffuse < dr_dtAdvect) and accelerates time advancements
        !!       upto the orders of dr_dtAdvect.
        !! Do not allow to turn on a switch for the super time stepping when there is
        !! no diffusion (viscosity, conductivity, and magnetic resistivity) used.

        !! No need to use the super time stepping when diffusion time step is larger than
        !! advection time step.
        if (dr_useSTS) then

           ! initialize
           dt_diffuse_temp = dr_dtDiffuse

           ! determine if want to use STS
           if (dr_dtDiffuse > dr_dtAdvect) then
              useSTS_local = .false.
           else
              useSTS_local = .true.
           endif
        endif
     endif


     !! If super time stepping is not used, the subcycling is meaningless.
     if (.not. useSTS_local) then
        nstepTotalSTS_local = 1
     else
        nstepTotalSTS_local = dr_nstepTotalSTS
     endif

     !! Check CFL condition when the STS is used to accelerate diffusion time,
     !! i.e., when dr_useSTSforDiffusion = .true.
     if ((useSTS_local) .and. (dr_useSTSforDiffusion)) then
        nstepSTS = 1
        call Driver_superTimeStep(dr_dtDiffuse,dr_nuSTS,nstepSTS,nstepTotalSTS_local,dr_dt_subSTS)

        !dt_diffuse_temp = dr_dtDiffuse
        !! Reduce dt_diffuse_temp, thereby reduce dr_dt_subSTS, if dr_dt_subSTS exceeds advection time scale.
        !! dt_diffuse_temp is pre-calculated here and will be used for the STS subcycling later.

        if (.not. dr_allowDtSTSDominate) then
           do while (dr_dt_subSTS > dr_dtAdvect)

              dt_diffuse_temp = dt_diffuse_temp * (dr_dtAdvect/dr_dt_subSTS)

              call Driver_superTimeStep(dt_diffuse_temp,dr_nuSTS,nstepSTS,nstepTotalSTS_local,dr_dt_subSTS)
           enddo
        endif
     endif



     !! Begin subcycling of super time stepping algorithm
     dr_dtSTS=0.
     do nstepSTS = 1,nstepTotalSTS_local
        if (useSTS_local) then
           !! DEV - this dtOld needs to be checked for the very first iteration.
           !!       The same check is needed in Driver_verifyInitDt for a restart case.

           if (dr_useSTSforDiffusion) then
              !! Use the pre-calculated dr_diffuse_temp
              call Driver_superTimeStep(dt_diffuse_temp,dr_nuSTS,nstepSTS,nstepTotalSTS_local,dr_dt_subSTS)

           else
              !! Note: If using the STS for overcoming CFL limited time step dr_dtAdvection, 
              !! then Shi, Li, and Liang (IEE, 2006, vol 153, pp 55-60), suggest to use:
              !!       (1) 0.15 <= dr_nuSTS <= 0.25
              !!       (2) 2 <= dr_nstepTotalSTS <= 5
              !! Note that the time advancing is go beyond CFL limit by using the STS algorithm.
              call Driver_superTimeStep(dr_dtAdvect,dr_nuSTS,nstepSTS,nstepTotalSTS_local,dr_dt_subSTS)
           endif

           dr_dt = dr_dt_subSTS

           !! This is for a slow start-up.
           dr_dt = min(dr_dt,dr_dtOld*dr_tstepChangeFactor )
           dr_dtSTS = dr_dtSTS + dr_dt
        endif


        !--------------------------------------------------------------------
        !- Start Physics Sequence
        !----
#ifdef DEBUG_DRIVER
        print*,'going into Hydro/MHD'  ! DEBUG
        print*,'going into hydro myPE=',dr_globalMe
#endif
        
        ! 1. Cosmology-Friedmann Eqn.
        call Timers_start("cosmology")
        call Cosmology_solveFriedmannEqn(dr_simTime, dr_dt)
        call Timers_stop("cosmology")
        call Driver_driftUnk(__FILE__,__LINE__,driftUnk_flags)
        
        dr_simTime = dr_simTime + dr_dt
        dr_simGeneration = 0

        ! 2. Hydro/MHD
        call Timers_start("Hydro")
        call Hydro(blockCount, blockList,   &
                   dr_simTime, dr_dt, dr_dtOld,  sweepDummy)
        call Timers_stop("Hydro")
        call Driver_driftUnk(__FILE__,__LINE__,driftUnk_flags)
#ifdef DEBUG_DRIVER
        print*, 'return from Hydro/MHD timestep'  ! DEBUG
        print*,'returning from hydro myPE=',dr_globalMe
#endif

        ! 3. Diffusive processes: 
        !    Radiation, viscosity, conduction, & magnetic registivity
        call RadTrans(blockCount, blockList, dr_dt)
        call Diffuse(blockCount, blockList, dr_dt)
        call Driver_driftUnk(__FILE__,__LINE__,driftUnk_flags)
#ifdef DEBUG_DRIVER
        print*, 'return from Diffuse'
#endif

        ! 4. Add source terms:
        !    Stirring, flame, burning, heating, heat exchange, cooling, ionization,
        !    energy deposition, & deleptonization
        call Timers_start("sourceTerms")
        call Driver_sourceTerms(blockCount, blockList, dr_dt)
        call Timers_stop("sourceTerms")
        call Driver_driftUnk(__FILE__,__LINE__,driftUnk_flags)
#ifdef DEBUG_DRIVER
        print*,'done source terms'  ! DEBUG
        print*,'return from Driver_sourceTerms '  ! DEBUG
#endif

        ! 5. Advance Particles
        call Timers_start("Particles_advance")
        call Particles_advance(dr_dtOld, dr_dt)
        call Driver_driftUnk(__FILE__,__LINE__,driftUnk_flags)
#ifdef DEBUG_DRIVER
        print*, 'return from Particles_advance '  ! DEBUG
#endif

        ! 6. Calculate gravitational potentials
        call Timers_stop("Particles_advance")
        call Gravity_potentialListOfBlocks(blockCount,blockList)
        call Driver_driftUnk(__FILE__,__LINE__,driftUnk_flags)
#ifdef DEBUG_DRIVER
        print*, 'return from Gravity_potential '  ! DEBUG
#endif

        ! 7. Cosmology-Redshift
        call Timers_start("cosmology")
        call Cosmology_redshiftHydro( blockCount, blockList)
        call Timers_stop("cosmology")
        call Driver_driftUnk(__FILE__,__LINE__,driftUnk_flags)

        !! save for old dt
        dr_dtOld = dr_dt

     enddo !end of subcycling of super time stepping

     !----
     !- End Physics Sequence
     !--------------------------------------------------------------------

     !output a plotfile before the grid changes
     call Timers_start("IO_output")

     if (.not. useSTS_local) then
        call IO_output(dr_simTime, &
             dr_dt, dr_nstep+1, dr_nbegin, endRun, PLOTFILE_AND_PARTICLEFILE)
     else
        call IO_output(dr_simTime, &
             dr_dtSTS, dr_nstep+1, dr_nbegin, endRun, PLOTFILE_AND_PARTICLEFILE)
     endif
     call Timers_stop("IO_output")


     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement(dr_nstep, dr_simTime, gridChanged)
     call Timers_stop("Grid_updateRefinement")

     if (gridChanged) dr_simGeneration = dr_simGeneration + 1

     ! backup needed old
     if (.not. useSTS_local) dr_dtOld = dr_dt

#ifdef DEBUG_DRIVER
     print*, 'going into Driver_computeDt '  ! DEBUG
#endif
     ! calculate new
     call Driver_computeDt(dr_nbegin,  dr_nstep,      &
                           dr_simTime, dr_dtOld, dr_dtNew)
#ifdef DEBUG_DRIVER
     print*, 'return from Driver_computeDt '  ! DEBUG
#endif

     ! store new
     if (.not. useSTS_local) dr_dt = dr_dtNew

     call Timers_start("IO_output")
     if (.not. useSTS_local) then
        call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin,endRun,&
             CHECKPOINT_FILE_ONLY)
     else
        call IO_output(dr_simTime,dr_dtSTS,dr_nstep+1,dr_nbegin,endRun,&
             CHECKPOINT_FILE_ONLY)
     endif
     call Timers_stop("IO_output")

     call Eos_logDiagnostics(.TRUE.)


     !!*****************************************************************************
     !!  Evolution Loop -- check termination conditions
     !!*****************************************************************************

     !Exit if this step was handled specially as the last step
     if(shortenedDt) exit
     !Exit if a .dump_restart or .kill was found during the last step
     if(endRun) exit

     !! the simulation ends before nend iterations if
     !!  (i)   the simulation time is greater than the maximum time (tmax)
     !!  (ii)  the redshift falls below the minimum redshift  
     !!        (also called redshiftFinal) 
     !!  (iii) the wall clock time is greater than the maximum 
     !!        (wall_clock_time_max)


     !!Update redshift from Driver's POV.  Need this for exit condition. -PR   
     !!old redshift needed for accurate restarts.                              
     dr_redshiftOld = dr_redshift
     call Cosmology_getRedshift(dr_redshift)

     if (dr_simTime >= dr_tmax) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached max SimTime"
        endif
        exit
     end if

     call Driver_getElapsedWCTime(dr_elapsedWCTime)
     if (dr_elapsedWCTime >  dr_wallClockTimeLimit) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached max wall clock time"
        endif
        exit
     end if

     if (dr_redshift < dr_redshiftfinal .and. dr_useRedshift) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached redshiftfinal"
        endif
        exit
     end if

  enddo
  !The value of dr_nstep after the loop is (dr_nend + 1) if the loop iterated for
  !the maximum number of times.  However, we need to retain the value that
  !dr_nstep had during the last loop iteration, otherwise the number for nstep
  !that will be stored in a final checkpoint file will be wrong.
  dr_nstep = min(dr_nstep,dr_nend)

  !!******************************************************************************
  !! End of Evolution Loop
  !!******************************************************************************

  ! Compute the analytic solution for Reinicke-Meyer simulation:
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  do lb = 1, blockCount
     call Simulation_computeAnalytical(blockList(lb), dr_simTime)
  end do

  ! Write a checkpoint that will contain the analytic solution:
  call IO_writeCheckpoint

  call Timers_stop("evolution")
  call Profiler_stop("FLASH_evolution")
  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')
  !if a file termination, this may already be done.
  if(.NOT.endRun) call IO_outputFinal()
  call Timers_getSummary( max(0,dr_nstep-dr_nbegin+1))
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return
  
end subroutine Driver_evolveFlash
