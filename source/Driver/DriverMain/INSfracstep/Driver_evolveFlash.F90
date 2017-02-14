!!****if* source/Driver/DriverMain/INSfracstep/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  call Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!!  This is the main global driver for simulations that are:
!!      Spatially refined, State form, strang split
!!
!!  DOC: Driver_evolveFlash needs more explanation 
!!
!! NOTES
!!
!!  variables that begin with "dr_" like, dr_globalMe or dr_dt, dr_beginStep
!!  are stored in the data fortran module for the Driver unit, Driver_data.
!!  The "dr_" is meant to indicate that the variable belongs to the Driver Unit.
!!  all other normally named variables i, j, etc are local variables.
!!
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs, dr_nbegin,       &
                         dr_nend, dr_dt,                        &
                         dr_tmax, dr_simTime, dr_redshift,      &
                         dr_nstep, dr_dtOld, dr_dtNew,          &
                                  dr_simGeneration,                      &
                         dr_restart
  use IncompNS_interface, ONLY : IncompNS, IncompNS_velomg2center
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt
  use Logfile_interface,ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
                               Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, Particles_dump
  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_updateRefinement
!  use Gravity_interface, ONLY :  Gravity_potentialListOfBlocks
  use IO_interface,      ONLY : IO_output,IO_outputFinal
!  use Cosmology_interface, ONLY:  Cosmology
  use Simulation_interface, ONLY: Simulation_adjustEvolution

  implicit none

#include "constants.h"
#include "Flash.h"

  integer   :: localNumBlocks

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: sweepDummy
  
  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr

  logical :: gridChanged
  logical :: endRunPl !Should we end our run on this iteration, based on conditions detected by the IO unit?
  logical :: endRun !Should we end our run on this iteration, based on conditions detected by the IO unit?
  logical :: endRunWallClock !Should we end our run on this iteration, based on wall clock time?

  endRunPl = .false.
  endRun = .false.

  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")


  ! Initial Timestep:
  ! backup needed old
  dr_dtOld = dr_dt

  ! calculate new
  call Driver_computeDt(dr_nbegin,  dr_nstep,      &
                        dr_simTime, dr_dtOld, dr_dtNew)
  ! store new
  dr_dt = dr_dtNew

  do dr_nstep = dr_nBegin, dr_nend
     
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
        
        write (numToStr(1:), "(1PE12.6)") dr_dt
        write (strBuff(3,1), "(A)") "dt"
        write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))
        
        call Logfile_stamp( strBuff, 3, 2, "step")
     end if

     call Simulation_adjustEvolution(blockCount, blockList, dr_nstep, dr_dt, dr_simTime)

     !--------------------------------------------------------------------
     !- Start Physics Sequence
     !----
     dr_simTime = dr_simTime + dr_dt
     dr_simGeneration = 0

#ifdef DEBUG_DRIVER
     print*, 'going into IncompNS'
#endif
     call Timers_start("IncompNS")
     call IncompNS(blockCount, blockList,   &
                   dr_simTime, dr_dt, dr_dtOld,  sweepDummy)
     call Timers_stop("IncompNS")
#ifdef DEBUG_DRIVER
  print*, 'return from IncompNS timestep'
#endif

     call Timers_start("Particles_advance")
     call Particles_advance(dr_dtOld, dr_dt)
#ifdef DEBUG_DRIVER
     print*, 'return from Particles_advance '
#endif
     call Timers_stop("Particles_advance")     

     !----
     !- End Physics Sequence
     !--------------------------------------------------------------------

     !output a plotfile before the grid changes
     call Timers_start("IO_output")

     call IO_output(dr_simTime, &
             dr_dt, dr_nstep+1, dr_nbegin, endRunPl, PLOTFILE_AND_PARTICLEFILE)
     call Timers_stop("IO_output")


     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement(dr_nstep, dr_simTime, gridChanged)
     call Timers_stop("Grid_updateRefinement")

     if (gridChanged) dr_simGeneration = dr_simGeneration + 1

     if (dr_globalMe .eq. MASTER_PE) then
        write(*,*) ' '
        write(*,'(I6,A,g16.8,A,g16.8)') dr_nstep,&
                ', TimeStep= ',dr_dt,', SimTime= ', dr_simTime
     endif

     if (dr_globalMe .eq. MASTER_PE) &
     write(*,*) '###############################################################################'

     ! Compute next step dt:
     ! backup needed old
     dr_dtOld = dr_dt

     ! calculate new
     call Driver_computeDt(dr_nbegin,  dr_nstep,      &
                           dr_simTime, dr_dtOld, dr_dtNew)
     ! store new
     dr_dt = dr_dtNew

     ! Velocities and Omg to Center variables
     ! In your Simulation Config set REQUIRES physics/IncompNS/IncompNSExtras
     ! Note this will add velocity and vorticity variables to your CENTER data structure. 
     ! Average Velocities and Vorticity to cell-centers
     call IncompNS_velomg2center(blocklist,blockcount)


     call Timers_start("io")
     call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin,endRun,&
             CHECKPOINT_FILE_ONLY)
     call Timers_stop("io")
     endRun = (endRunPl .OR. endRun)

     !!*****************************************************************************
     !!  Evolution Loop -- check termination conditions
     !!*****************************************************************************

     !Exit if a .dump_restart or .kill was found during the last step
     if(endRun) exit

     !! the simulation ends before nend iterations if
     !!  (i)   the simulation time is greater than the maximum time (tmax)
     !!  (ii)  the redshift falls below the minimum redshift  
     !!        (also called zfinal)
     !!  (iii) the wall clock time is greater than the maximum 
     !!        (wall_clock_time_max)

     if (dr_simTime >= dr_tmax) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached max SimTime"
        endif
        exit
     end if

     call dr_wallClockLimitExceeded(endRunWallClock)
     if (endRunWallClock) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached max wall clock time"
        endif
        exit
     end if

  enddo
  dr_nstep = min(dr_nstep,dr_nend)

  !!******************************************************************************
  !! End of Evolution Loop
  !!******************************************************************************

  call Timers_stop("evolution")
  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')
  if(.NOT.endRun) call IO_outputFinal( )
  call Timers_getSummary( max(0,dr_nstep-dr_nbegin+1))
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return
  
end subroutine Driver_evolveFlash
