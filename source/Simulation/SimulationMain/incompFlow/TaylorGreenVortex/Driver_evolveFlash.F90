!!****if* source/Simulation/SimulationMain/incompFlow/TaylorGreenVortex/Driver_evolveFlash
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
!!$#define DEDUG_ALL

#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_nbegin,       &
                         dr_nend, dr_dt, dr_wallClockTimeLimit, &
                         dr_tmax, dr_simTime, dr_redshift,      &
                         dr_nstep, dr_dtOld, dr_dtNew,          &
                                  dr_simGeneration,                      &
                         dr_restart, dr_elapsedWCTime
  use IncompNS_interface, ONLY : IncompNS, IncompNS_velomg2center
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime
  use Logfile_interface,ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
                               Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, Particles_dump
  use Grid_interface,    ONLY : Grid_getLocalNumBlks,  &
                                Grid_getListOfBlocks,  &
                                Grid_updateRefinement, &
                                Grid_getBlkPtr,        &
                                Grid_releaseBlkPtr
  use Gravity_interface, ONLY :  Gravity_potentialListOfBlocks


  use IO_interface,      ONLY : IO_output,IO_outputFinal
  use Simulation_interface, ONLY: Simulation_adjustEvolution

#ifdef WRITE_TECPLOT
  use IO_data , ONLY : io_plotFileNumber, IO_plotFileIntervalStep, IO_plotFileIntervalTime
  use IncompNS_data, only : ins_cflflg
#endif

  implicit none

#include "constants.h"
#include "Flash.h"

  ! Local Variables:
  integer :: localNumBlocks
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: sweepDummy

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr

  logical :: endRun

  logical :: tecplot_flg
  integer :: lb, blockID, count, firstfileflag

  integer :: ibd,NumBodies

!-----------------------------------------------------------------------------------------

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

  if (dr_globalMe == MASTER_PE) write(*,*) 'dr_dt ===',dr_dt

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

#ifdef WRITE_TECPLOT
  count = io_plotFileNumber-1
  if (.not. dr_restart) then
     firstfileflag = 0
     dr_nstep = 0
     count    = 0
     call outtotecplot(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
                       0.0,blockList,blockCount,firstfileflag)
  endif
  firstfileflag = 1
#endif

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


     call Timers_start("sourceTerms")
     call Driver_sourceTerms(blockCount, blockList, dr_dt)
     call Timers_stop("sourceTerms")
#ifdef DEBUG_DRIVER
     print*,'done source terms'
     print*, 'return from Drivers_sourceTerms '
#endif
     call Timers_start("Particles_advance")
     call Particles_advance(dr_dtOld, dr_dt)
#ifdef DEBUG_DRIVER
     print*, 'return from Particles_advance '
#endif
     call Timers_stop("Particles_advance")

     call Gravity_potentialListOfBlocks(blockCount,blockList)
#ifdef DEBUG_DRIVER
     print*, 'return from Gravity_potential '
#endif

     !----
     !- End Physics Sequence
     !--------------------------------------------------------------------

!!$     call Timers_start("Grid_updateRefinement")
!!$     call Grid_updateRefinement( dr_nstep, dr_simTime)
!!$     call Timers_stop("Grid_updateRefinement")

     !----
     !- Output results and plot files:
     !--------------------------------------------------------------------
#ifdef WRITE_TECPLOT
     ! Output to Tecplot
     if (ins_cflflg .eq. 1) then ! Constant cfl
       if (dr_nstep .gt. 1) then
       tecplot_flg = (1/IO_plotFileIntervalTime*MOD(dr_simtime,IO_plotFileIntervalTime) .le. &
                      dr_dt/IO_plotFileIntervalTime)
       else
       tecplot_flg = .false.
       endif
     else                        ! Constant timestep
       tecplot_flg = (MOD(dr_nstep,IO_plotFileIntervalStep) .eq. 0)
     endif

     if (tecplot_flg) then
        ! Write to Grid to Tecplot:
        count = count + 1
        call outtotecplot(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
                          0.0,blockList,blockCount,firstfileflag)

        if (count .gt. 0) firstfileflag = 1
     endif
#endif

     ! Average Velocities and Vorticity to cell-centers
     call IncompNS_velomg2center(blocklist,blockcount)

     ! Flash Output Routine:
     call Timers_start("io")
     call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin,endRun)
     call Timers_stop("io")

     !--------------------------------------------------------------------
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

     if(dr_nstep .eq. dr_nend) then
        ! Check Taylor Green Vortex errors:=
        call Grid_getListOfBlocks(LEAF,blockList,blockCount)
        call compute_errors(dr_globalMe,blockCount,blockList)
     endif

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

     call Driver_getElapsedWCTime(dr_elapsedWCTime)
     if (dr_elapsedWCTime >  dr_wallClockTimeLimit) then
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
