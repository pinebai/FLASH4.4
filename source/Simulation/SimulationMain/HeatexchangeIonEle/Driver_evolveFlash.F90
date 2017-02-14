!!****if* source/Simulation/SimulationMain/HeatexchangeIonEle/Driver_evolveFlash
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
!! This routine implements the Strang splitting scheme for time
!! advancement. A single step in the this driver 
!! includes two sweeps, the first one in order XYZ, and
!! the second one in order ZYX. This driver works with directionally
!! split operators only. The routine also controls the regridding of
!! the mesh if necessary and the simulation output.
!!
!!  This is a modification of the standard Driver_evolveFlash for testing 
!!  ion-electron heat exchange.
!!  
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


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif


subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe,  dr_nbegin, &
       dr_nend, dr_dt, dr_wallClockTimeLimit, &
       dr_tmax, dr_simTime, dr_simGeneration, dr_fSweepDir, dr_rSweepDir,&
       dr_nstep, dr_dtOld, dr_dtNew, dr_restart, dr_elapsedWCTime, &
       dr_redshiftInitial, dr_redshiftFinal, dr_redshift, dr_redshiftOld, &
       dr_useRedshift
  use Simulation_data, ONLY : sim_rhoInit, &
       sim_schemeOrder,&
       sim_maxTolCoeff0,sim_maxTolCoeff1,sim_maxTolCoeff2,sim_maxTolCoeff3

  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use Particles_interface, ONLY : Particles_advance
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getListOfBlocks, Grid_updateRefinement, Grid_computeVarDiff
  use Hydro_interface, ONLY : Hydro
  use Gravity_interface, ONLY :  Gravity_potentialListOfBlocks
  use IO_interface, ONLY :IO_output,IO_outputFinal
  use Cosmology_interface, ONLY : Cosmology_redshiftHydro, &
    Cosmology_solveFriedmannEqn, Cosmology_getRedshift
  use Simulation_interface, ONLY :   Simulation_computeAnalytical
  implicit none

#include "constants.h"
#include "Flash.h"

  integer   :: localNumBlocks, i

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(4,2) :: strBuff
  character(len=15) :: numToStr
  
  ! unitTest-like tracking of imperfection...
  logical ::  perfect = .true.

  logical :: gridChanged
  character(len=20) :: fileName
  integer, parameter        :: fileUnit = 2
  integer,dimension(4) :: prNum
  integer :: temp
  real :: maxActualDt = 0.0
  real :: maxError

  logical :: endRun !Should we end our run on this iteration?
  logical :: shortenedDt !Is the last timestep being shortened to reach dr_tmax?
  integer :: blk
  real    :: eNorm1, eNorm2, eNorm, temp1Norm, temp2Norm
  real    :: tNumNorm1, tNumNorm2, dNorm

  ! stays true if no errors are found
  perfect = .true.
  maxError = 0.0

  temp = dr_globalMe
  
  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))
  
  open(fileUnit,file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe

  endRun = .false.
  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')


  call Timers_start("evolution")

!!******************************************************************************
!! Start of Evolution Loop
!!******************************************************************************

  do dr_nstep = dr_nbegin, dr_nend


     call dr_shortenLastDt(dr_dt, dr_simTime, dr_tmax, shortenedDt, 2)
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

           write (numToStr(1:), "(1PE12.6)") dr_dt
           write (strBuff(3,1), "(A)") "dt"
           write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))

           call Logfile_stamp( strBuff(1:3,:), 3, 2, "step")

        else

           write (numToStr(1:), "(F8.3)") dr_redshift
           write (strBuff(3,1), "(A)") "z"
           write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))

           write (numToStr(1:), "(1PE12.6)") dr_dt
           write (strBuff(4,1), "(A)") "dt"
           write (strBuff(4,2), "(A)") trim(adjustl(NumToStr))
           
           call Logfile_stamp( strBuff, 4, 2, "step")

        endif



     end if

     !!--------------------------------------------------------------------
     !!- Start Physics Sequence
     !!--------------------------------------------------------------------
#ifdef DEBUG_DRIVER
     print*, 'going into Hydro/MHD'
#endif

     call Timers_start("cosmology")
     call Cosmology_solveFriedmannEqn(dr_simTime, dr_dt)
     call Timers_stop("cosmology")

     maxActualDt = max(maxActualDt, dr_dt)

     dr_simTime = dr_simTime + dr_dt
     dr_simGeneration = 0

     call Timers_start("hydro")
#ifdef DEBUG_DRIVER
     print*,'going into hydro'
#endif
     call Hydro( blockCount, blockList, &
                dr_simTime, dr_dt, dr_dtOld, dr_fSweepDir)

     call Timers_stop("hydro")

     
#ifdef DEBUG_DRIVER
     print*, 'return from Hydro/MHD timestep'
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

     call Timers_start("cosmology")
     call Cosmology_redshiftHydro( blockCount, blockList)
     call Timers_stop("cosmology")

!!******************************************************************************
!!Second "half-step" of the evolution loop
!!******************************************************************************

     call Timers_start("cosmology")
     call Cosmology_solveFriedmannEqn(dr_simTime, dr_dt)
     call Timers_stop("cosmology")

     dr_simTime = dr_simTime + dr_dt
     dr_simGeneration = 0
     call Timers_start("hydro")
     call Hydro( blockCount, blockList, &
                dr_simTime, dr_dt, dr_dtOld, dr_rSweepDir)
     call Timers_stop("hydro")
  
     call Timers_start("sourceTerms")
     call Driver_sourceTerms(blockCount, blockList, dr_dt)
     call Timers_stop("sourceTerms")

     call Timers_start("Particles_advance")
     call Particles_advance(dr_dt, dr_dt)
     call Timers_stop("Particles_advance")
     
     call Gravity_potentialListOfBlocks(blockCount,blockList)

     call Timers_start("cosmology")
     call Cosmology_redshiftHydro( blockCount, blockList)
     call Timers_stop("cosmology")

     !--------------------------------------------------------------------
     !- End Physics Sequence -- Start Simulation Bookkeeping
     !--------------------------------------------------------------------
     do blk=1, localNumBlocks
      call Simulation_computeAnalytical(blk,  dr_simTime)
     end do
   
     call Grid_computeVarNorm(0, 2, T1AN_VAR,temp1Norm, 1)
     call Grid_computeVarNorm(0, 2, T2AN_VAR,temp2Norm, 1)

     call Grid_computeVarDiff(-1, T1AN_VAR, TION_VAR, T1ER_VAR)
     call Grid_computeVarDiff(-1, T2AN_VAR, TELE_VAR, T2ER_VAR)

     call Grid_computeVarNorm(0, 2, T1ER_VAR,eNorm1, 1)
     call Grid_computeVarNorm(0, 2, T2ER_VAR,eNorm2, 1)

     call Grid_computeVarNorm(0, 2, TION_VAR,tNumNorm1, 1)
     call Grid_computeVarNorm(0, 2, TELE_VAR,tNumNorm2, 1)
     call Grid_computeVarNorm(0, 2, DENS_VAR,dNorm,  1)
     dNorm = dNorm / sim_rhoInit

     eNorm = sqrt(eNorm1**2 + eNorm2**2)

!!$     print*,'eNormL',eNorm, eNorm1, eNorm2

     eNorm1 = eNorm1 / temp1Norm
     eNorm2 = eNorm2 / temp2Norm
     eNorm = sqrt(eNorm1**2 + eNorm2**2)

     maxError = max(maxError,eNorm)
!!$     print*,'rNormL',eNorm, eNorm1, eNorm2
999  format(1x,6(1PG23.16))
!!$     print*, dr_simTime,eNorm, temp1Norm/dNorm, temp2Norm/dNorm,tNumNorm1/dNorm, tNumNorm2/dNorm
     write(fileUnit,999) dr_simTime,eNorm, temp1Norm/dNorm, temp2Norm/dNorm,tNumNorm1/dNorm, tNumNorm2/dNorm

     !output a plotfile before the grid changes
     call Timers_start("IO_output")
     call IO_output( dr_simTime, &
          dr_dt, dr_nstep+1, dr_nbegin, endRun, PLOTFILE_AND_PARTICLEFILE)
     call Timers_stop("IO_output")

     !!if (itemp_limit) .eq. 1) call Hydro_timstepPrecompute()

     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement( dr_nstep, dr_simTime, gridChanged)
     call Timers_stop("Grid_updateRefinement")

     if (gridChanged) dr_simGeneration = dr_simGeneration + 1

     dr_dtOld = dr_dt                     ! backup needed old 
     ! calculate new
     
     call Timers_start("compute dt")
     call Driver_computeDt(dr_nbegin, dr_nstep, &
                         dr_simTime, dr_dtOld, dr_dtNew)
     call Timers_stop("compute dt")
     dr_dt = dr_dtNew                                        ! store new


     

     !!-----------------------------------------------------------------
     !! Output for current step in evolution
     !!-----------------------------------------------------------------

     call Timers_start("IO_output")
     call IO_output( dr_simTime, &
          dr_dt, dr_nstep+1, dr_nbegin, endRun, CHECKPOINT_FILE_ONLY)
     call Timers_stop("IO_output")

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
        end if
        exit
     end if
     
     call Driver_getElapsedWCTime(dr_elapsedWCTime)
     if (dr_elapsedWCTime >  dr_wallClockTimeLimit) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached max wall clock time"
        end if
        exit
     end if

     if (dr_redshift < dr_redshiftfinal .and. dr_useRedshift) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached redshiftfinal"
        end if
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



  call Timers_stop("evolution")

  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')

  !if a file termination, this may already be done.
  if(.NOT.endRun) call IO_outputFinal( )

  call Timers_getSummary( max(0,dr_nstep-dr_nbegin+1))

  !finish unit test write out file

!  call sim_getComputedError(maxError)

  print*,dr_globalMe,maxError,sim_maxTolCoeff0,&
       sim_maxTolCoeff1 * maxActualDt   ,&
       sim_maxTolCoeff2 * maxActualDt**2,&
       sim_maxTolCoeff3 * maxActualDt**3
  if (sim_schemeOrder == 0) then
     if (maxError > sim_maxTolCoeff0) perfect = .FALSE.
  end if
  if (sim_schemeOrder == 1) then
     if (maxError > sim_maxTolCoeff1 * maxActualDt) perfect = .FALSE.
  end if
  if (sim_schemeOrder == 2) then
     if (maxError > sim_maxTolCoeff2 * maxActualDt**2) perfect = .FALSE.
  end if
  if (sim_schemeOrder == 3) then
     if (maxError > sim_maxTolCoeff3 * maxActualDt**3) perfect = .FALSE.
  end if

  if (perfect) then
    write(fileUnit,'("Heatexchange unitTest PASSED!  all results conformed with expected values.")')
    write(*,'("Heatexchange unitTest PASSED!  all results conformed with expected values.")')
    call Logfile_stamp( "Heatexchange unitTest PASSED!")
  else
    write(fileUnit,'("Heatexchange unitTest FAILED!")')
    write(*,'("Heatexchange unitTest FAILED!")')
    call Logfile_stamp( "Heatexchange unitTest FAILED!")
  endif

  close(fileUnit)

  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")

  call Logfile_close()


  return
  
end subroutine Driver_evolveFlash



