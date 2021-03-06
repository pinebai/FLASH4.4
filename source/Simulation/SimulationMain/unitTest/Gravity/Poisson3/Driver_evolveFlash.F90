!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3/Driver_evolveFlash
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

  use Driver_data, ONLY: dr_globalMe, dr_nbegin, &
       dr_nend, dr_dt, dr_wallClockTimeLimit, &
       dr_simTime, dr_redshift, &
       dr_nstep, dr_dtOld, dr_dtNew, dr_restart, dr_elapsedWCTime
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, Particles_dump
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getListOfBlocks, Grid_updateRefinement
  use Hydro_interface, ONLY : Hydro
  use Gravity_interface, ONLY :  Gravity_potentialListOfBlocks, Gravity_unitTest
  !use IO_data, ONLY: io_justCheckpointed 
  use IO_interface, ONLY :IO_output,IO_outputFinal

  implicit none

#include "constants.h"
#include "Flash.h"

  integer   :: localNumBlocks
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr

  ! for unit test
  logical,save :: perfect = .true.
  character(len=20) :: fileName
  integer, parameter        :: fileUnit = 2
  integer,dimension(4) :: prNum
  integer :: temp,i
  
  temp = dr_globalMe
  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))
  open(fileUnit,file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe
  ! ------------ end of unitTest setup ---------------------------------------
  
  call Logfile_stamp( 'Entering evolution routine' , '[Driver_evolveFlash]')


  call Timers_start("evolution")

  if (dr_nend .GE. dr_nbegin) then


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
     !--------------------------------------------------------------------
     !- Start Physics Sequence
     !----
#ifdef DEBUG_DRIVER
     print*, 'going into Hydro/MHD'
#endif

  
     dr_simTime = dr_simTime + dr_dt

     call Timers_start("hydro")
#ifdef DEBUG_DRIVER
     print*,'going into hydro'
#endif
     call Hydro(blockCount, blockList, &
                dr_simTime, dr_dt, dr_dtOld, SWEEP_XYZ)

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

     dr_simTime = dr_simTime + dr_dt
     call Timers_start("hydro")
     call Hydro( blockCount, blockList, &
                dr_simTime, dr_dt, dr_dtOld, SWEEP_ZYX)
     call Timers_stop("hydro")



     call Timers_start("sourceTerms")
     call Driver_sourceTerms(blockCount, blockList, dr_dt)
     call Timers_stop("sourceTerms")

     call Timers_start("Particles_advance")
     call Particles_advance(dr_dt, dr_dt)
     call Timers_stop("Particles_advance")
     
     call Gravity_potentialListOfBlocks(blockCount,blockList)

     !----
     !- End Physics Sequence
     !--------------------------------------------------------------------
  end if

  ! Gravity unitTest calculations-------------------------------------
  call Gravity_unitTest(fileUnit,perfect)
  if (perfect) then
     write(fileUnit,'("all results conformed with expected values.")')
  else
     write(fileUnit,'("Failure in Gravity unitTest at time",G10.4)')dr_simTime
  end if


!! Eliminted all code beyond here, not needed for unit test -PMR
!! NO!  We'd actually like to SEE what was calculated. LBR
    !io_justCheckpointed = .false.
     
  ! ------------------------------- Gravity unitTest output
   close (fileUnit)   ! for Gravity_unitTest
  ! --------------------------------

  call Timers_stop("evolution")

  call Logfile_stamp( 'Exiting evolution routine' , '[Driver_evolveFlash]')

  call IO_outputFinal()

  call Timers_getSummary( dr_nstep)


  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")

  call Logfile_close()


  return
  
end subroutine Driver_evolveFlash



