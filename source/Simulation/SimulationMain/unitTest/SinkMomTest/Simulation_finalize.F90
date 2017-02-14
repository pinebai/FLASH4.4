!!****if* source/Simulation/SimulationMain/unitTest/SinkMomTest/Simulation_finalize
!!
!! NAME
!!  Simulation_finalize
!!
!! SYNOPSIS
!!
!!  call Simulation_finalize()
!!
!! DESCRIPTION
!!
!!  This function xan clean up the Simulation unit, deallocates memory, etc.
!!  However, as nothing needs to be done, only this stub is included.
!!
!!  For this unitTest, the final determination of success or failure is
!!  processed here.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

#include "constants.h"

subroutine Simulation_finalize()

  use Logfile_interface, ONLY : Logfile_stamp
   use Simulation_data, ONLY : sim_globalMe, &
        sim_testInitialized, &
        sim_testRefVals, &
        sim_testLastVals, &
        sim_testWorstVals, &
        sim_testWorstWhen, &
        sim_momXTol, sim_momYTol, sim_momZTol, sim_massTol

  implicit none

  character (len = 4                ) :: charProcessorID
  character (len = MAX_STRING_LENGTH) :: fileName

  logical :: perfect, imperfect

  integer :: i
  integer :: fileUnit
  integer :: ut_getFreeFileUnit
  real    :: tols(4)

  if (sim_globalMe .GE. 16) return ! avoid writing lots of (empty) files if using many procs...

  perfect = .TRUE.

  tols(1) = sim_massTol
  tols(2) = sim_momXTol
  tols(3) = sim_momYTol
  tols(4) = sim_momZTol

!   ...Open the (processor specific) indicator file. This is a file that will contain
!      the success (failure) status of the unit test.
!
!
  write (charProcessorID,'(I4.4)') sim_globalMe

  fileUnit = ut_getFreeFileUnit ()
  fileName = "unitTest_" // charProcessorID

  open (fileUnit, file = fileName)
  if (sim_globalMe==MASTER_PE) then
     do i=1,min(4,1+N_DIM)
        imperfect = (abs(sim_testWorstVals(i)-sim_testRefVals(i)) > tols(i))
        perfect = (perfect .AND. .NOT. imperfect)
        if (.NOT. perfect) then
           select case (i)
           case(1)
              call Logfile_stamp(sim_testWorstWhen(i), &
                   '[Simulation_finalize] Total mass conservation - worst at')
           case(2)
              call Logfile_stamp(sim_testWorstWhen(i), &
                   '[Simulation_finalize] Total X momentum conservation - worst at')
           case(3)
              if (imperfect .OR. N_DIM>1) call Logfile_stamp(sim_testWorstWhen(i), &
                   '[Simulation_finalize] Total Y momentum conservation - worst at')
           case(4)
              if (imperfect .OR. N_DIM>2) call Logfile_stamp(sim_testWorstWhen(i), &
                   '[Simulation_finalize] Total Z momentum conservation - worst at')
           end select
        end if
        if (imperfect) then
           call Logfile_stamp(sim_testWorstVals(i)-sim_testRefVals(i), &
             '[Simulation_finalize] Test failed,    Error')
           call Logfile_stamp(sim_testRefVals(i), &
             '[Simulation_finalize] Test failed, reference value')
           call Logfile_stamp(sim_testLastVals(i), &
             '[Simulation_finalize] Test failed, last value')
           call Logfile_stamp(sim_testWorstVals(i), &
             '[Simulation_finalize] Test failed, worst value')
        else if (.NOT. perfect .AND. i .LE. 1+N_DIM) then
           call Logfile_stamp(sim_testWorstVals(i)-sim_testRefVals(i), &
             '[Simulation_finalize] Test succeeded, Error')
           call Logfile_stamp(sim_testRefVals(i), &
             '[Simulation_finalize] Test succeeded, reference value')
           call Logfile_stamp(sim_testLastVals(i), &
             '[Simulation_finalize] Test succeeded, last value')
           call Logfile_stamp(sim_testWorstVals(i), &
             '[Simulation_finalize] Test succeeded, worst value')
        end if
     end do
  end if


!   ...Final chores. The exact phrase 'all results conformed with expected values.' must
!      be included to the indicator files to ensure recognition of a successful unit test
!      run by the 'flashTest/lib/flashModule.py' script.
!
!  
  if (perfect) then
      write (fileUnit,'(a)') 'SUCCESS all results conformed with expected values.'
  else
      write (fileUnit,'(a)') 'FAILURE'
  end if

99 format(1x,A,":",4(2x,1PG14.7))
  if (sim_globalMe==MASTER_PE) then
     write (fileUnit,99) "Error is ",sim_testWorstVals-sim_testRefVals
     write (fileUnit,99) "Reference",sim_testRefVals
     write (fileUnit,99) "Last val ",sim_testLastVals
     write (fileUnit,99) "Worst val",sim_testWorstVals
     write (fileUnit,99) "Worst at ",sim_testWorstWhen
  end if
  close (fileUnit)


end subroutine Simulation_finalize
