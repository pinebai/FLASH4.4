!!****if* source/Simulation/SimulationMain/unitTest/IO/IOTypes/Driver_evolveFlash
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
!! This is a very simple version of the Driver_evolveFlash routine,
!! that is meant to be used exclusively with IO Unit
!! testing. There is no time advancement involved here.
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

subroutine Driver_evolveFlash()

  use Driver_data, ONLY:   dr_nbegin,  dr_restart, dr_initialSimTime, &
       dr_globalMe, dr_globalNumProcs
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use physicaldata, ONLY : unk
  use IO_interface, ONLY : IO_writeCheckpoint, IO_writePlotfile, &
                           IO_writeParticles
  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: iOut, var, fail
  integer, parameter :: comm = MESH_COMM

  iOut = 2



  !initialize the fake grid variable with dummy values
  do var = 1, NUNK_VARS
     
     unk(var, :,:,:,:) = var*1.0

  end do
     
  call Timers_start("IO_writeCheckpoint")
  call IO_writeCheckpoint()
  call Timers_stop("IO_writeCheckpoint")

  call Timers_start("IO_writePlotfile")
  call IO_writePlotfile()
  call Timers_stop("IO_writePlotfile")

  call Timers_start("IO_writeParticles")
  call IO_writeParticles( .false.)
  call Timers_stop("IO_writeParticles")


  call Timers_getSummary(0)


  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")

  call Logfile_close()


  !I read the values using a single processor so only do the 
  !comparison if a small number of processors.
  if (dr_globalNumProcs <= 4) then
     call check_mesh_values(dr_globalMe, fail)
  else
     fail = 0
  end if


  if (dr_globalMe == 0) then
     open(iOut,file='unitTest_0000')
     if (fail == 0) then
        write(iOut,'("all results conformed with expected values.")')
     else
        write(iOut,'("test failed.")')
     end if
     close(iOut)
  end if

end subroutine Driver_evolveFlash


subroutine check_mesh_values(me, fail)
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  implicit none
  integer, intent(IN) :: me
  integer, intent(OUT) :: fail
  integer :: i, j, ptr
  character (len=4) :: unklabels(UNK_VARS_BEGIN:UNK_VARS_END), label
  character, dimension(4*NUNK_VARS) :: packedLabels
  integer, external :: sim_expected_mesh_values

  ptr = 1
  do i=UNK_VARS_BEGIN,UNK_VARS_END
     !The index i contains the expected attribute values.
     call Simulation_mapIntToStr(i, unklabels(i), MAPBLOCK_UNK)
     label = unkLabels(i)
     do j = 1, 4
        packedLabels(ptr) = label(j:j)
        ptr = ptr + 1
     end do
  end do
  fail = sim_expected_mesh_values(me, packedLabels)
end subroutine check_mesh_values
