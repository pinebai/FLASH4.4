!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/sim_printBlockData
!!
!!  NAME 
!!
!!   sim_printBlockData
!!
!!  SYNOPSIS
!!
!!   sim_printBlockData ()
!!
!!  DESCRIPTION
!!
!!   This routine prints the values of the relevant block variables used for the simulation.
!!   This routine is meant for checking purposes only. It only Each processor generates its own
!!   output file 
!!
!! ARGUMENTS
!!
!!***

subroutine sim_printBlockData ()

  use Simulation_data, ONLY : sim_baseName, &
                              sim_globalMe

  use ed_interface,    ONLY : ed_printBlockVariable
  use Grid_interface,  ONLY : Grid_getListOfBlocks

  implicit none

#include "Flash.h"
#include "constants.h"

  character (len = 4                ) :: charPID
  character (len = MAX_STRING_LENGTH) :: fileName

  integer :: block
  integer :: blockCount
  integer :: blockID
  integer :: fileUnit
  integer :: ut_getFreeFileUnit

  integer :: blockList (1:MAXBLOCKS)
!
!
!   ...Open the (processor specific) printout file.
!
!
  write (charPID,'(I4.4)') sim_globalMe

  fileUnit = ut_getFreeFileUnit ()
  fileName = trim (sim_baseName) // "BlockVariables" // charPID // ".txt"

  open (fileUnit, file = fileName)
!
!
!     ...Get the list of blocks on current processor.
!
!
  call Grid_getListOfBlocks (LEAF,   blockList, blockCount)
!
!
!     ...Lets have a look what density and temperature(s) are in each cell.
!
!
  do block = 1,blockCount

     blockID = blockList (block)

     call ed_printBlockVariable (blockID, DENS_VAR, fileUnit)
     call ed_printBlockVariable (blockID, TELE_VAR, fileUnit)
     call ed_printBlockVariable (blockID, TEMP_VAR, fileUnit)

  end do
!
!
!   ...Close the printout file.
!
!
  close (fileUnit)
!
!
!     ...Ready!
!
!
  return
end subroutine sim_printBlockData
