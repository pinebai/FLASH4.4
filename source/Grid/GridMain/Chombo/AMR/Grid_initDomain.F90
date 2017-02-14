!!****if* source/Grid/GridMain/Chombo/AMR/Grid_initDomain
!!
!! NAME
!!
!!  Grid_initDomain
!!
!!
!! SYNOPSIS
!!
!!  Grid_initDomain(logical(IN)  :: restart,
!!                  logical(INOUT) :: particlesInitialized)
!!
!!
!! DESCRIPTION
!!
!!  Create the mesh, initialize all the mesh data structures
!!  and apply initial conditions
!!
!!  Initially very few blocks are created (number supplied at runtime).
!!  then user-defined refinment critera is applied to determine the 
!!  blocks that need to be refined and derefined.  
!!
!!  After the refinement, the newly created child blocks are filled via
!!  prolongation from the coarse parents.  This prolongation step can use
!!  prolongation routine supplied with paramesh or defined by the user.
!!
!!  Once the prolongation is done, the guardcells are filled.  Finally, the
!!  EOS is called on the block interiors to make them thermodynamically
!!  consistent.
!!
!! ARGUMENTS
!!
!!  restart : is true if the execution is starting from a checkpoint
!!            file, otherwise false.
!!  particlesInitialized : is true if particle positions were initialized before returning
!!                         from this routine
!!
!! NOTES
!!  When restarting from a checkpoint file, block interiors are assumed to
!!  have been filled when this interface is called. The EOS is not called on
!!  the block interiors in this implementation for use with Paramesh. It is
!!  assumed that data is already thermodynamically consistent, because
!!  that is how data are written to checkpoint files.
!!
!!***


subroutine Grid_initDomain(restart,particlesInitialized)

  use Grid_interface, ONLY : Grid_fillGuardCells
  use Simulation_interface, ONLY : Simulation_initRestart
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_eosModeNow, gr_eosMode
  implicit none

#include "Flash.h"
#include "constants.h"


  logical, intent (IN) :: restart
  logical, intent(INOUT) :: particlesInitialized

  integer :: i
  integer :: oldLocalNumBlocks !need this if running with particles

  integer :: blockID

  integer :: iblk

  integer ,dimension(MAXBLOCKS) :: blkList
  integer :: blkCount

  call gr_createDomain()
  if(.not.restart) then
     call gr_expandDomain(particlesInitialized)
     gr_eosModeNow = gr_eosMode !may be different from gr_eosModeInit
  else
     call Simulation_initRestart()
  end if

  call Grid_fillGuardCells(CENTER_FACES, ALLDIR)

  call gr_solversInit()
  !print *, 'finished with Grid_initDomain'
  
end subroutine Grid_initDomain
