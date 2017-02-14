!!****if* source/Grid/GridMain/Chombo/UG/Grid_initDomain
!!
!! NAME
!!
!!  Grid_initDomain
!!
!!
!! SYNOPSIS
!!
!!  call Grid_initDomain(integer(IN)  :: myPE,
!!                       integer(IN)  :: numProcs,
!!                       logical(IN)  :: restart,
!!                       logical(INOUT) :: particlesInitialized)
!!
!!
!! DESCRIPTION
!! 
!!  Create the mesh, initialize all the mesh data structures
!!  and apply initial conditions
!!
!!
!! ARGUMENTS
!!
!!  myPE  : my Processor number
!!  numProcs : number of processors in the run
!!  restart : is true if the execution is starting from a checkpoint
!!            file, otherwise false.
!!  particlesInitialized : is true if particle positions were initialized before returning
!!                         from this routine
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

!!REORDER(4): solnData

subroutine Grid_initDomain(restart,particlesInitialized)

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_fillGuardCells, Grid_getBlkPtr, Grid_releaseBlkPtr

  use Grid_data, ONLY : gr_eosModeInit, gr_eosMode, gr_eosModeNow, gr_globalMe,&
       gr_globalNumProcs
  use Eos_interface, ONLY : Eos_wrapped
  use Simulation_interface, ONLY : Simulation_initBlock, &
    Simulation_initRestart
  use RadTrans_interface, ONLY : RadTrans_sumEnergy
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  logical, intent(in) :: restart
  logical, intent(inout) :: particlesInitialized

  integer :: blockID=1
  
  integer, dimension(2, MDIM) :: blkLimits,blkLimitsGC
  real, pointer:: solnData(:,:,:,:)
  
  call gr_createDomain()
  if(.not.restart) then
     !  zero data in case we don't initialize all data in
     !  Simulation_initBlock... in particular the total vs. internal
     !  energies can cause problems in the eos call that follows
     call Grid_getBlkPtr(blockID, solnData,CENTER)
     solnData = 0.0
     call Grid_releaseBlkPtr(blockID, solnData,CENTER)

     call Simulation_initBlock(blockID)
     
#ifdef ERAD_VAR
     ! Sum radiation energy density over all meshes. This call is
     ! needed for mesh replication.
     call RadTrans_sumEnergy(ERAD_VAR, 1, (/blockID/))
#endif

     call Grid_getBlkIndexLimits(1,blkLimits,blkLimitsGC)

     call Eos_wrapped(gr_eosModeInit, blkLimits, blockID)

  ! Do no call the EOS here any more when restarting, since those calls
  ! may (depending on the Eos implementation, or whether Eos has
  ! been called explicitly before the checkpoint was written)
  ! introduce small data differences.
  ! As long as checkpoint data is written in a state where the
  ! solution data is thermodynamically consistent, there is no
  ! need for the call here. - KW
  else
     call Simulation_initRestart()
  end if

  call Grid_fillGuardCells(CENTER_FACES, ALLDIR)

  gr_eosModeNow = gr_eosMode
  call gr_solversInit()


end subroutine Grid_initDomain
