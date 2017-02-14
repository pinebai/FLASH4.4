!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_potentialListOfBlocks
!!
!!  NAME 
!!
!!     TreeRay_potentialListOfBlocks
!!
!!  SYNOPSIS
!!
!!  TreeRay_potentialListOfBlocks(integer(IN) :: blockCount,
!!                                integer(IN) :: blockList(blockCount))
!!
!!  DESCRIPTION 
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to calculate potential
!!
!! SIDE EFFECTS
!!
!!
!!
!!***


subroutine TreeRay_potentialListOfBlocks(blockCount,blockList)


  use TreeRay_data, ONLY : tr_boundary, &
       tr_bhUseTreeRay, tr_comm, tr_useGravity
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       GRID_PDE_BND_ISOLATED, GRID_PDE_BND_DIRICHLET, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_solvePoisson
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList


  real, POINTER, DIMENSION(:,:,:,:) :: solnVec

  integer       :: ierr

  real          :: invscale
  integer       :: bcTypes(6)
  real          :: bcValues(2,6) = 0.


!=========================================================================

  if (.not.tr_bhUseTreeRay) return

  ! if gravity is used, 
  if (tr_useGravity) return


  invscale = 1.

  bcTypes = tr_boundary
  where (bcTypes == PERIODIC)
     bcTypes = GRID_PDE_BND_PERIODIC
  elsewhere (bcTypes == ISOLATED)
     bcTypes = GRID_PDE_BND_ISOLATED
  elsewhere (bcTypes == DIRICHLET)
     bcTypes = GRID_PDE_BND_DIRICHLET
  elsewhere (bcTypes == OUTFLOW)
     bcTypes = GRID_PDE_BND_NEUMANN
  end where
  bcValues = 0.

  call Timers_start("treeray Barrier")
  call MPI_Barrier (tr_comm, ierr)
  call Timers_stop("treeray Barrier")

  call Timers_start("treeray")

  ! 0 is used for ipotvar instead of GPOT_VAR
  ! GPOT_VAR is not defined if TreeRay is used without Gravity
  call Grid_solvePoisson (0, DENS_VAR, bcTypes, bcValues, &
       invscale)

  call Timers_stop ("treeray")
  
  return
end subroutine TreeRay_potentialListOfBlocks

