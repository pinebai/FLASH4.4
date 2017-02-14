!!****if* source/Grid/GridMain/Chombo/AMR/gr_expandDomain
!!
!!  NAME
!!     gr_expandDomain
!!
!!  SYNOPSIS
!!     call gr_expandDomain(logical,(OUT) :: particlesInitialized)
!!
!!  DESCRIPTION
!!
!!    The grid is initialized in gr_createDomain, with specified
!!    number of blocks. Typically single block. This routine
!!    refines appropriate portions of the initialized physical 
!!    domain according to the given criterion, and applies initial 
!!    conditions to the AMR domain. 
!!
!!  ARGUMENTS
!!    particlesInitialized : is true if this routine initialized particles positions
!!
!!***

!!REORDER(4):solnData

#include "constants.h"
#include "flash_bool.h"
#include "Flash.h"
#include "Eos.h"

subroutine gr_expandDomain (particlesInitialized)
  use iso_c_binding, ONLY : c_int
  use chombo_f_c_interface, ONLY : ch_is_initial_refinement_done, &
       ch_refine_initial_grid, ch_finalize_initial_grid
  use Grid_interface, ONLY : Grid_markRefineDerefine
  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs
  implicit none
  include 'Flash_mpi.h'
  logical, intent(out) :: particlesInitialized
  integer(c_int) :: isRefinementDone

  particlesInitialized = .false.


  !The code in gr_expandDomain is laid out in the same way as the grid
  !intialization in Chombo's AMR::initialGrid.  It is not possible to use
  !this AMR class because FLASH needs to set initial conditions
  !(Simulation_initBlock) and perform custom tagging (Grid_markRefineDerefine).
  !The remaining functionality (grid generation and refinement) in AMR class
  !has been copied into our custom Chombo_adaptive_grid object.
  !The logic that tells us which levels to initialize and whether
  !to continue refining is complicated, so it is hidden inside
  !Chombo_adaptive_grid and queried by gr_expandDomain.

  call ch_is_initial_refinement_done(isRefinementDone)

  do while (isRefinementDone == FLASH_FALSE)
     !loop over all levels and initialize grids (1),
     !then loop over all levels and initialize data (2),
     !then loop over levels and generate tags (3)
     !do this in three separate steps to handle
     !the case where initial data is generated
     !using a multilevel operation (for instance,
     !computing initial velocity from initial
     !vorticity through a multilevel elliptic solve)
     !DFM(11/28/2000)

     !(1) and (2)
     call gr_initializeExistingLevels(gr_meshMe)

     !(3)
     !Grid_markRefineDerefine corresponds to AMRLevel::tagCellsInit.
     !IMPORTANT! We must exchange guardcells before tagging!!!!
     call Grid_markRefineDerefine()

     !ch_refine_initial_grid takes the tags and refines / derefines the grid.
     call ch_refine_initial_grid()

     call ch_is_initial_refinement_done(isRefinementDone)
  end do

  !(1) and (2)
  call gr_initializeExistingLevels(gr_meshMe)

  !Performs any post initialization cleanup.  We may need to perform 
  !a coarse average so that coarser levels contain interpolated data
  !from finer levels where possible.
  call ch_finalize_initial_grid()

end subroutine gr_expandDomain


subroutine gr_initializeExistingLevels(myPE)
  use chombo_f_c_interface, ONLY : ch_build_initial_grid
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getListOfBlocks, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped
  use Simulation_interface, ONLY : Simulation_initBlock
  use Grid_data, ONLY : gr_eosModeInit
  use RadTrans_interface, ONLY: RadTrans_sumEnergy
  implicit none
  integer, intent(IN) :: myPE

  real, pointer:: solnData(:,:,:,:)
  integer ,dimension(MAXBLOCKS) :: blkList
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: blkCount, i

  !--------------------------------------------------------------------
  !(1) ch_build_initial_grid wraps calls to AMRLevel::initialGrid.
  !AMRLevel::initialGrid load balances the boxes at this level,
  !allocates the FArrayBox solution arrays and initializes the
  !LevelFluxRegister, FineInterp, CoarseAverage objects for each level.
  !On the first entry it uses the boxes defined in gr_createDomain
  !and on subsequent entries it uses the soon to be refined boxes.
  call ch_build_initial_grid()
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !(2) The code here corresponds to AMRLevel::initialData.
  call Grid_getListOfBlocks(ALL_BLKS, blkList, blkCount)

  do i = 1, blkCount
     !  We need to zero data in case we reuse blocks from previous levels
     !  but don't initialize all data in Simulation_initBlock... in particular
     !  the total vs. internal energies can cause problems in the eos call that 
     !  follows.
     call Grid_getBlkPtr(blkList(i), solnData)
     solnData = 0.0
     call Grid_releaseBlkPtr(blkList(i), solnData)
     !      Now reinitialize the solution on the new grid so that it has
     !      the exact solution.
     call Simulation_initBlock (blkList(i))
  end do

#ifdef ERAD_VAR
     ! Sum radiation energy density over all meshes. This call is
     ! needed for mesh replication.
     call RadTrans_sumEnergy(ERAD_VAR, blkCount, blkList)
#endif

  ! This is here for safety, in case the user did not take care to make things
  ! thermodynamically consistent in the initial state.- KW
  call Timers_start("eos")
  do i = 1, blkCount
     call Grid_getBlkIndexLimits(blkList(i), blkLimits, blkLimitsGC)
     call Eos_wrapped(gr_eosModeInit,blkLimits,blkList(i))
  end do
  call Timers_stop("eos")
  !--------------------------------------------------------------------
end subroutine gr_initializeExistingLevels
