!!****if* source/Grid/GridMain/Chombo/AMR/gr_updateRefinement
!!
!! NAME
!!
!!  gr_updateRefinement
!!
!!
!! SYNOPSIS
!!  
!!  call gr_updateRefinement(OPTIONAL, logical(OUT) :: gridChanged)
!!
!!
!! DESCRIPTION
!!
!!  Internal routine to Grid_updateRefinement which handles much of the 
!!  housekeeping for the routine.
!!
!!  This routine calls Chombo to actually
!!  carry out the refinements.  During this stage, the blocks are
!!  redistributed across processors (if needed).
!!
!!  After the refinement, the newly created child blocks are filled via
!!  prolongation from the coarse parents.
!!
!!  Finally, the
!!  EOS is called on the block interiors if this may be necessary to make them
!!  thermodynamically consistent.
!!
!!  Note that this implementation does not fill guardcells after prolongation;
!!  the contents of guard cells should be considered undefined when this
!!  routine returns. The EOS is called on the block interiors if this may be
!!  necessary to make them thermodynamically consistent.
!!
!!  This routine also calls gr_updateParticleRefinement to move the particles
!!  to the correct block after the grid refines.
!!
!! ARGUMENTS
!!
!!  gridChanged : indicates if the grid changed as a result of this call
!!
!!
!! NOTES
!!   This routine was broken off from Grid_updateRefinement to allow the user 
!!   to call this sequence of routines without code duplication.  The user can
!!   now refine/derefine the grid in any user defined way (by using Grid_markRefineDerefine)
!!   etc. and then call gr_updateRefinement to do the housekeeping.
!!
!!   This implementation does not take care of updating guard cells, either
!!   before or after (de)refinement and prolongation. The caller is responsible
!!   for having guard cells filled. At least some layers of valid guard cells
!!   are needed on entry, in order for prolongation to work correctly!
!!
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_updateRefinement(gridChanged)

  use Grid_data, ONLY : gr_blkList, gr_convertToConsvdForMeshCalls, &
       gr_convertToConsvdInMeshInterp, gr_eosMode
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkIndexLimits,&
       Grid_getLocalNumBlks
  use Eos_interface, ONLY : Eos_wrapped
  use Particles_interface, ONLY : Particles_updateRefinement
  use chombo_f_c_interface, ONLY : ch_regrid
  use iso_c_binding, ONLY : c_int
  implicit none

  logical, intent(out),OPTIONAL :: gridChanged

  integer,dimension(2,MDIM) :: blkLimitsGC, blklimits
  integer, dimension(MDIM) :: layers
  integer(c_int) :: baseLevel
  integer :: i, ivar, blkID, count, oldLocalNumBlocks
  logical :: grid_changed, writeBlockInfo

#ifdef DEBUG_CHOMBO
  writeBlockInfo = .true.
#else
  writeBlockInfo = .false.
#endif

  baseLevel = 1 !Unit-based.  Conversion happens in chombo_f_c_api.
  layers = 0
  grid_changed = .true.


  call Timers_start("tree") !2 of 2 (split into 2 so valid to TAU)

  !store the local number of blocks before refinement
  call Grid_getLocalNumBlks(oldLocalNumBlocks)


  ! If using the old logic for conserved variables, convert primitive
  ! variables to conserved in all blocks.
  if (gr_convertToConsvdForMeshCalls) then
     call Grid_getListOfBlocks(ALL_BLKS, gr_blkList, count)
     call gr_primitiveToConserve(gr_blkList, count)
  endif
  

  !DEV CD: Add a grid_changed argument to ch_regrid.
  !Not essential for now as this is only an optimization.
  call ch_regrid(baseLevel)

  if (writeBlockInfo) then
     call gr_writeBlockInfo()
  end if
  


  if (gr_convertToConsvdInMeshInterp) then
     call Grid_getListOfBlocks(ALL_BLKS, gr_blkList, count)
     call gr_sanitizeDataAfterInterp(gr_blkList, count, 'after amr_prolong', layers)
  end if


  ! If using conserved variables (old logic), convert blocks
  ! back from conserved form now.
  if (gr_convertToConsvdForMeshCalls) then
     call Grid_getListOfBlocks(ALL_BLKS, gr_blkList, count)
     call gr_conserveToPrimitive(gr_blkList, count, .FALSE.)
  endif


  !! Do not fill guard cells here any more. All physics units as well
  !! as infrastructure units are required to fill guardcells as they
  !! need them.  That includes particle IO (where guard cells may be needed
  !! for the interpolation stencil for particle properties).

  call Timers_stop("tree") !2 of 2

  ! Call the EOS to make sure the energy and pressure are consistent in the
  ! interiors of blocks.  We only need to do this for
  !  (a)  new leaf (child) blocks that have just been created,
  !  (b)  blocks that were parents but are now leaf blocks because their children
  !       have just been removed.
  ! There is no simple way to check block by block whether these conditions are
  ! true. However, if the grid has not changed at all (i.e., if grid_changed
  ! has not been set to 1 by PARAMESH in the amr_refine_derefine call), we can
  ! be sure that there are no blocks where (a) or (b) applies.

  if (grid_changed) then
     call Timers_start("eos")
     call Grid_getListOfBlocks(ALL_BLKS, gr_blkList, count)

     do i = 1, count
        blkID = gr_blkList(i)
        call Grid_getBlkIndexLimits(blkID, blkLimits, blkLimitsGC)
        call Eos_wrapped(gr_eosMode, blkLimits, blkID)
     end do
     call Timers_stop("eos")
  end if

  
  ! Make sure the particles get moved to the correct block.
  ! If particles are not included this will simply be a stub (empty) routine.
  
  
  call Timers_start("updateParticleRefinement")
  call Particles_updateRefinement(oldLocalNumBlocks)
  call Timers_stop("updateParticleRefinement")
  
  if (present(gridChanged)) gridChanged = grid_changed
  
  return
end subroutine gr_updateRefinement
