!!****if* source/Grid/GridMain/paramesh/gr_updateRefinement
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
!!  This routine calls the PARAMESH routine amr_refine_derefine to actually
!!  carry out the refinements.  During this stage, the blocks are
!!  redistributed across processors (if needed).
!!
!!  After the refinement, the newly created child blocks are filled via
!!  prolongation from the coarse parents.  This prolongation step can use
!!  prolongation routines from paramesh, or defined by the user.
!!
!!  Once the prolongation is done, the guardcells are filled.  Finally, the
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

#define DEBUG_POSITIVITY 0

subroutine gr_updateRefinement( gridChanged)

#include "Flash.h"
  use Grid_data, ONLY : gr_blkList, gr_convertToConsvdForMeshCalls,&
       gr_convertToConsvdInMeshInterp, gr_eosMode, gr_meshMe, gr_gcellsUpToDate
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_getListOfBlocks,Grid_getBlkIndexLimits
  use Simulation_interface, ONLY : Simulation_customizeProlong
  
  use tree, ONLY : newchild, nodetype, lnblocks, grid_changed
#ifndef FLASH_GRID_PARAMESH2
  use physicaldata, ONLY: mpi_pattern_id, no_permanent_guardcells
#endif
  use paramesh_interfaces, ONLY : amr_refine_derefine, &
                                  amr_prolong
  use Eos_interface, ONLY : Eos_wrapped
  use Particles_interface, ONLY : Particles_updateRefinement
  implicit none


#include "constants.h"

  logical, intent(out),OPTIONAL :: gridChanged

  integer :: i, ivar
  integer :: blkListAncestors(MAXBLOCKS)

  integer,dimension(2,MDIM) :: blkLimitsGC, blklimits
  integer :: count, numLeafBlocks, numAncestorBlocks
  integer :: oldLocalNumBlocks !need this if running with particles
  integer, dimension(MDIM) :: layers
!=============================================================================

  call Timers_start("tree") !2 of 2 (split into 2 so valid to TAU)
  !store the local number of blocks before refinement
  oldLocalNumBlocks = lnblocks
  layers =0 
  call Grid_getBlkIndexLimits(1,blkLimits,blkLimitsGC)
  
  grid_changed = 0              ! will be 1 after amr_refine_derefine if the grid actually changed

    
  ! Now perform the indicated refinements and derefinements.  First blocks are
  ! destroyed and new children are created, then blocks are redistributed.

  ! Paramesh2: During redistribution, the block interiors and all
  ! guard cells are moved. (NOTE that in FLASH2 it was possible that
  ! only a single perimeter of guardcells were moved, depending on the
  ! choice of amr_redist_blk implementation. The copy of amr_redist_blk
  ! used in the Paramesh2 code of FLASH3 does copy all guard cells.)

  ! Paramesh 3 and later: During redistribution, only data in block interiors are moved.
  call Timers_start("amr_refine_derefine")
  call amr_refine_derefine()
  call Timers_stop("amr_refine_derefine")
#ifndef FLASH_GRID_PARAMESH2
  if (grid_changed .NE. 0) mpi_pattern_id = -abs(mpi_pattern_id) !make it different from recognized values
#endif           
  
  ! update the grid coordinates for the new mesh
  call Timers_start("updateData")
  call gr_updateData()
  call Timers_stop("updateData")
  
  ! If using the old logic for conserved variables, convert primitive
  ! variables to conserved in all blocks that are not new children.
  ! Include even ancestor blocks, since (w/ PM3) the data of some of them
  ! may be needed for interpolation in amr_prolong - if they are ancestors
  ! that were actually parents up to the amr_refine_derefine call.
  ! The prolonging will stuff the new children.
  if (gr_convertToConsvdForMeshCalls) then
      count = 0
      numAncestorBlocks = 0
      do i = 1, lnblocks
         if ((nodetype(i)==LEAF .AND. .not. newchild(i)) .OR. nodetype(i)==PARENT_BLK) then
            count = count+1
            gr_blkList(count)=i
         else if (nodetype(i)==ANCESTOR) then
            numAncestorBlocks = numAncestorBlocks+1
            blkListAncestors(numAncestorBlocks)=i
         endif
      enddo
     call gr_primitiveToConserve(gr_blkList,count)
     call gr_primitiveToConserve(blkListAncestors,numAncestorBlocks)
  endif
  
  
  ! Initialize the data in the newly created children by prolonging the data
  ! from the parent to the children.

  ! OLD COMMENT: the prolongation stencil can only include a single
  ! coarse zone on either side of the parent of the new children,
  ! since only one guardcell is guaranteed to be valid.
  ! ACTUAL FLASH3 BEHAVIOR:
  ! with Paramesh2: the parents of the new children still have guard
  ! cell values as they were on entry to this subroutine.
  ! with Paramesh3 and later: the values of variables in guard cells of the
  ! parents of new children may have become invalid when blocks were
  ! redistributed in amr_refine_derefine. However, within amr_prolong
  ! processing, amr_1blk_guardcell will be called to fill in guard
  ! cells of such blocks as needed.
  ! In both cases, the prolongation stencil (interpolation stencil)
  ! can make use of all guard cell layers available (i.e., up to NGUARD).
  
  call Simulation_customizeProlong(BEFORE)
  call amr_prolong (gr_meshMe, 1, NGUARD)
  call Simulation_customizeProlong(AFTER)


#ifndef FLASH_GRID_PARAMESH2
  ! the following condition should not be true for PM2 anyway
  if (gr_convertToConsvdInMeshInterp) then
     call Grid_getListOfBlocks(ACTIVE_BLKS, gr_blkList, numLeafBlocks)
     call gr_sanitizeDataAfterInterp(gr_blkList, numLeafBlocks, 'after amr_prolong', layers)
  end if

  if (.NOT. no_permanent_guardcells) then
     call gr_freeCommRecvBuffer
  end if
#endif

  ! If using conserved variables (old logic), convert parent and leaf blocks
  ! back from conserved form now.
  ! We do this for all blocks, even ancestors, to keep data in a sane range
  ! even for blocks whose data won't have any impact on further propagation
  ! (but may still be written out to plotfiles etc.).
  if (gr_convertToConsvdForMeshCalls) then
     call Grid_getListOfBlocks(ACTIVE_BLKS, gr_blkList, numLeafBlocks)
     call gr_conserveToPrimitive(gr_blkList,numLeafBlocks, .FALSE.)
     call gr_conserveToPrimitive(blkListAncestors,numAncestorBlocks, .TRUE.)
  endif

#ifndef FLASH_GRID_PARAMESH2
  if (grid_changed .NE. 0) call gr_checkGridConsistency()
#endif

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

  if (grid_changed .NE. 0) then
     call Timers_start("eos")
     call Grid_getListOfBlocks(LEAF, gr_blkList,count)

     do i = 1, count
          call Eos_wrapped(gr_eosMode ,blkLimits, gr_blkList(i))
     end do
     call Timers_stop("eos")
  end if

  
  ! Make sure the particles get moved to the correct block.
  ! If particles are not included this will simply be a stub (empty) routine.
  
  
  call Timers_start("updateParticleRefinement")
  call Particles_updateRefinement(oldLocalNumBlocks)
  call Timers_stop("updateParticleRefinement")
  
  if (present(gridChanged)) gridChanged = (grid_changed .NE. 0)

  if (gr_gcellsUpToDate) then
     !The guard cells will no longer be up to date if the grid has changed.
     if (grid_changed .NE. 0) then
        gr_gcellsUpToDate = .false.
     end if
  end if

  return
end subroutine gr_updateRefinement
