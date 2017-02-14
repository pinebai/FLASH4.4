!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgBndry
!!
!! NAME
!!
!!  gr_mgBndry
!!
!! SYNOPSIS
!!
!!  call gr_mgBndry(integer(in) :: level,
!!                  integer(in) :: ivar,
!!                  integer(in) :: nlayers,
!!                  integer(in) :: leaf_only,
!!                  integer(in) :: iopt,
!!                  integer(in) :: call)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   level : 
!!
!!   ivar : 
!!
!!   nlayers : 
!!
!!   leaf_only : 
!!
!!   iopt : 
!!
!!   call : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_bndry()

!  Description: Update boundary zones on a given level for a given variable.

!  Parameters:  level       Level to update.
!               ivar        Index of variable to update.
!               leaf_only   If nonzero, variable is only defined on leaf-node
!                           blocks; obtain boundary information by restriction
!                           from finer neighbors.
!               iopt        COPY_UNK_TO_WORK: copy unk to work, exchange work's guardcells
!                           EXCHANGE_WORK: just exchange work's guardcells without unk copy
!               call        BEGIN_SERIES: first call in a series
!                           CONTINUE_SERIES: intermediate call in a series
!                           END_SERIES: cleanup call (don't exchange)
!                           CONTINUE_SERIES: do first call, intermediate, and cleanup stuff

subroutine gr_mgBndry(level, ivar, nlayers, leaf_only, iopt, call)

!===============================================================================
#include "Flash.h"

  use gr_mgData, ONLY: nodetype_save, ili, jli, kli, iui, jui, kui, &
                       mesh_lrefmax,solvelevel
  
  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr

  use tree, only : nodetype,lrefine
  use workspace, ONLY: work
  use paramesh_dimensions, only : nguard_work
  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs
  use Driver_data, ONLY: dr_simTime

implicit none
#include "Multigrid.h"
#include "constants.h"


integer, intent(in) :: level, ivar, nlayers, leaf_only, iopt, call

integer :: i, j, k, lb, lnblocks, lrefinei
real    :: time

real,    pointer, dimension(:,:,:,:), save :: solnData

logical, save :: first_call=.true.

integer, parameter :: ndim = NDIM
integer, parameter :: k2d = K2D
integer, parameter :: k3d = K3D


!===============================================================================

!call timer_start("mg_bndry")

!call timer_start("in_mg_bndry_before_gc")

time     = dr_simTime

call Grid_getLocalNumBlks(lnblocks)

if ((call == MG_BEGIN_SERIES).or.(call == MG_CONTINUE_SERIES).or.(call == MG_STANDALONE)) then
!   call timer_start("work copy")
   
   if (iopt.eq.MG_COPY_UNK_TO_WORK) then

      do lb = 1, lnblocks
         lrefinei = lrefine(lb)
         if ((lrefinei == level) .or. &
              (lrefinei == level+1) .or. (lrefinei == level-1)) then

            ! Point to blocks center vars:
            call Grid_getBlkPtr(lb,solnData,CENTER)
            
            do k = kli, kui
               do j = jli, jui
                  do i = ili, iui
                     work(i,j,k,lb,1) = solnData(ivar,i,j,k)
                  enddo
               enddo
            enddo

            ! Release pointers:
            call Grid_releaseBlkPtr(lb,solnData,CENTER)
         endif
      enddo

   endif

!   call timer_stop("work copy")

endif

! Temporarily re-mark the blocks' node types.  The PARAMESH guard cell routine
! works by first restricting leaf-node data to parents, then trading boundary
! data at the same level of refinement, then finally interpolating boundary
! data for nodetype 2 blocks to their leaf-node children.  We alter the
! nodetype and neigh_type arrays (with the help of get_tree_nodetypes) in order
! to fool amr_guardcell into doing what we want.

! What we want is:  if the ivar variable is defined everywhere on a level (not
! just on leaf nodes), we do not want to clobber valid data in the restriction
! step.  So we treat all blocks on the chosen level as leaf nodes, all blocks
! above it as temporarily nonexistent, and blocks below it as leaf or parent,
! depending on what they normally are.  This will clobber valid data in
! 'parent' blocks below the chosen level.

if ((call == MG_BEGIN_SERIES).or.(call == MG_STANDALONE)) then
   if (leaf_only == 0) then    
      
   do lb = 1, lnblocks
      nodetype(lb) = nodetype_save(lb)
      if ( lrefine(lb) > level)  nodetype(lb) = -1
      if ( lrefine(lb) == level) nodetype(lb) = 1
      if ((lrefine(lb) == level-1) .and. (nodetype(lb) /= 1)) &
           nodetype(lb) = 2
   enddo

   call amr_get_new_nodetypes (gr_meshNumProcs, gr_meshMe, level)

   ! If the ivar variable is defined only on leaf nodes, then it's OK to clobber
   ! parent data.  We probably don't need to do anything here (since at most a
   ! jump of one level is supposed to be guaranteed), but we mark all blocks at
   ! the next higher level as leaf-node just to make sure they provide boundary
   ! data to our chosen level.
   
   else
      
   do lb = 1, lnblocks
      nodetype(lb) = nodetype_save(lb)
      if (lrefine(lb) > level+1)  nodetype(lb) = -1
      if (lrefine(lb) == level+1) nodetype(lb) = 1
      if ((lrefine(lb) == level) .and. (nodetype(lb) /= 1)) &
           nodetype(lb) = 2
   enddo

   ! If not Pfft solving at Largest level covering the whole domain
   if ((level .lt. mesh_lrefmax)) &
   call amr_get_new_nodetypes (gr_meshNumProcs, gr_meshMe, level+1)

   endif
end if

!call timer_stop("in_mg_bndry_before_gc")

! Now call the PARAMESH guard cell routine (on the work array).

!call timer_start("mg_guardcell_in_mg_bndry")

if (call == MG_STANDALONE) then
   call gr_mgGuardcell(gr_meshMe, ivar, nlayers, time, 1, 0)
else if ((call == MG_BEGIN_SERIES).or.(call == MG_CONTINUE_SERIES)) then
   call gr_mgGuardcell(gr_meshMe, ivar, nlayers, time, 0, 0)
end if

!call timer_stop("mg_guardcell_in_mg_bndry")
!call timer_stop("mg_bndry")
!===============================================================================

return
end
