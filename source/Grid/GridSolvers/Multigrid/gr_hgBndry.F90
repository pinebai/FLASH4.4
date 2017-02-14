!!****if* source/Grid/GridSolvers/Multigrid/gr_hgBndry
!!
!! NAME
!!  gr_hgBndry
!!
!! SYNOPSIS
!!
!!  gr_hgBndry(integer, intent(in) :: level,
!!             integer, intent(in) :: ivar,
!!             integer, intent(IN) :: nlayers,
!!             integer, intent(in) :: leafOnly,
!!             integer, intent(in) :: iopt,
!!             integer, intent(in) :: call,
!!             logical, intent(in) :: extrap)
!!
!! DESCRIPTION
!! 
!!  Copy a given variable in and out of work, updating boundary zones.
!!  Level-setting used to be handled here, but it's now in gr_hgSetMaxLevel.
!!  This routine copies in-or-out the blocks on level (0 for all levels), and
!!  fills the guardcells.
!!
!! ARGUMENTS
!!
!!  level     - the level to copy, or all levels if 0
!!  ivar      - the variable to operate on
!!  nlayers   - minimum number of guard cell layers to fill.
!!  leafOnly  - only operate on the leaf blocks. Only used when level=0.
!!              Underlying calls may act on all blocks anyway.
!!  iopt      - this can be
!!              MG_COPY_UNK_TO_WORK copy unk to work, exchange guardcells
!!              MG_UPDATE_UNK       do as above, but copy back as well
!!              MG_EXCHANGE_WORK    just exchange guardcells; no copy
!!  call      - this can be
!!              MG_BEGIN_SERIES first call in a series
!!              MG_CONTINUE_SERIES intermediate call in a series
!!              MG_END_SERIES cleanup (no guardcell fill)
!!              MG_STANDALONE do everything
!!  extrap    - passed to gr_hgGuardCell.
!!                .false. => fill the boundary as if the function goes to zero at the edge
!!                .true.  => extrapolate the values outwards instead
!!
!! NOTES
!!
!!  In this version of the algorithm this routine is less important in how
!!  it is structured.  Copying in and out routinely is probably a good idea
!!  for coherent dataflow, but turning off the copies is a reasonable 
!!  optimization.
!!
!!  Generally gr_hgSolveLevel calls this routine with nlayers = 1.  The number
!!  of guard cell layers actually filled by PARAMESH will be adjusted (upward)
!!  in gr_hgGuardCell, as required by monotonic interpolation routines etc.
!!
!! SEE ALSO
!!
!!  gr_hgGuardCell
!!***

!!REORDER(5): unk
!!REORDER(4): solnData

subroutine gr_hgBndry (level, ivar, nlayers, leafOnly, iopt, call, extrap)

  use physicaldata, ONLY : unk
  use workspace, ONLY : work
  use tree, ONLY : nodetype, lrefine, lnblocks, newchild

  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_blklist
  use gr_hgData, ONLY: gr_hgsaveNodetype, gr_hgSaveNewchild,&
       hg_ili, hg_jli, hg_kli, hg_iui, hg_jui, hg_kui, &
       hg_ile, hg_jle, hg_kle, hg_iue, hg_jue, hg_kue, &
       gr_hgSolnIndex

  use Driver_interface, ONLY : Driver_getElapsedWCTime
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_hgInterface, ONLY: gr_hgGuardCell
  use Grid_interface, ONLY: Grid_getBlkPtr, Grid_releaseBlkPtr

  implicit none

#include "Multigrid.h"
#include "Flash.h"
#ifdef Grid_releaseBlkPtr
! disabling drift per-block logging, see: drift
#undef Grid_releaseBlkPtr
#endif

#include "constants.h"

  integer, intent(IN) :: level, ivar, nlayers, leafOnly, iopt, call
  logical, intent(IN) :: extrap
  real,pointer :: solnData(:,:,:,:)
  integer, save :: myPE, numPEs, numProcs
  integer :: i, j, k, lb



  call Timers_start("gr_hgBndry")
  myPE=gr_meshMe
  numProcs = gr_meshNumProcs
  if ((call == MG_BEGIN_SERIES) .or. (call == MG_CONTINUE_SERIES) .or. &
    (call == MG_STANDALONE)) then
   
     call Timers_start("work copy")

     if ((iopt == MG_COPY_UNK_TO_WORK) .or. (iopt == MG_UPDATE_UNK)) then
        do lb = 1, lnblocks
          call Grid_getBlkPtr(lb,solnData)
           if (((level == 0) .and. (leafOnly==0 .OR. nodetype(lb) == LEAF)) .or. &
                (lrefine(lb) == level) .or. &
                (lrefine(lb) == level+1) .or. (lrefine(lb) == level-1)) then
              do k = hg_kli, hg_kui
                 do j = hg_jli, hg_jui
                    do i = hg_ili, hg_iui
                       work(i,j,k,lb,1) = solnData(ivar,i,j,k)
                    enddo
                 enddo
              enddo
           endif
          call Grid_releaseBlkPtr(lb, solnData)
        enddo
     endif
     
     call Timers_stop("work copy")
     
  endif
  
  call gr_hgGuardCell(myPE, nlayers, extrap)
  
  call Timers_start("work copy")
  
  if (iopt == MG_UPDATE_UNK) then
     do lb = 1, lnblocks
        if (((level == 0) .and. (leafOnly==0 .OR. nodetype(lb) == LEAF)) .or. &
             (lrefine(lb) == level) .or. &
             (lrefine(lb) == level+1) .or. (lrefine(lb) == level-1)) then
           unk(ivar,hg_ile:hg_ili-1,:,:,lb) = work(hg_ile:hg_ili-1,:,:,lb,1)
           unk(ivar,hg_iui+1:hg_iue,:,:,lb) = work(hg_iui+1:hg_iue,:,:,lb,1)
           if (NDIM >= 2) then
              unk(ivar,:,hg_jle:hg_jli-1,:,lb) = work(:,hg_jle:hg_jli-1,:,lb,1)
              unk(ivar,:,hg_jui+1:hg_jue,:,lb) = work(:,hg_jui+1:hg_jue,:,lb,1)
           endif
           if (NDIM == 3) then
              unk(ivar,:,:,hg_kle:hg_kli-1,lb) = work(:,:,hg_kle:hg_kli-1,lb,1)
              unk(ivar,:,:,hg_kui+1:hg_kue,lb) = work(:,:,hg_kui+1:hg_kue,lb,1)
           endif
        endif
     enddo
  endif

   call Timers_stop("work copy")
   call Timers_stop("gr_hgBndry")
!================================================================

  return
end subroutine gr_hgBndry
