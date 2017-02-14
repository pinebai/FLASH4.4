!!****if* source/Simulation/SimulationMain/Jeans/Grid_markRefineDerefine
!!
!! NAME
!! 
!!  Grid_markRefineDerefine
!!
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine()
!!
!! ARGUMENTS
!!
!! DESCRIPTION
!!  This routine is specific to Jeans.  It overrides the default
!!  routine located in Grid/common/paramesh
!!
!!  This routine is the main interface to the refining criteria in FLASH.
!!  By default it will check the second derivative condition for the variables
!!  specified by refine_var_[1-4], and set the refine and derefine logical
!!  flags accordingly.  
!! 
!!  After the second derivatives are checked, it loops through all of the 
!!  blocks and makes sure that they will still have a refinement level in
!!  the range lrefine_min < block's refinement level < lrefine_max.
!!
!!  Any user defined refinement should be put in this routine.  For instance,
!!  to refine at circular region of the grid, you can loop over all the 
!!  blocks, compute the distance from the center of the circle, and set
!!  refine(blockID) = .TRUE. and derefine(blockID) = .FALSE. if the block
!!  falls within the circle.  Much of this functionality is provided by the
!!  MarkGridRef module, whose routines should be called from here.
!!
!!***

!!REORDER(4): solnData

subroutine Grid_markRefineDerefine()

  use Grid_data, ONLY : gr_blkList, gr_numRefineVars,&
                        gr_refineOnParticleCount,&
                        gr_enforceMaxRefinement, gr_maxRefine,&
                        gr_lrefineMaxRedDoByTime,&
                        gr_lrefineMaxRedDoByLogR,&
                        gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK, &
                        gr_meshMe
  use Driver_interface, ONLY : Driver_getSimTime
  use Grid_interface, ONLY : Grid_fillGuardCells, &
    Grid_getListOfBlocks, Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_releaseBlkPtr
  use Grid_data, ONLY : gr_refine_cutoff,gr_derefine_cutoff, &
       gr_refine_filter, gr_numRefineVars, gr_refine_var

  use Tree, ONLY : derefine, refine, lrefine, nodetype, stay, lnblocks
  use Tree, ONLY : newchild, lrefine_min, lrefine_max, child
  use Tree, ONLY : parent, nchild

  use Simulation_data, ONLY: sim_deltaRef, sim_deltaDeRef, sim_refDensity

  use Eos_interface, ONLY : Eos_wrapped

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer :: idir = ALLDIR
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: mode = MODE_DENS_EI

  integer :: nsend, nrecv
  integer :: reqr(MAXBLOCKS), reqs(MAXBLOCKS*nchild)
  integer :: statr(MPI_STATUS_SIZE, MAXBLOCKS)
  integer :: stats(MPI_STATUS_SIZE, MAXBLOCKS*nchild)
  
  integer :: ierr

  integer :: blockCount

  real          :: time

  real, pointer, dimension(:,:,:,:) :: solnData

  integer       :: i, j, k, l, iref, lb
  real          :: ref_cut, deref_cut, ref_filter

  real          :: delta_max(MAXBLOCKS)
  real          :: delta_max_par(MAXBLOCKS)


  if(gr_lrefineMaxRedDoByTime) then
     call gr_markDerefineByTime()
  end if

  call Grid_fillGuardCells( CENTER, idir)
  call Grid_getListOfBlocks(ACTIVE_BLKS,gr_blkList,blockCount)

  do i = 1, blockCount
     call Grid_getBlkIndexLimits(gr_blkList(i),blkLimits,blkLimitsGC)
     call Eos_wrapped(mode,blkLimitsGC,gr_blkList(i))
  end do

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
  end do

  delta_max(1:lnblocks) = -huge(1.0)

  do l = 1, blockCount

     lb = gr_blkList(l)

     call Grid_getBlkPtr(lb,solnData)
     call Grid_getBlkIndexLimits(lb,blkLimits,blkLimitsGC)

     do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
        do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
           do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)
              delta_max(lb) = max(delta_max(lb), &
                   solnData(DENS_VAR,i,j,k)/sim_refDensity-1.)
           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(lb,solnData)

  end do

  delta_max_par(1:lnblocks) = 0.
  nrecv = 0
  do lb = 1,lnblocks
     if(parent(1,lb) .gt. -1) then
        if (parent(2,lb) .ne. gr_meshMe) then
           nrecv = nrecv + 1
           call MPI_IRecv(delta_max_par(lb),    & 
                1, & 
                FLASH_REAL, & 
                parent(2,lb), & 
                lb, & 
                MPI_COMM_WORLD, & 
                reqr(nrecv), & 
                ierr)
        else
           delta_max_par(lb) = delta_max(parent(1,lb))
        end if
     end if
  end do
  
! parents send error to children

  nsend = 0
  do lb = 1,lnblocks
     do j = 1,nchild
        if(child(1,j,lb) .gt. -1) then
           if (child(2,j,lb) .ne. gr_meshMe) then
              nsend = nsend + 1
              call MPI_ISend(delta_max(lb), & 
                   1, & 
                   FLASH_REAL, & 
                   child(2,j,lb), &  ! PE TO SEND TO
                   child(1,j,lb), &  ! THIS IS THE TAG
                   MPI_COMM_WORLD, & 
                   reqs(nsend), & 
                   ierr)
           end if
        end if
     end do
  end do
  
  if (nsend.gt.0) then
     call MPI_Waitall (nsend, reqs, stats, ierr)
  end if
  if (nrecv.gt.0) then
     call MPI_Waitall (nrecv, reqr, statr, ierr)
  end if
  
  do lb = 1,lnblocks
         
     if (nodetype(lb).eq.1) then            

! test for derefinement

        if (.not.refine(lb).and..not.stay(lb) & 
             .and.delta_max(lb) .le. sim_deltaDeRef & 
             .and.delta_max_par(lb) .le. sim_deltaDeRef) then
           derefine(lb) = .TRUE.
        else
           derefine(lb) = .FALSE.
        end if

! test for refinement

        if (delta_max(lb) .gt. sim_deltaRef) then
           derefine(lb) = .FALSE.
           refine(lb) = .TRUE.
        end if

        if (delta_max(lb) .gt. sim_deltaDeRef & 
             .or. delta_max_par(lb) .gt. sim_deltaDeRef) & 
             stay(lb) = .TRUE.
            
     end if
         
  end do
  
  if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()

  !------------------------------------------------------------------------------

  ! Make sure that the refinement level of the blocks is between lrefine_min
  ! and lrefine_max

  do i = 1, lnblocks

     if ((lrefine(i) == lrefine_min) .and. (nodetype(i) == 1)) & 
          derefine(i) = .FALSE.

     if ((lrefine(i) < lrefine_min) .and. (nodetype(i) == 1)) then
        refine(i) = .TRUE.
        derefine(i) = .FALSE.
     end if

     if ((lrefine(i) == lrefine_max) .and. (nodetype(i) == 1)) & 
          refine(i) = .FALSE.

  end do

  if(gr_enforceMaxRefinement) call gr_enforceMaxRefine(gr_maxRefine)

  if(gr_lrefineMaxRedDoByLogR) &
       call gr_unmarkRefineByLogRadius(gr_lrefineCenterI,&
       gr_lrefineCenterJ,gr_lrefineCenterK)
  return
end subroutine Grid_markRefineDerefine
