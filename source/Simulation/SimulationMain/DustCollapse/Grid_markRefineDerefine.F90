!!****if* source/Simulation/SimulationMain/DustCollapse/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!!  This routine is normally called by the implementation of
!!  Grid_updateRefinement.
!!
!! ARGUMENTS
!!
!!  none
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

!!DEV: Should there be a REORDER(...) solnData here? - KW
!!DEV: Yes, there should. - PR
!!REORDER(4): solnData

subroutine Grid_markRefineDerefine()

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount,&
                        gr_enforceMaxRefinement, gr_maxRefine,&
                        gr_lrefineMaxByTime,&
                        gr_lrefineMaxRedDoByTime,&
                        gr_lrefineMaxRedDoByLogR,&
                        gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK,&
                        gr_blkList, gr_eosModeNow, gr_eosMode, &
                        gr_meshMe, gr_meshComm

  use tree, ONLY : newchild, refine, derefine, stay, nodetype,&
       lrefine,lrefine_max, parent, nchild,child
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Grid_interface, ONLY : Grid_fillGuardCells
  use Particles_interface, only: Particles_sinkMarkRefineDerefine

  use Simulation_data, ONLY : sim_initDens, sim_ictr,sim_jctr,&
       sim_kctr, sim_initRad
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getListOfBlocks,&
                             Grid_releaseBlkPtr,Grid_fillGuardCells
  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"

  
  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,iref,blkCount,lb,i,j
  logical,save :: gcMaskArgsLogged = .FALSE.
  integer,save :: eosModeLast = 0
  integer,parameter :: maskSize = NUNK_VARS

  real, dimension(:,:,:,:), pointer :: solnData
  real :: maxdens(MAXBLOCKS),maxdens_parent(MAXBLOCKS)
  integer :: nsend,nrecv,ierr

  integer, dimension(MAXBLOCKS) :: reqr
  integer, dimension(MAXBLOCKS*nchild) :: reqs
  integer, dimension(MPI_STATUS_SIZE,MAXBLOCKS) :: statr
  integer, dimension(MPI_STATUS_SIZE,MAXBLOCKS*nchild) :: stats

  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  logical :: gcMask(NUNK_VARS)

  if(gr_lrefineMaxRedDoByTime) then
     call gr_markDerefineByTime()
  end if
  
  if(gr_lrefineMaxByTime) then
     call gr_setMaxRefineByTime()
  end if

  if (gr_eosModeNow .NE. eosModeLast) then
     gcMaskArgsLogged = .FALSE.
     eosModeLast = gr_eosModeNow
  end if


  gcMask = .FALSE.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do
  gcMask(DENS_VAR) = .TRUE.


  if (.NOT.gcMaskArgsLogged) then
     call Logfile_stampVarMask(gcMask, .true., '[Grid_markRefineDerefine]', 'gcArgs')
  end if

  call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.true.,&
       maskSize=NUNK_VARS, mask=gcMask, makeMaskConsistent=.true.,doLogMask=.NOT.gcMaskArgsLogged,&
       selectBlockType=ACTIVE_BLKS)
  gcMaskArgsLogged = .TRUE.

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


!------------------------------------------------------------------------------
!
! Apply problem-specific refinement criteria.

! Dust collapse problem:  refine center of cloud.  _Don't_ refine blocks
! that are in the "fluff" (max density less than 0.5*starting density of cloud).

  call gr_markInRadius(sim_ictr, sim_jctr, sim_kctr, sim_initRad, lrefine_max)

  call Grid_getListOfBlocks(ACTIVE_BLKS, gr_blkList,blkCount)

  do i = 1, blkCount
     lb = gr_blkList(i)
     call Grid_getBlkPtr(lb,solnData,CENTER)
     maxdens(lb) = maxval(solnData(DENS_VAR,:,:,:))
     call Grid_releaseBlkPtr(lb,solnData)
  end do

! Communicate maxdens of parents to their leaf children.
! Maximally refined children collect messages from parents.

  maxdens_parent(:) = 0.0
  nrecv = 0
  do i = 1, blkCount
     lb = gr_blkList(i)
     if (nodetype(lb) == LEAF .AND. lrefine(lb) == lrefine_max) then
        if(parent(1,lb).gt.-1) then
           if (parent(2,lb).ne.gr_meshMe) then
              nrecv = nrecv + 1
              call MPI_IRecv(maxdens_parent(lb),1, &
                   FLASH_REAL, &
                   parent(2,lb), &
                   lb, &
                   gr_meshComm, &
                   reqr(nrecv), &
                   ierr)
           else
              maxdens_parent(lb) = maxdens(parent(1,lb))
           end if
        end if
     end if
  end do

  ! parents send maxdens to children

  nsend = 0
  do i = 1, blkCount
     lb = gr_blkList(i)
     if (nodetype(lb) == PARENT_BLK .AND. lrefine(lb) == lrefine_max-1) then
        do j = 1,nchild
           if(child(1,j,lb).gt.-1) then
              if (child(2,j,lb).ne.gr_meshMe) then
                 nsend = nsend + 1
                 call MPI_ISend(maxdens(lb), &
                      1, &
                      FLASH_REAL, &
                      child(2,j,lb), &  ! PE TO SEND TO
                      child(1,j,lb), &  ! THIS IS THE TAG
                      gr_meshComm, &
                      reqs(nsend), &
                      ierr)
              end if
           end if
        end do
     end if
  end do

  if (nsend.gt.0) then
     call MPI_Waitall (nsend, reqs, stats, ierr)
  end if
  if (nrecv.gt.0) then
     call MPI_Waitall (nrecv, reqr, statr, ierr)
  end if

!!  maxdens_parent(:) = 0.0  ! <-- uncomment line for previous behavior
  do i = 1, blkCount
     lb = gr_blkList(i)
     if (nodetype(lb) == LEAF) then
        if (maxdens(lb) < 0.5*sim_initDens) then
           refine(lb)   = .false.
!           if (maxdens_parent(lb) < 0.5*sim_initDens .AND. .NOT. stay(lb)) derefine(lb)   = .true.
           if (maxdens_parent(lb) < 0.5*sim_initDens) derefine(lb)   = .true.
        endif
     else if (nodetype(lb) == PARENT_BLK) then
        if (maxdens(lb) < 0.5*sim_initDens)  refine(lb)   = .false.
     end if
  enddo
  

  if(gr_enforceMaxRefinement) call gr_enforceMaxRefine(gr_maxRefine)

  if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()

  if(gr_lrefineMaxRedDoByLogR) &
       call gr_unmarkRefineByLogRadius(gr_lrefineCenterI,&
       gr_lrefineCenterJ,gr_lrefineCenterK)

  call Particles_sinkMarkRefineDerefine()


end subroutine Grid_markRefineDerefine
