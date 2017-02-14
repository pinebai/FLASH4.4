!!****if* source/Grid/GridSolvers/Multigrid/gr_hgResidual
!!
!! NAME
!!  gr_hgResidual
!!
!! SYNOPSIS
!!  
!!  gr_hgResidual(integer, intent(in) :: level
!!                integer, intent(in) :: gr_iSource
!!                integer, intent(in) :: gr_iSoln
!!                integer, intent(in) :: ires)
!!
!! DESCRIPTION
!!  
!!  Compute the current residual for the solution to the Poisson equation.
!!  This routine uses the work array and boundary filling on the work array
!!  as the boundary values are used in the computation of the residual one
!!  block in.  Only the leaf block residuals are calculated.
!!
!!  This is done by applying the first-order laplacian stencil as defined 
!!  by hg_cx, hg_cy, and hg_cz in gr_hgData and gr_hgInit.
!!
!! ARGUMENTS
!!
!!  level        - the level (of leaf blocks) to take the residual for
!!  gr_iSource - the density variable
!!  gr_iSoln   - the potential variable
!!  ires         - the residual variable
!!
!! RESULT
!!
!!  The residual between the source and solutions at the blocks at level 
!!  is placed the variable ires.
!!
!! NOTES
!!
!!  The original description was somewhat dishonest in that extrapolations
!!  are done for the exterior, however the work boundary fill is quite real.
!!  A potential optimization would be to do this ONCE on all leaf blocks,
!!  saving us a number of boundary fills.
!!
!!  Note that gr_iSource is copied into work on the leaf blocks here.
!!
!!***

!!REORDER(5): unk

subroutine gr_hgResidual(level, gr_iSource, gr_iSoln, ires, dt,chi,theta)

!==============================================================================
#include "Flash.h"
#include "Multigrid.h"
#include "constants.h"

  use tree, ONLY : lnblocks,lrefine_min,lrefine,bsize
  use Grid_interface, ONLY : Grid_getDeltas
  use physicaldata, ONLY : unk
  use workspace, ONLY : work
  use gr_hgData, ONLY : hg_ili, hg_iui, hg_jli, hg_jui, hg_kli, hg_kui, &
                        gr_hgBndTypes, hg_cx, hg_cy, hg_cz, & 
                        gr_hgSaveNodetype

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_hgInterface, ONLY: gr_hgBndry
  use Grid_data, ONLY : gr_meshComm

  implicit none

  include 'Flash_mpi.h'

  integer, intent(in)          :: level, gr_iSource, gr_iSoln, ires
  
  integer                      :: b, i, j, k, ii, jj, kk, ierr

  real, dimension(MDIM)        :: deltas
  real                         :: avg, sum, lsum, vol, lvol, nbinv, bvol, cvol, bsum
  real                         :: dxinv2, dyinv2, dzinv2
  real,intent(IN),OPTIONAL     :: dt, chi, theta
  
  !=======================================================================

  call Timers_start("gr_hgResidual")
  
  
  
  !call mg_bndry(level, gr_iSoln, 1, 0, COPY_UNK_TO_WORK, STANDALONE)
  ! Use (EXCHANGE_WORK, CONTINUE_SERIES) under the assumption that we're
  ! calling gr_hgResidual immediately after a call to hg_solveLevel --
  ! performance savings by not filling work() again
  call gr_hgBndry(level, gr_iSoln, NGUARD, 0, &
                MG_COPY_UNK_TO_WORK, MG_STANDALONE, .false.)  !haven't checked this
  
  do b = 1, lnblocks
     if (lrefine(b) == level) then
                
        do k = NGUARD*K3D+1, NGUARD*K3D+NZB           ! working on interior only
           do j = NGUARD*K2D+1, NGUARD*K2D+NYB
              do i = NGUARD+1, NGUARD+NXB
                 unk(ires,i,j,k,b) = unk(gr_iSource, i,j,k,b)
              enddo
           enddo
        enddo
        
        call Grid_getDeltas(b,deltas)
        
        if (NDIM == 1) then
           
           dxinv2 = 1.0/deltas(IAXIS)**2
           do i = NGUARD+1, NGUARD+NXB
              do ii = -1, 1  ! over the stencil width
                 ! Do the del-squared operator on work, which has been filled with gr_iSoln
                 unk(ires,i,1,1,b) = &
                      unk(ires,i,1,1,b) - dxinv2*hg_cx(ii,i-NGUARD)*work(i+ii,1,1,b,1)
              enddo
           enddo
           
        else if (NDIM == 2) then
           
           dxinv2 = 1.0/deltas(IAXIS)**2
           dyinv2 = 1.0/deltas(JAXIS)**2
           
           do j = NGUARD+1, NGUARD+NYB
              do i = NGUARD+1, NGUARD+NXB
                 do ii = -1, 1
                    unk(ires,i,j,1,b) = &
                         unk(ires,i,j,1,b) - dxinv2*hg_cx(ii,1)*work(i+ii,j,1,b,1)
                 enddo
                do jj = -1, 1
                    unk(ires,i,j,1,b) = &
                         unk(ires,i,j,1,b) - dyinv2*hg_cy(jj,1)*work(i,j+jj,1,b,1)
                 enddo
              enddo
           enddo
           
        else ! NDIM == 3
           
           dxinv2 = 1.0/deltas(IAXIS)**2
           dyinv2 = 1.0/deltas(JAXIS)**2
           dzinv2 = 1.0/deltas(KAXIS)**2
           
           do k = NGUARD+1, NGUARD+NZB
              do j = NGUARD+1, NGUARD+NYB
                 do i = NGUARD+1, NGUARD+NXB
                    do ii = -1, 1
                       unk(ires,i,j,k,b) = &
                            unk(ires,i,j,k,b) &
                            - dxinv2*hg_cx(ii,1)*work(i+ii,j,k,b,1)
                    enddo
                    do jj = -1, 1
                       unk(ires,i,j,k,b) = &
                            unk(ires,i,j,k,b) &
                            - dyinv2*hg_cy(jj,1)*work(i,j+jj,k,b,1)
                    enddo
                    do kk = -1, 1
                       unk(ires,i,j,k,b) = &
                            unk(ires,i,j,k,b) & 
                             - dzinv2*hg_cz(kk,1)*work(i,j,k+kk,b,1)
                    enddo
                 enddo
              enddo
           enddo
        endif
     endif

  enddo

  ! Recenter variable for periodic so average is zero going into fft
  if ( ALL(gr_hgBndTypes(1:2*NDIM) == MG_BND_PERIODIC .or. gr_hgBndTypes(1:2*NDIM) == MG_BND_NEUMANN) ) then
     
     lvol = 0.
     lsum = 0.
     nbinv = 1. / real(NXB)
     if (NDIM >= 2) nbinv = nbinv / real(NYB)
     if (NDIM == 3) nbinv = nbinv / real(NZB)
     
     do b = 1, lnblocks
        if ((lrefine(b) == level) .or. &
             ((gr_hgSaveNodetype(b) == LEAF) .and. (lrefine(b) < level))) then
           bvol = bsize(1,b)
           if (NDIM >= 2) bvol = bvol * bsize(2,b)
           if (NDIM == 3) bvol = bvol * bsize(3,b)
           cvol = bvol * nbinv
           lvol = lvol + bvol
           bsum = 0.
           do k = hg_kli, hg_kui
              do j = hg_jli, hg_jui
                 do i = hg_ili, hg_iui
                    bsum = bsum + unk(ires,i,j,k,b)
                 enddo
              enddo
           enddo
           lsum = lsum + bsum * cvol
        endif
     enddo
     
     call mpi_allreduce ( lsum, sum, 1, FLASH_REAL, &
          MPI_SUM, gr_meshComm, ierr )
     call mpi_allreduce ( lvol, vol, 1, FLASH_REAL, &
          MPI_SUM, gr_meshComm, ierr )
     
     avg = sum / vol
     
     do b = 1, lnblocks
        if (lrefine(b) == level) then
           unk(ires,:,:,:,b) = unk(ires,:,:,:,b) - avg
        endif
     enddo
     
  endif
  
  call Timers_stop("gr_hgResidual")

  !=======================================================================
  
  return
end subroutine gr_hgResidual
