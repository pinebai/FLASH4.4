!!****if* source/Grid/GridSolvers/Multigrid/gr_hgRestrict
!!
!! NAME
!!  gr_hgRestrict
!!
!! SYNOPSIS
!!
!!  gr_hgRestrict(integer, intent(in) :: level,
!!                integer, intent(in) :: ito,
!!                integer, intent(in) :: ifrom)
!!
!! DESCRIPTION
!!  
!!  Restrict the interior data from blocks on a given level to their parents.
!!  This is used for defining the gravity source and residual on all levels.
!!  It works by zeroing the parent blocks, and then averaging 2^DIM child
!!  cells into a single parent cell, handling on and off-processor blocks
!!  separately.
!!
!! ARGUMENTS
!!
!!  level - the child level being restricted from
!!  ito   - the grid variable to restrict into (on the parents)
!!  ifrom - the grid variable to restrict from (on the children)
!!
!! NOTES
!!  
!!  This routine is significantly faster than the alternatives like
!!  gr_restrictTree
!!
!! SEE ALSO
!!
!! gr_hgProlong, gr_restrictTree
!!
!!***

!==============================================================================

!!REORDER(5): unk
!!REORDER(4): solnData

subroutine gr_hgRestrict(level, ito, ifrom)

  use Grid_data, ONLY : gr_meshComm, gr_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
  use tree, ONLY : lrefine, lnblocks,child,nchild,parent,nodetype
  use gr_hgdata, ONLY: send_restrict_data, recv_restrict_data, &
                       send_restrict_req, nbbuf_restrict,&
                       hg_restrict_n1, hg_restrict_n2, hg_restrict_n3,&
                       hg_ili, hg_iui, hg_jli, hg_jui, hg_kli, hg_kui
  use physicaldata, ONLY: unk

  implicit none

#include "constants.h"
#include "Flash.h"  
#include "Flash_mpi.h"
  
  integer, intent(in)          :: ifrom, ito, level
  
  integer                      :: b, c, p, h, i, j, k
  integer                      :: ierr1, ierr2, ierr3
  integer                      :: i1, i2, j1, j2, k1, k2, ichild
  integer                      :: ierr, nsent
  real, pointer                :: solnData(:,:,:,:)
  integer                      :: status(MPI_STATUS_SIZE)
  integer                      :: send_status(MPI_STATUS_SIZE,nbbuf_restrict)
  
  !=======================================================================
  !call timer_start("gr_hgRestrict")


  nsent = 0
  h = 1
  do b = 1, lnblocks
     if (lrefine(b) == level) then

        ! First, compute restricted values for this block.
        
        if (NDIM == 1) then
           
           do i = 1, hg_restrict_n1
              i1 = hg_ili + 2*i - 2
              i2 = i1 + 1
              send_restrict_data(i,1,1,h)  = 0.5 * (unk(ifrom,i1,1,1,b) + &
                   unk(ifrom,i2,1,1,b))
           enddo
           
        else if (NDIM == 2) then
           
           do j = 1, hg_restrict_n2
              j1 = hg_jli + 2*j - 2
              j2 = j1 + 1
              do i = 1, hg_restrict_n1
                 i1 = hg_ili + 2*i - 2
                 i2 = i1 + 1
                 send_restrict_data(i,j,1,h) = 0.25*(unk(ifrom,i1,j1,1,b)+&
                      unk(ifrom,i2,j1,1,b)+&
                      unk(ifrom,i1,j2,1,b)+&
                      unk(ifrom,i2,j2,1,b))
              enddo
           enddo
           
        else ! NDIM == 3
           
           do k = 1, hg_restrict_n3
              k1 = hg_kli + 2*k - 2
              k2 = k1 + 1
              do j = 1, hg_restrict_n2
                 j1 = hg_jli + 2*j - 2
                 j2 = j1 + 1
                 do i = 1, hg_restrict_n1
                    i1 = hg_ili + 2*i - 2
                    i2 = i1 + 1
                    send_restrict_data(i,j,k,h) = 0.125*(unk(ifrom,i1,j1,k1,b)+&
                         unk(ifrom,i2,j1,k1,b)+&
                         unk(ifrom,i1,j2,k1,b)+&
                         unk(ifrom,i2,j2,k1,b)+&
                         unk(ifrom,i1,j1,k2,b)+&
                         unk(ifrom,i2,j1,k2,b)+&
                         unk(ifrom,i1,j2,k2,b)+&
                         unk(ifrom,i2,j2,k2,b))
                 enddo
              enddo
           enddo
        endif
        
        
        ! Next, if parent is on this processor, copy the restricted values directly.
        ! If parent is off-processor, send the data to the owning processor
        ! (non-blocking).
        if (parent(2,b) == gr_meshMe) then ! local parent
           p = parent(1,b)
           do ichild = 1, nchild
              if (child(1,ichild,p) == b) then
                 c = ichild
                 exit
              endif
              if (ichild == nchild) then
                call Driver_abortFlash("[gr_hgRestrict] could not find child block");
              endif
           enddo
           if ((c == 1) .or. (c == 3) .or. (c == 5) .or. (c == 7)) then
              i1 = hg_ili
              i2 = hg_ili + hg_restrict_n1 - 1
           else
              i1 = hg_iui - hg_restrict_n1 + 1
              i2 = hg_iui
           endif
           if ((c == 1) .or. (c == 2) .or. (c == 5) .or. (c == 6)) then
              j1 = hg_jli
              j2 = hg_jli + hg_restrict_n2 - 1
           else
              j1 = hg_jui - hg_restrict_n2 + 1
              j2 = hg_jui
           endif
           if ((c == 1) .or. (c == 2) .or. (c == 3) .or. (c == 4)) then
              k1 = hg_kli
              k2 = hg_kli + hg_restrict_n3 - 1
           else
              k1 = hg_kui - hg_restrict_n3 + 1
              k2 = hg_kui
           endif
           unk(ito,i1:i2,j1:j2,k1:k2,p) = send_restrict_data(:,:,:,h)
        else                          ! remote parent
           nsent = nsent + 1
           call mpi_issend(send_restrict_data(1,1,1,h), &
                hg_restrict_n1*hg_restrict_n2*hg_restrict_n3, &
                FLASH_REAL, parent(2,b), b, &
                gr_meshComm, send_restrict_req(nsent), ierr)
           h = h + 1
        endif
        
        if (h > nbbuf_restrict) &
             call Driver_abortFlash("Buffer space exceeded in gr_hgRestrict")
        
     endif
  enddo
  
  do b = 1, lnblocks
     if ((lrefine(b) == level-1) .and. (nodetype(b) > LEAF)) then
        do c = 1, nchild
           if (child(2,c,b) /= gr_meshMe) then
              
              ! If child is on another processor, receive the restricted data from
              ! the child (blocking).
              
              call mpi_recv(recv_restrict_data(1,1,1), &
                   hg_restrict_n1*hg_restrict_n2*hg_restrict_n3, &
                   FLASH_REAL, child(2,c,b), child(1,c,b), &
                   gr_meshComm, status, ierr)
              if ((c == 1) .or. (c == 3) .or. (c == 5) .or. (c == 7)) then
                 i1 = hg_ili
                 i2 = hg_ili + hg_restrict_n1 - 1
              else
                 i1 = hg_iui - hg_restrict_n1 + 1
                 i2 = hg_iui
              endif
              if ((c == 1) .or. (c == 2) .or. (c == 5) .or. (c == 6)) then
                 j1 = hg_jli
                 j2 = hg_jli + hg_restrict_n2 - 1
              else
                 j1 = hg_jui - hg_restrict_n2 + 1
                 j2 = hg_jui
              endif
              if ((c == 1) .or. (c == 2) .or. (c == 3) .or. (c == 4)) then
                 k1 = hg_kli
                 k2 = hg_kli + hg_restrict_n3 - 1
              else
                 k1 = hg_kui - hg_restrict_n3 + 1
                 k2 = hg_kui
              endif
              
              unk(ito,i1:i2,j1:j2,k1:k2,b) = recv_restrict_data(:,:,:)
              
           endif
        enddo
     endif
  enddo
  
  call mpi_waitall(nsent, send_restrict_req, send_status, ierr)
  
  !call timer_stop("gr_hgRestrict")
  
  !==============================================================================
  
  return
end subroutine gr_hgRestrict
