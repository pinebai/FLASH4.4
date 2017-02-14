!!****if* source/Grid/GridMain/paramesh/Grid_receiveInputData
!!
!! NAME
!!
!!  Grid_receiveInputData
!!
!! SYNOPSIS
!!
!!  call Grid_receiveInputData(integer(IN) :: localNumBlocks,
!!                             integer(IN) :: alnblocks,
!!                             integer(IN) :: xx)
!!
!! DESCRIPTION 
!!
!!  Initializes grid arrays from arrays read by the I/O unit.
!!
!! ARGUMENTS  
!!
!!  localNumBlocks : the number of blocks on my processor.
!!
!!  alnblocks : the approximate number of local blocks on each
!!              processor if we give each processor an equal
!!              number of blocks.  Calculated from
!!              int(globalNumBlocks/meshNumProcs) + 1.
!!
!!  xx : an integer representing a cutoff point.  Processors
!!       less than this value are assigned alnblocks blocks and
!!       processors greater than or equal to this value are
!!       assigned lnblocks-1 blocks.
!!
!!***

subroutine Grid_receiveInputData(localNumBlocks, alnblocks, xx)

#include "constants.h"
#include "Flash.h"

  use tree, ONLY : nfaces, nchild, neigh, parent, child
#ifdef FLASH_GRID_PARAMESH3OR4
  use Grid_data, ONLY : gr_gsurr_blks, gr_is_gsurr_blks_initialized
#endif
#ifdef FLASH_GRID_PARAMESH4DEV_SURR_BLKS_OPTION
  use tree, ONLY : surr_blks
  use physicaldata, ONLY : use_flash_surr_blks_fill, surr_blks_valid
#endif
  use Grid_data, ONLY : gr_gid, gr_globalMe, gr_meshNumProcs

  implicit none
  integer, intent(IN) :: localNumBlocks, alnblocks, xx
  integer :: div, blockID, i, j, k, ngid

  !-----------------------------------------------------------------------------
  ! build the tree information from the gid array and the number of blocks
  ! on each processor
  !-----------------------------------------------------------------------------

  div = xx*alnblocks

  do blockID = 1,localNumBlocks

     ! neighbor data
     ngid = 0
     do j = 1,nfaces
        ngid = ngid + 1

        if (gr_gid(ngid,blockID).gt.0) then

           if (gr_gid(ngid,blockID).le.div) then
              neigh(2,j,blockID) = int((gr_gid(ngid,blockID)-1)/alnblocks)

              if (neigh(2,j,blockID).gt.gr_meshNumProcs-1)  &
                   neigh(2,j,blockID) = gr_meshNumProcs - 1

              neigh(1,j,blockID) = gr_gid(ngid,blockID) -  &
                   (alnblocks*neigh(2,j,blockID))
           else
              neigh(2,j,blockID) = &
                   int((gr_gid(ngid,blockID)-1-div)/(alnblocks-1)) + xx

              if (neigh(2,j,blockID).gt.gr_meshNumProcs-1)  &
                   neigh(2,j,blockID) = gr_meshNumProcs - 1

              neigh(1,j,blockID) = gr_gid(ngid,blockID) - div - &
                   ((alnblocks-1)*(neigh(2,j,blockID)-xx))
           end if
        else
           neigh(1,j,blockID) = gr_gid(ngid,blockID)
           neigh(2,j,blockID) = gr_gid(ngid,blockID)
        end if
     end do

     ! parent data
     ngid = ngid + 1
     if (gr_gid(ngid,blockID).gt.0) then
        if (gr_gid(ngid,blockID).le.div) then
           parent(2,blockID) = int((gr_gid(ngid,blockID)-1)/alnblocks)
           if (parent(2,blockID).gt.gr_meshNumProcs-1)  &
                parent(2,blockID) = gr_meshNumProcs - 1
           parent(1,blockID) = gr_gid(ngid,blockID) -  &
                (alnblocks*parent(2,blockID))
        else
           parent(2,blockID) = &
                int((gr_gid(ngid,blockID)-1-div)/(alnblocks-1)) + xx
           if (parent(2,blockID).gt.gr_meshNumProcs-1)  &
                parent(2,blockID) = gr_meshNumProcs - 1
           parent(1,blockID) = gr_gid(ngid,blockID) - div - &
                ((alnblocks-1)*(parent(2,blockID)-xx))
        end if
     else
        parent(1,blockID) = gr_gid(ngid,blockID)
        parent(2,blockID) = gr_gid(ngid,blockID)
     end if

     ! children data
     do j = 1,nchild
        ngid = ngid + 1
        if (gr_gid(ngid,blockID).gt.0) then
           if (gr_gid(ngid,blockID).le.div) then
              child(2,j,blockID) = int((gr_gid(ngid,blockID)-1)/alnblocks)
              if (child(2,j,blockID).gt.gr_meshNumProcs-1)  &
                   child(2,j,blockID) = gr_meshNumProcs - 1
              child(1,j,blockID) = gr_gid(ngid,blockID) -  &
                   (alnblocks*child(2,j,blockID))
           else
              child(2,j,blockID) = &
                   int((gr_gid(ngid,blockID)-1-div)/(alnblocks-1)) + xx
              if (child(2,j,blockID).gt.gr_meshNumProcs-1)  &
                   child(2,j,blockID) = gr_meshNumProcs - 1
              child(1,j,blockID) = gr_gid(ngid,blockID) - div - &
                   ((alnblocks-1)*(child(2,j,blockID)-xx))
           end if
        else
           child(1,j,blockID) = gr_gid(ngid,blockID)
           child(2,j,blockID) = gr_gid(ngid,blockID)
        end if
     end do

  end do


#ifdef FLASH_GRID_PARAMESH4DEV_SURR_BLKS_OPTION
  ! surrounding neighbor data
  if (gr_is_gsurr_blks_initialized) then
     if (gr_globalMe == MASTER_PE) &
          print *, "read 'gsurr_blks' dataset - applying pm4dev optimization."

     do blockID = 1,localNumBlocks
        do k = 1, 1+(K3D*2)
           do j = 1, 1+(K2D*2)
              do i = 1, 1+(K1D*2)
                 !Store the node type first of all.
                 surr_blks(3,i,j,k,blockID) = &
                      gr_gsurr_blks(2,i,j,k,blockID)
                 if (gr_gsurr_blks(1,i,j,k,blockID).gt.0) then
                    if (gr_gsurr_blks(1,i,j,k,blockID).le.div) then
                       surr_blks(2,i,j,k,blockID) = &
                            int((gr_gsurr_blks(1,i,j,k,blockID)-1)/alnblocks)
                       if (surr_blks(2,i,j,k,blockID).gt.gr_meshNumProcs-1)  &
                            surr_blks(2,i,j,k,blockID) = gr_meshNumProcs - 1
                       surr_blks(1,i,j,k,blockID) = &
                            gr_gsurr_blks(1,i,j,k,blockID) - &
                            (alnblocks*surr_blks(2,i,j,k,blockID))
                    else
                       surr_blks(2,i,j,k,blockID) = &
                            int((gr_gsurr_blks(1,i,j,k,blockID)-1-div)/&
                            (alnblocks-1)) + xx
                       if (surr_blks(2,i,j,k,blockID).gt.gr_meshNumProcs-1)  &
                            surr_blks(2,i,j,k,blockID) = gr_meshNumProcs - 1
                       surr_blks(1,i,j,k,blockID) = &
                            gr_gsurr_blks(1,i,j,k,blockID) - div - &
                            ((alnblocks-1)*(surr_blks(2,i,j,k,blockID)-xx))
                    end if
                 else
                    surr_blks(1,i,j,k,blockID) = gr_gsurr_blks(1,i,j,k,blockID)
                    surr_blks(2,i,j,k,blockID) = gr_gsurr_blks(1,i,j,k,blockID)
                    surr_blks(3,i,j,k,blockID) = gr_gsurr_blks(2,i,j,k,blockID)
                 end if
              end do
           end do
        end do
     end do
     surr_blks_valid = use_flash_surr_blks_fill
  else
     surr_blks_valid = .false.
  end if
#endif

end subroutine Grid_receiveInputData
