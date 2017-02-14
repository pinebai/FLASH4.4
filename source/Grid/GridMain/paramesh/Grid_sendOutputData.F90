!!****if* source/Grid/GridMain/paramesh/Grid_sendOutputData
!!
!! NAME
!!  Grid_sendOutputData
!!
!! SYNOPSIS
!!
!!  call Grid_sendOutputData()
!!  
!! DESCRIPTION 
!!
!!  This routine allows the Grid unit to checkpoint any scalar data
!!  stored in the Grid_data Fortran modules and is called by the
!!  routine IO_updateScalars before checkpointing.  To send data to
!!  the IO unit this routine calls IO_setScalar.  In addition this
!!  routine may prepare any other data the IO unit needs for
!!  checkpointing which the Grid unit owns.
!!
!!  For example, the Grid unit owns the variable which defines the
!!  grid geometry gr_geometry. This value needs to be checkpointed so
!!  that visualization tools can determine if the simulation ran with
!!  cartesian, spherical, cylindrical etc. coordinates.  To send
!!  scalar data such as gr_geometry to the IO unit to be checkpointed,
!!  the routine calls IO_setScalar("geometry", gr_geometry)
!!
!!  Other data which is sent to the IO unit is globalNumBlocks, globalOffset
!!  and globalNumBlocks
!!
!!  This routine also prepares the gid data structure which is also
!!  needed for visualization purposes. gid((nfaces + nchildren + 1 parent)
!!
!!
!!  ARGUMENTS  
!!
!!
!!  SEE ALSO
!!   
!!   IO_setScalar, IO_updateScalars
!!  
!!***

subroutine Grid_sendOutputData()

#include "constants.h"
#include "Flash.h"

  use Grid_data, ONLY :gr_str_geometry, gr_globalNumBlocks, gr_nToLeft, gr_globalOffset, &
       gr_gid, gr_meshComm, gr_meshMe, gr_meshNumProcs
  use Grid_interface, ONLY : Grid_getLocalNumBlks

  use tree, ONLY : neigh, child, parent, nchild, nfaces
#ifdef FLASH_GRID_PARAMESH3OR4
  use Grid_data, ONLY : gr_gsurr_blks
  use tree, ONLY : surr_blks
#endif

  use IO_interface, ONLY : IO_setScalar

  implicit none

 include "Flash_mpi.h"

  integer :: localNumBlocks, i, ierr, blockID, j, ngid, k

  !Put NXB, NYB and NZB into a saved variable to prevent problems with
  !an xlf "feature."
  integer,parameter :: local_nxb = NXB
  integer,parameter :: local_nyb = NYB
  integer,parameter :: local_nzb = NZB
  integer,parameter :: dimensionality = NDIM


  !set the scalars for the grid unit
  call IO_setScalar("nxb", local_nxb)
  call IO_setScalar("nyb", local_nyb)
  call IO_setScalar("nzb", local_nzb)
  
  call IO_setScalar("geometry", gr_str_geometry)
  call IO_setScalar("dimensionality", dimensionality)
  
  
  
  
  !! Get the local number of blocks from everybody
  !! to find total number of blocks
  call Grid_getLocalNumBlks(localNumBlocks)  

  call MPI_Allgather(localNumBlocks, 1, MPI_INTEGER, gr_nToLeft, & 
       1, MPI_INTEGER, gr_meshComm, ierr)
  
  gr_globalNumBlocks = 0
  
  do i = 0,gr_meshNumProcs-1
     gr_globalNumBlocks = gr_globalNumBlocks + gr_nToLeft(i)         
  end do
     
  
  call IO_setScalar("globalNumBlocks", gr_globalNumBlocks)
  
  ! compute the number of processors to the left of a processor
  do i = gr_meshNumProcs-1,1,-1
     gr_nToLeft(i) = gr_nToLeft(i-1)
  end do
  
  
  gr_nToLeft(0) = 0
  do i = 2,gr_meshNumProcs-1
     gr_nToLeft(i) = gr_nToLeft(i) + gr_nToLeft(i-1)
  end do
  
  gr_globalOffset = gr_nToLeft(gr_meshMe)
     



  !-----------------------------------------------------------------------------
  ! compute the global id -- this is a single array which stores the 
  ! neighbor block numbers, the parent, and the children of a given block
  !-----------------------------------------------------------------------------
  do blockID = 1,localNumBlocks
     ngid = 0

     ! loop over the faces and store the neighbors
     do j = 1,nfaces
        ngid = ngid + 1

        ! if the neighbor exists, then store the block number of the neighbor
        ! -- take into account the number of blocks below the processor 
        ! that the
        ! neighbor is on, so the block number is global
        if (neigh(1,j,blockID).gt.0) then
           gr_gid(ngid,blockID) = neigh(1,j,blockID) +  & 
                gr_nToLeft(neigh(2,j,blockID))
        else

           ! the neighbor is either a physical boundary or does not exist at that 
           ! level of refinement
           gr_gid(ngid,blockID) = neigh(1,j,blockID)
        end if
     end do
        
     ! store the parent of the current block
     ngid = ngid + 1
     if (parent(1,blockID).gt.0) then
        gr_gid(ngid,blockID) = parent(1,blockID) +  & 
             gr_nToLeft(parent(2,blockID))
     else
        gr_gid(ngid,blockID) = parent(1,blockID)
     end if
     
     ! store the children of the current block
     do j = 1,nchild
        ngid = ngid + 1
        if (child(1,j,blockID).gt.0) then
           gr_gid(ngid,blockID) = child(1,j,blockID) +  & 
                gr_nToLeft(child(2,j,blockID))
        else
           gr_gid(ngid,blockID) = child(1,j,blockID)
        end if
     end do

#ifdef FLASH_GRID_PARAMESH3OR4
     ! store the surrounding neighbors of the current block
     do k = 1, 1+(K3D*2)
        do j = 1, 1+(K2D*2)
           do i = 1, 1+(K1D*2)
              if (surr_blks(1,i,j,k,blockID).gt.0) then
                 gr_gsurr_blks(1,i,j,k,blockID) = surr_blks(1,i,j,k,blockID) + &
                      gr_nToLeft(surr_blks(2,i,j,k,blockID))
              else
                 gr_gsurr_blks(1,i,j,k,blockID) = surr_blks(1,i,j,k,blockID)
              end if
              ! store the node type in index 2 (no need to store proc ID)
              gr_gsurr_blks(2,i,j,k,blockID) = surr_blks(3,i,j,k,blockID)
           end do
        end do
     end do
#endif

  end do


  !need to make sure this happens before checkpointing
  call gr_ptFillBlkParticleInfo()


end subroutine Grid_sendOutputData
