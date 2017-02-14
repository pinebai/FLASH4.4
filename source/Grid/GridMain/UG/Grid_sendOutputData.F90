!!****if* source/Grid/GridMain/UG/Grid_sendOutputData
!!
!! NAME
!!  Grid_sendOutputData
!!
!! SYNOPSIS
!!
!!  call Grid_sendOutputData()
!!  
!! DESCRIPTION 
!!   This routine prepares the Grid_IO data that is needed 
!!   in order to write the data to a checkpoint or plotfile.
!!   Data includes, ngid (nfaces + nchildren + 1 parent), globalOffset
!!   globalNumBlocks and other data that is specific to the grid type, 
!!   Paramesh, UG or other
!!
!!  ARGUMENTS  
!!
!!***

subroutine Grid_sendOutputData()

  use Grid_data, ONLY : gr_axisNumProcs, gr_str_geometry, &
       gr_globalNumBlocks, gr_gid, gr_globalOffset, gr_meshMe, gr_meshNumProcs
  use IO_interface, ONLY : IO_setScalar

  implicit none
#include "Flash.h"


  integer :: leftNeigh, rightNeigh, blockID

  !Put NXB, NYB and NZB into a saved variable to prevent problems with
  !an xlf "feature."  These don't mean a whole lot in nofbs mode.
  integer,parameter :: local_nxb = NXB
  integer,parameter :: local_nyb = NYB
  integer,parameter :: local_nzb = NZB
  integer,parameter :: dimensionality = NDIM

  gr_globalNumBlocks = gr_meshNumProcs

  gr_globalOffset = mod(gr_meshMe, gr_globalNumBlocks)



  !set the scalars for the grid unit
  call IO_setScalar("nxb", local_nxb)
  call IO_setScalar("nyb", local_nyb)
  call IO_setScalar("nzb", local_nzb)
  call IO_setScalar("dimensionality", dimensionality)

  
  call IO_setScalar("globalNumBlocks", gr_globalNumBlocks)
  
  call IO_setScalar("geometry", gr_str_geometry)

  
  !!Write the global ID
  
  !-------------------------------------------------------------------------
  ! compute the global id -- this is a single array which stores the 
  ! neighbor block numbers, the parent, and the children of a given block
  !-------------------------------------------------------------------------
  

  !! get the neighbor blocks - the flash3 UG way is to just to number
  !! them in fortran ordering , gid(ngid, blockID)
  !! in UG there is always one block per proc so blockID = 1
  
  blockID = 1
  
  
  if (NDIM >= 1) then
     leftNeigh = gr_meshMe - 1 
     rightNeigh = gr_meshMe + 1
     !! gr_axisNumProcs is the number of procs in the x dim     
     if (mod(gr_meshMe, gr_axisNumProcs(1)) == 0) then
        leftNeigh = -21
     endif
     
     if (mod(gr_meshMe + 1, gr_axisNumProcs(1)) == 0) then
        rightNeigh = -21
     endif
     
          
     gr_gid(1, blockID) = leftNeigh     
     gr_gid(2, blockID) = rightNeigh
     !! non existent for parent and children
     gr_gid(3, blockID) = -21
     gr_gid(4, blockID) = -21
     gr_gid(5, blockID) = -21
     
  end if
  
  
#if NDIM >= 2
  
  leftNeigh = gr_meshMe + gr_axisNumProcs(1) 
  rightNeigh = gr_meshMe - gr_axisNumProcs(1) 
  
  if(leftNeigh >= gr_meshNumProcs) then
     leftNeigh = -21
  end if
  
  if(rightNeigh < 0) then
     rightNeigh = -21
  end if
  
  gr_gid(3, blockID) = leftNeigh   
  gr_gid(4, blockID) = rightNeigh
  
  !! set parents and children to non existent
  gr_gid(5, blockID) = -21  
  gr_gid(6, blockID) = -21
  gr_gid(7, blockID) = -21
  gr_gid(8, blockID) = -21
  gr_gid(9, blockID) = -21
  
#endif

#if NDIM > 2
  
  leftNeigh = gr_meshMe + (gr_axisNumProcs(1) * gr_axisNumProcs(2)) 
  rightNeigh = gr_meshMe - (gr_axisNumProcs(1) * gr_axisNumProcs(2)) 
  
  if(leftNeigh >= gr_meshNumProcs) then
     leftNeigh = -21
  end if
  
  if(rightNeigh < 0) then
     rightNeigh = -21
  end if
  
  gr_gid(5, blockID) = leftNeigh   
  gr_gid(6, blockID) = rightNeigh
  
  !! set parents and children to non existent
  gr_gid(7, blockID) = -21  
  gr_gid(8, blockID) = -21
  gr_gid(9, blockID) = -21
  gr_gid(10, blockID) = -21
  gr_gid(11, blockID) = -21
  gr_gid(12, blockID) = -21
  gr_gid(13, blockID) = -21
  gr_gid(14, blockID) = -21
  gr_gid(15, blockID) = -21
  
#endif

  





end subroutine Grid_sendOutputData
