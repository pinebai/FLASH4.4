!!****if* source/Grid/GridMain/paramesh/Grid_getBlkNeighLevels
!!
!! NAME
!!
!!  Grid_getBlkNeighLevels
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkNeighLevels(integer(IN)  :: blockID,
!!                              integer(OUT) :: levels(-1:1,-K2D:K2D,-K3D:K3D) )
!!
!! DESCRIPTION
!!
!!  For a given block, return refinement level information for its
!!  neigboring blocks in all directions, including diagonal.
!!
!!  It is assumed that blockID indicates a leaf block.
!!
!!  The level information returned can be used to determine the leaf
!!  block resolution from which data in the various guard cell regions
!!  of block blockID is derived when Grid_fillGuardCells is called.
!!  In particular, the caller can determine from this info for each
!!  guard cell region whether the guard cell data has been derived from
!!  interpolation or averaging, or by copying of unmodified leaf block
!!  data only.
!!
!! ARGUMENTS
!!
!!   blockID : ID of block in current processor
!!   levels : returns refinement level for the neigboring regions if
!!            the information is easily available on the executing MPI task,
!!            or one of the following indicators:
!!            * 1 .. lrefine_max     a refinement level
!!            * -2                   unknown or mixed refinement level
!!                                   <= lrefine(blockID)
!!            * -3                   unknown or mixed refinement level
!!            * 10002                unknown or mixed refinement level
!!                                   >= lrefine(blockID)
!!            * 10003                unknown or mixed refinement level
!!                                   >  lrefine(blockID)
!!            * 0                    no information available
!!            The argument levels is an integer array of shape (3,1,1) for
!!            NDIM=1, (3,3,1) for NDIM=2, or (3,3,3) for NDIM=3.
!!            If its bounds are declared as dimension(-1:1,-K2D:K2D,-K3D:K3D),
!!            then
!!            * levels(0,0,0) always refers to the inner cells of block blockID
!!              itself and will always be set to the refinement level of blockID.
!!              The functionality of Grid_getBlkNeighLevels thus includes that
!!              of Grid_getBlkRefineLevel in all implementations.
!!            * levels(-1,0,0) and (+1,0,0) refer to left (West) and right (East)
!!              neighbors in the I direction, respectively.
!!            * levels(0,-1,0) and (0,+1,0) refer to lower (North) and right (South)
!!              neighbors in the J direction, respectively.
!!            * levels(0,0,-1) and (0,0,+1) refer to above and below
!!              neighbors in the K direction, respectively.
!!            * other elements refer to neighbors in diagonal directions.
!!
!! NOTES
!!
!!  With a PARAMESH 4 Grid implementation, refinement level information
!!  is taken from the PARAMESH private array surr_blks and is always
!!  available for all directions.
!!  With the PARAMESH 2 Grid implementation, refinement level information
!!  is taken from the PARAMESH private array neigh and is currently only
!!  available for face directions, not diagonal directions.
!!  With a uniform Grid implementation, refinement level 1 is returned
!!  for all directions since all blocks are considered to be at refinement
!!  level 1.
!!  With the Chombo AMR Grid, this interface is not yet implemented as of FLASH 4.2.
!!
!!  The actual argument corresponding to the level dummy argument may be
!!  declared in the caller in various ways, for example
!!     integer :: levels(-1:1, -K2D:K2D , -K3D:K3D)
!!  or
!!     integer :: levels(3   , 1:1+2*K2D, 1:1+2*K3D)
!!  can be used equivalently (independent of whether NDIM = 1,2, or 3).
!!
!!  NDIM,K2D,K3D are defined in Flash.h or constants.h.
!!
!! SEE ALSO
!!
!!  Grid_getBlkRefineLevel
!!  gr_findAllNeghID
!!***

subroutine Grid_getBlkNeighLevels(blockID, levels)

#include "Flash.h"
#include "constants.h"

  use tree, ONLY : nfaces, neigh, lrefine

  implicit none

  integer,intent(in)  :: blockID
  integer,intent(OUT) :: levels(-1:1, -K2D:K2D , -K3D:K3D)

  integer :: myRefine,iface,i,j,k
  integer,parameter,dimension(3,6) :: &
       nei2surr = &
       RESHAPE( (/ -1,0,0 , 1,0,0 , 0,-1,0 , 0,1,0 , 0,0,-1 , 0,0,1 /), &
                SHAPE=(/3,6/) )

  myRefine = lrefine(blockID)

  levels(:,:,:) = 0
  levels(0,0,0) = myRefine

  do iface=1,nfaces
     i=nei2surr(1,iface)
     j=nei2surr(2,iface)
     k=nei2surr(3,iface)
     if (neigh(1,iface,blockID) == -1) then
        levels(i,j,k) = myRefine - 1
     else if (neigh(1,iface,blockID) <= PARAMESH_PHYSICAL_BOUNDARY) then
        levels(i,j,k) = myRefine
     else
        levels(i,j,k) = 10002
     end if
  end do

end subroutine Grid_getBlkNeighLevels
