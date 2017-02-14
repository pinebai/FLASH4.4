!!****if* source/Grid/GridMain/paramesh/Paramesh2/Grid_getBlkNeighLevels
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
!!                                   <  lrefine(blockID)
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

  use Grid_data, ONLY: gr_meshMe
  use tree, ONLY : nfaces, neigh, neigh_type, lrefine, nodetype, lnblocks

  implicit none

  integer,intent(in)  :: blockID
  integer,intent(OUT) :: levels(-1:1, -K2D:K2D , -K3D:K3D)

  integer :: myRefine,myType,iface,i,j,k
  integer :: block1, iface2,i2,j2,k2
  integer :: block2, iface3,i3,j3,k3
! Using PARAMETER instead of SAVE in the next statement sometimes leads to internal
! compiler errors with GNU Fortran (GCC) 4.1.2 20080704 (Red Hat 4.1.2-54) - KW
  integer, save, dimension(3,6) :: &
       nei2surr = &
       RESHAPE( (/ -1,0,0 , 1,0,0 , 0,-1,0 , 0,1,0 , 0,0,-1 , 0,0,1 /), &
                SHAPE=(/3,6/) )

  myRefine = lrefine(blockID)
  myType = nodetype(blockID)

  levels(:,:,:) = 0
  levels(0,0,0) = myRefine

  if (myType .NE. LEAF) then

     do iface=1,nfaces
        i=nei2surr(1,iface)
        j=nei2surr(2,iface)
        k=nei2surr(3,iface)
        if (neigh(1,iface,blockID) == -1) then
           levels(i,j,k) = myRefine - 1
        else if (neigh(1,iface,blockID) <= PARAMESH_PHYSICAL_BOUNDARY) then
           levels(i,j,k) = myRefine
        else if (neigh_type(iface,blockID) == LEAF) then
           levels(i,j,k) = myRefine
        else
           levels(i,j,k) = 10002
        end if
     end do

  else
     
     do iface=1,nfaces
        i=nei2surr(1,iface)
        j=nei2surr(2,iface)
        k=nei2surr(3,iface)
        if (neigh(1,iface,blockID) == -1) then
           levels(i,j,k) = myRefine - 1
        else if (neigh(1,iface,blockID) <= PARAMESH_PHYSICAL_BOUNDARY) then
           levels(i,j,k) = myRefine
        else if (neigh_type(iface,blockID) == 2) then
           levels(i,j,k) = myRefine + 1
        else
           levels(i,j,k) = myRefine
        end if
        if (NDIM > 1 .AND. &
            neigh(1,iface,blockID) > 0  .AND. & 
            neigh(1,iface,blockID) .LE. lnblocks  .AND. & 
            neigh(1,iface,blockID) .NE. blockID  .AND. & 
            neigh(2,iface,blockID) == gr_meshMe) then
           block1 = neigh(1,iface,blockID)
           do iface2=1,nfaces
              if (.NOT. ALL(abs(nei2surr(:,iface2))==abs(nei2surr(:,iface)))) then
                 i2=nei2surr(1,iface2)
                 j2=nei2surr(2,iface2)
                 k2=nei2surr(3,iface2)
                 if (levels(i+i2,j+j2,k+k2) == 0) then
                    if (neigh(1,iface2,block1) == -1) then
                       levels(i+i2,j+j2,k+k2) = myRefine - 1
                    else if (neigh(1,iface2,block1) <= PARAMESH_PHYSICAL_BOUNDARY) then
                       levels(i+i2,j+j2,k+k2) = myRefine
                    else if (neigh_type(iface2,block1) == 2) then
                       levels(i+i2,j+j2,k+k2) = myRefine + 1
                    else
                       levels(i+i2,j+j2,k+k2) = myRefine
                    end if
                    if (NDIM > 2 .AND. &
                         neigh(1,iface2,block1) > 0  .AND. & 
                         neigh(1,iface2,block1) .LE. lnblocks  .AND. & 
                         neigh(1,iface2,block1) .NE. blockID  .AND. & 
                         neigh(1,iface2,block1) .NE. block1  .AND. & 
                         neigh(2,iface2,block1) == gr_meshMe) then
                       block2 = neigh(1,iface2,block1)
                       do iface3=1,nfaces
                          i3=nei2surr(1,iface3)
                          j3=nei2surr(2,iface3)
                          k3=nei2surr(3,iface3)
                          if (abs(product( (/i+i2+i3,j+j2+j3,k+k2+k3/) )) == 1) then
                             if (levels(i+i2+i3,j+j2+j3,k+k2+k3) == 0) then
                                if (neigh(1,iface3,block2) == -1) then
                                   levels(i+i2+i3,j+j2+j3,k+k2+k3) = myRefine - 1
                                else if (neigh(1,iface3,block2) <= PARAMESH_PHYSICAL_BOUNDARY) then
                                   levels(i+i2+i3,j+j2+j3,k+k2+k3) = myRefine
                                else if (neigh_type(iface3,block2) == 2) then
                                   levels(i+i2+i3,j+j2+j3,k+k2+k3) = myRefine + 1
                                else
                                   levels(i+i2+i3,j+j2+j3,k+k2+k3) = myRefine
                                end if
                             end if
                          end if
                       end do
                    end if
                 end if
              end if
           end do
        end if
     end do


  end if

end subroutine Grid_getBlkNeighLevels
