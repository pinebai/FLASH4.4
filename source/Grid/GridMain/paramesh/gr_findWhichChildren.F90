!!****if* source/Grid/GridMain/paramesh/gr_findWhichChildren
!!
!! NAME
!!
!!  gr_findWhichChildren
!!
!!
!! SYNOPSIS
!!
!!  gr_findWhichChildren(integer(IN)  :: numNegh,
!!                       integer(IN)  :: Negh(MDIM),
!!                       integer(OUT) :: whichChildren(numNegh)
!!
!! DESCRIPTION
!!   
!!    Given the face/edge/corner of a neighboring block, this routine
!!    finds out which of the children of this block share the
!!    face/edge/corner with the neighbor block. This routine makes use
!!    of the child numbering as provided by Paramesh, which is as follows
!!
!!         Front face                      back face
!!       ********************           ********************
!!       *        *         *           *        *         *
!!       *   3    *    4    *           *   7    *    8    *
!!       *        *         *           *        *         *
!!       ********************           ********************
!!       *        *         *           *        *         *
!!       *   1    *    2    *           *   5    *    6    *
!!       *        *         *           *        *         *
!!       ********************           ********************
!!
!! ARGUMENTS
!!
!!   numNegh - The number of children that are neighbors, for a 3D problem
!!             there will 4 per face, 2 per edge and 1 per corner. For a 
!!             2D problem, 2 per edge and 1 per corner
!!   Negh    - the definition of face/edge/corner. It is an MDIM sized
!!             1D array, with valid values in first NDIM entries. Those
!!             can be LEFT_EDGE, RIGHT_EDGE or CENTER.
!!   whichChildren - id of the children blocks of interest, as defined by
!!             Paramesh.
!!
!! EXAMPLE
!!
!!    Consider two blocks A and B in a 2 dimensional situation, with the
!!    child identifiers being 1, 2, 3, and 4 as shown below.
!!
!!            block# 1        block# 2
!!       *******************************
!!       *       *      *              *
!!       *   3   *  4   *              *
!!       *       *      *              *
!!       ****************              *
!!       *       *      *              *
!!       *   1   *  2   *              *
!!       *       *      *              *
!!       *******************************
!!
!!   In this example, we are interested in finding out the children of
!!   interest along the Left edge of Block 2. The call should be
!!   made with numNegh=2, Negh(IAXIS)=LEFT_EDGE, Negh(JAXIS)=CENTER.
!!   The returned values in the whichChildren will be 2 and 4.
!!       
!!
!!***
subroutine gr_findWhichChildren(numNegh,Negh,whichChildren)

#include "constants.h"
#include "Flash.h"

  implicit none

  integer,intent(IN) :: numNegh
  integer, dimension(MDIM),intent(IN) :: Negh
  integer, intent(OUT) :: whichChildren(numNegh)
  integer,dimension(NDIM,LOW:HIGH) :: pos
  integer :: i,j,k,n,m
  
  pos=0
  do i = 1, NDIM
     
     !! If the negh is on left edge, 
     !!that is the right half of the neighbor's parent block
     if(Negh(i)==LEFT_EDGE)pos(i,1)=1
     
     !! If the negh is on right edge, 
     !! that is the left half of the neighbor's parent block
     if(Negh(i)==RIGHT_EDGE)pos(i,1)=-1
     
     !! if center, we are interested in both halves.
     !! n here is keeping track of count of the children.
     if(Negh(i)==CENTER)then
        pos(i,LOW)=-1; pos(i,HIGH)=1
     end if
  end do
  
  !! Now translate those positions into the child places within the parent
  !! block
  whichChildren(:)=1
  k=1
  m=2
  do i = 1,NDIM
     if(pos(i,LOW)>=0)whichChildren(:)=whichChildren(:)+k
     if(pos(i,HIGH)>0) then
        if(m==2) then
           whichChildren(m:numNegh)=whichChildren(m:numNegh)+k
        else
           whichChildren(m:numNegh)=whichChildren(LOW:HIGH)+k
        end if
        m=m+1
     end if
     k=k*2
  end do
  return
end subroutine gr_findWhichChildren

