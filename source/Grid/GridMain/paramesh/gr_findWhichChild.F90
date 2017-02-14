!!****if* source/Grid/GridMain/paramesh/gr_findWhichChild
!!
!! NAME
!!
!!  gr_findWhichChild
!!
!!
!! SYNOPSIS
!!
!!  gr_findWhichChild(real(IN)     :: pos(MDIM),
!!                    real(IN)     :: bndBox(LOW:HIGH,MDIM),
!!                    integer(IN)  :: negh(MDIM),
!!                    integer(OUT) :: whichChild)
!!
!! DESCRIPTION
!!   
!!    Given a block bounding box and the co-ordinates of a physical point on the 
!!    domain, this routine finds out which of the children of this block
!!    is in the neighbourhood of the point. This function is useful when
!!    a Lagrangian particle is moving out of the current block, to another 
!!    block which is at a higher refinement. In such a sitution, the current
!!    block can find out the parent of the destination block, but still doesn't
!!    know which of that parent's children is the destination. Paramesh 
!!    designates children of a block in a predetermined order. This routine
!!    determines the appropriate child identifier and returns that value
!!
!!
!! ARGUMENTS
!!   pos     - co-ordinates of position of interest
!!   bndBox - bounding box of current block
!!   negh    - For current block, the location of the neighboring parent block
!!   whichChild - identification of the child of the neighboring parent block
!!                that includes the given position
!!
!! EXAMPLE
!!
!!    Consider two blocks A and B in a 2 dimensional situation, with the
!!    child identifiers being 1, 2, 3, and 4 as shown below.
!!
!!            block# 1        block# 2
!!       *******************************
!!       *       *      *              *
!!       *   4   *  3   *              *
!!       *       *    + *              *
!!       ****************              *
!!       *       *      *              *
!!       *   1   *  2   *              *
!!       *       *      *              *
!!       *******************************
!!       
!!   In this example, point of interest is marked by the "+" sign.
!!   The input values will be the coordinates of the point at the +
!!   sign in the argument "pos", blockID=2, negh(IAXIS)=LEFT_EDGE
!!   negh(JAXIS)=CENTER, and the returned value will be 3.
!!
!!***
subroutine gr_findWhichChild(pos,bndBox,negh,whichChild)

#include "constants.h"
#include "Flash.h"
  implicit none

  real,dimension(MDIM), intent(IN) :: pos
  real,dimension(LOW:HIGH,MDIM),intent(IN) :: bndBox
  integer, dimension(MDIM),intent(IN) :: negh
  integer, intent(OUT) :: whichChild
  real,dimension(MDIM) :: delta
  real :: delDiff,midpoint
  integer :: i,k

  delta(1:MDIM) = bndBox(HIGH,1:MDIM)-bndBox(LOW,1:MDIM)
  whichChild=1
  k=1
  do i = 1,NDIM
     delDiff=0.0
     if(negh(i)==LEFT_EDGE)delDiff=-delta(i)/2.0
     if(negh(i)==RIGHT_EDGE)delDiff=3.0*delta(i)/2.0
     if(negh(i)==CENTER)delDiff=delta(i)/2.0
     midPoint=bndBox(LOW,i)+delDiff
     if(pos(i)>=midPoint)whichChild=whichChild+k
     k=k*2
  end do

end subroutine gr_findWhichChild
