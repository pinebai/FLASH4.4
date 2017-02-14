!!****if* source/Grid/GridMain/paramesh/gr_updateData
!!
!! NAME
!!
!!  gr_updateData
!!
!! SYNOPSIS
!!
!!  call gr_updateData()
!!
!! DESCRIPTION
!!
!!
!! This is private subroutine which should be called from the function interface
!! for grid refinement, immediately after the mesh package finished refining. 
!! It should be called only when there is re-gridding going on. For each 
!! block on a processor, this routine gets the bounding box and refinement level
!! information, using these quantities the data structure for each block 
!! is calculated and filled.
!!
!!
!!***


subroutine gr_updateData()

  use tree, ONLY : bnd_box,lnblocks,lrefine,lrefine_max,nodetype
  use Grid_data, ONLY:  gr_oneBlock, gr_delta, gr_imin, gr_jmin,gr_kmin
  use Grid_data, ONLY:  gr_boxContainingLeafNodes

  implicit none
#include "constants.h"
#include "Flash.h"

  integer :: i,j,stride,indLeft,indRight
  real :: indCenter,temp
  integer :: nBlocks
  real :: bnd_box_x,bnd_box_y,bnd_box_z
  real,dimension(NDIM) :: half_delta


  gr_boxContainingLeafNodes(LOW ,IAXIS) =  HUGE(bnd_box_x)
  gr_boxContainingLeafNodes(HIGH,IAXIS) = -HUGE(bnd_box_x)
#if NDIM > 1
  gr_boxContainingLeafNodes(LOW ,JAXIS) =  HUGE(bnd_box_x)
  gr_boxContainingLeafNodes(HIGH,JAXIS) = -HUGE(bnd_box_x)
#else
  gr_boxContainingLeafNodes(LOW ,JAXIS) =  gr_jmin
  gr_boxContainingLeafNodes(HIGH,JAXIS) =  gr_jmin
#endif
#if NDIM > 2
  gr_boxContainingLeafNodes(LOW ,KAXIS) =  HUGE(bnd_box_x)
  gr_boxContainingLeafNodes(HIGH,KAXIS) = -HUGE(bnd_box_x)
#else
  gr_boxContainingLeafNodes(LOW ,KAXIS) =  gr_kmin
  gr_boxContainingLeafNodes(HIGH,KAXIS) =  gr_kmin
#endif

  ! the variable delta stores dx, dy, and dz values for
  ! each level of refinement given an initial domain size
  ! set in Grid_init
  nBlocks = lnblocks
  half_delta = gr_delta(1:NDIM,lrefine_max)/2.0
  do i = 1,nBlocks
     bnd_box_x = bnd_box(1,1,i)
     bnd_box_y = bnd_box(1,2,i)
     bnd_box_z = bnd_box(1,3,i)
     if (nodetype(i) == 1) then ! for LEAF blocks
        gr_boxContainingLeafNodes(LOW ,IAXIS) = min(bnd_box_x,     gr_boxContainingLeafNodes(LOW,IAXIS))
        gr_boxContainingLeafNodes(HIGH,IAXIS) = max(bnd_box(2,1,i),gr_boxContainingLeafNodes(HIGH,IAXIS))
#if NDIM > 1
        gr_boxContainingLeafNodes(LOW ,JAXIS) = min(bnd_box_y,     gr_boxContainingLeafNodes(LOW,JAXIS))
        gr_boxContainingLeafNodes(HIGH,JAXIS) = max(bnd_box(2,2,i),gr_boxContainingLeafNodes(HIGH,JAXIS))
#endif
#if NDIM > 2
        gr_boxContainingLeafNodes(LOW ,KAXIS) = min(bnd_box_z,     gr_boxContainingLeafNodes(LOW,KAXIS))
        gr_boxContainingLeafNodes(HIGH,KAXIS) = max(bnd_box(2,3,i),gr_boxContainingLeafNodes(HIGH,KAXIS))
#endif
     end if
     stride = 2**(lrefine_max - lrefine(i))

     !cornerID is a unique integer value
     !assigned to each block
     gr_oneBlock(i)%cornerID(1)= (bnd_box_x-gr_imin+half_Delta(1))/&
                                  gr_delta(1,lrefine_max) + 1

     ! calculate the left most index, use this to calculate coords
     indLeft = gr_oneBlock(i)%cornerID(1)-stride*NGUARD -1
     do j=1,NXB+2*NGUARD
        indCenter = indLeft + stride/2.0
        indRight = indLeft + stride
        gr_oneBlock(i)%firstAxisCoords(LEFT_EDGE,j) = gr_imin+&
             gr_delta(1,lrefine_max)*indLeft
        gr_oneBlock(i)%firstAxisCoords(CENTER,j) = gr_imin + &
             gr_delta(1,lrefine_max)*indCenter
        gr_oneBlock(i)%firstAxisCoords(RIGHT_EDGE,j) = gr_imin+&
             gr_delta(1,lrefine_max)*indRight
        indLeft = indRight
     end do
#if NDIM > 1
     gr_oneBlock(i)%cornerID(2)= (bnd_box_y-gr_jmin+half_delta(2))/&
                                  gr_delta(2,lrefine_max) + 1
     indLeft = gr_oneBlock(i)%cornerID(2)-stride*NGUARD -1
     do j=1,NYB+2*NGUARD
        indCenter = indLeft + stride/2.0
        indRight = indLeft + stride
        gr_oneBlock(i)%secondAxisCoords(LEFT_EDGE,j) = gr_jmin+&
             gr_delta(2,lrefine_max)*indLeft
        gr_oneBlock(i)%secondAxisCoords(CENTER,j) = gr_jmin+&
             gr_delta(2,lrefine_max)*indCenter
        gr_oneBlock(i)%secondAxisCoords(RIGHT_EDGE,j) = gr_jmin+&
             gr_delta(2,lrefine_max)*indRight
        indLeft = indRight
     end do
#else
     gr_oneBlock(i)%cornerID(2) = 1
     gr_oneBlock(i)%secondAxisCoords(LEFT_EDGE,1) = gr_jmin
     gr_oneBlock(i)%secondAxisCoords(CENTER,1) = gr_jmin
     gr_oneBlock(i)%secondAxisCoords(RIGHT_EDGE,1) = gr_jmin
#endif
#if NDIM > 2
     gr_oneBlock(i)%cornerID(3)= (bnd_box_z-gr_kmin+half_delta(3))/&
                                  gr_delta(3,lrefine_max) + 1
     indLeft = gr_oneBlock(i)%cornerID(3)-stride*NGUARD -1
     do j=1,NZB+2*NGUARD
        indCenter = indLeft + stride/2.0
        indRight = indLeft + stride
        gr_oneBlock(i)%thirdAxisCoords(LEFT_EDGE,j) = gr_kmin+&
             gr_delta(3,lrefine_max)*indLeft
        gr_oneBlock(i)%thirdAxisCoords(CENTER,j) = gr_kmin+&
             gr_delta(3,lrefine_max)*indCenter
        gr_oneBlock(i)%thirdAxisCoords(RIGHT_EDGE,j) = gr_kmin+&
             gr_delta(3,lrefine_max)*indRight
        indLeft = indRight
     end do
#else
     gr_oneBlock(i)%cornerID(3) = 1
     gr_oneBlock(i)%thirdAxisCoords(LEFT_EDGE,1) = gr_kmin
     gr_oneBlock(i)%thirdAxisCoords(CENTER,1) = gr_kmin
     gr_oneBlock(i)%thirdAxisCoords(RIGHT_EDGE,1) = gr_kmin
#endif

  end do
  
end subroutine gr_updateData


