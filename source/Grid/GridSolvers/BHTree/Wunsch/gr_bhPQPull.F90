!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhPQPull
!!
!! NAME
!!
!!  gr_bhPQPull
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhMAC(
!!              integer(out)  :: tr,
!!              integer(out)  :: cpu,
!!              integer(out)  :: btp,
!!              integer(out)  :: int_type
!!              real(out)     :: dr(MDIM+2),
!!              real(out)     :: perr
!!             )
!!
!! DESCRIPTION
!!
!!  Pulls (returns and deletes) the first element (the one with the highest
!!  error) from the priority queue.
!!
!! ARGUMENTS
!!
!!  tr          : block index of the pulled tree node
!!  cpu         : cpu of the pulled tree node
!!  btp         : position of the node of the pulled tree node
!!  int_type    : interaction type (see constants in gr_bhData)
!!  dr(MDIM+2)  : position vector between the point-of-calculation and the
!!                pulled node
!!  perr        : estimated error of the contribution of the pulled node
!!
!!
!!***

subroutine gr_bhPQPull(tr, cpu, btp, int_type, dr, perr)
  use gr_bhData, ONLY : gr_bhPriorityQueue, gr_bhPQSize, gr_bhTypePQElement, &
    gr_bhPQNull, gr_bhPQSumSquare, gr_bhTreeMyPE
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(OUT)                :: tr, cpu, btp, int_type
  real, dimension(MDIM+2),intent(OUT) :: dr
  real, intent(OUT)                   :: perr
  integer :: i_cur, i_left, i_right, i_larger
  type(gr_bhTypePQElement) :: dummy

  ! pick the first element
  tr       = gr_bhPriorityQueue(1)%tr       
  cpu      = gr_bhPriorityQueue(1)%cpu      
  btp      = gr_bhPriorityQueue(1)%btp      
  int_type = gr_bhPriorityQueue(1)%int_type 
  dr       = gr_bhPriorityQueue(1)%dr       
  perr     = gr_bhPriorityQueue(1)%perr     

  ! subtract the element perr from SumSquare
  gr_bhPQSumSquare = gr_bhPQSumSquare - perr*perr
  !print *, "PQP: subtracted: ", gr_bhTreeMyPE, perr**2

  ! replace it with the last one
  gr_bhPriorityQueue(1) = gr_bhPriorityQueue(gr_bhPQSize)
  gr_bhPriorityQueue(gr_bhPQSize) = gr_bhPQNull
  gr_bhPQSize = gr_bhPQSize - 1

  ! Max-heapify the queue
  i_cur = 1
  do
    ! get indeces of the two children of i_cur
    i_left  = 2 * i_cur
    i_right = i_left + 1
    if (i_right > gr_bhPQSize) then
      if (i_left > gr_bhPQSize) then
        exit
      else
        i_larger = i_left
      endif
    else
      ! find the larger one
      if (gr_bhPriorityQueue(i_left)%perr > gr_bhPriorityQueue(i_right)%perr) then
        i_larger = i_left
      else
        i_larger = i_right
      endif
    endif
    
    if (i_larger <= gr_bhPQSize) then
      if (gr_bhPriorityQueue(i_cur)%perr < gr_bhPriorityQueue(i_larger)%perr) then
        dummy = gr_bhPriorityQueue(i_larger)
        gr_bhPriorityQueue(i_larger) = gr_bhPriorityQueue(i_cur)
        gr_bhPriorityQueue(i_cur) = dummy
        i_cur = i_larger
      else
        ! correct order => exit
        exit
      endif
    else
      ! no children => exit
      exit
    endif
  enddo

  return
end subroutine gr_bhPQPull



