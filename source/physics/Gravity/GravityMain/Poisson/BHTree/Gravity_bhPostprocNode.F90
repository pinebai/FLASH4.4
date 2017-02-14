!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhPostprocNode
!!
!! NAME
!!
!!  Gravity_bhPostprocNode
!!
!!
!! SYNOPSIS
!!
!!  call Gravity_bhPostprocNode(
!!             real(in)                       :: ndSize,
!!             real(inout)                    :: node(:)
!!             )
!!
!! DESCRIPTION
!!
!!  Called after the tree is built for each node. With certain MACs,
!!  calculates the minimum distance needed for this node to be accepted.
!!
!! ARGUMENTS
!!
!!  ndSize      : physical size of the node (the largest extent of the node)
!!  node        : array of the tree node which is processed
!!
!!
!!***

subroutine Gravity_bhPostprocNode(ndSize, node)
  use Gravity_data, ONLY : useGravity, grv_bhIDMIN, grv_bhNODE5, grv_bhN5_DMIN, &
    grv_bhIB2 , grv_bhAccErr, grv_bhNewton, grv_bhMAC_APE, grv_bhMAC_MPE, &
    grv_bhMPDegree, grv_bhMPDegree_p2, grv_bhMACNum, grv_bhIM
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, intent(IN) :: ndSize
  real, dimension(:), intent(INOUT)  :: node

  if (.not. useGravity) return

  ! convert the second order moment B_2 into the minimum distance
  ! for the node acceptance given by formula 9 in S&W94
  if (grv_bhNODE5 == grv_bhN5_DMIN) then

    select case(grv_bhMACNum)
      case(grv_bhMAC_APE) ! "Gadget" MAC

        node(grv_bhIDMIN) = (grv_bhNewton*node(grv_bhIM) * ndSize**grv_bhMPDegree &
        &                / grv_bhAccErr)**(1.0/grv_bhMPDegree_p2)

      case(grv_bhMAC_MPE) ! Maximum Partial Error

       !print *, "PPN1: ", ndSize, node
       node(grv_bhIDMIN) = 0.5*ndSize + sqrt(0.25*ndSize*ndSize &
       & + sqrt(3.0*grv_bhNewton*node(grv_bhIB2)/grv_bhAccErr))
       !print *, "PPN2: ", ndSize, node
      case default
        ! shouldn't get here at all
        call Driver_abortFlash("Gravity_bhPostProcNode: never should get here.")
    end select

  endif

  return
end subroutine Gravity_bhPostprocNode

