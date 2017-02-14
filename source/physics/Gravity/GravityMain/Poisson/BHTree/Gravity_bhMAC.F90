!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhMAC
!!
!! NAME
!!
!!  Gravity_bhMAC
!!
!!
!! SYNOPSIS
!!
!!   logical res = Gravity_bhMAC(
!!                          real(in)       :: node(:),
!!                          real(in)       :: ndSize2,
!!                          real(in)       :: dr(MDIM+2),
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Multipole Acceptance Criterion. Determines whether the contribution of the
!!  node to the potential at the point of calculation will have a nacessary
!!  accuracy.
!!
!! ARGUMENTS
!!
!!  node        : array of the node tested
!!  ndSize2     : square of the physical size of the node (the largest 
!!                extent of the node)
!!  dr          : (1:MDIM) - position vector from the point-of-calculation to the node
!!                (MDIM+1) - square of the magnitude of the position vector
!!                (MDIM+2) - inverted magnitude of the position vector
!!  blockno     : number of block into which the point-of-calculation belongs
!!  point       : indeces of the point-of-calculation in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!! RESULT
!!  Returns TRUE if the node is accepted for the calculation of the potential.
!!  Otherwise returns FALSE.
!!
!!
!!***

logical function Gravity_bhMAC(node, ndSize2, dr, blockno, point, blkLimits, solnData)
  use Gravity_data, ONLY : grv_bhNewton, grv_bhMeanBlockAccErrInv, grv_bhMACNum &
    & , grv_bhMPDegree, grv_bhMPDegree_p2, grv_bhMPDegree_half &
    & , grv_bhMPDegree_halfp1, grv_bhMAC_APE, grv_bhMAC_MPE, grv_bhMAC_SS &
    & , grv_bhAccErrInv, grv_bhNODE5, grv_bhUseRelAccErr, grv_bhIDMIN &
    & , grv_bhN5_DMIN, grv_bhIM
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize2
  integer, dimension(MDIM), intent(IN) :: point
  real, dimension(MDIM+2), intent(IN) :: dr
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real :: dist_lim2, dist_mpdp2, dist_lim_mpdp2, gacc_err_i, la
  logical :: node_accepted

  !print *, "Gravity MAC: ", node, ndSize2, point, dr
  node_accepted = .false.
  ! if minimum node distance is present in the tree use it
  if (grv_bhNODE5 == grv_bhN5_DMIN) then
    !print *, "MAC: ", dr(MDIM+1), node(grv_bhIDMIN)**2
    if (dr(MDIM+1) > node(grv_bhIDMIN)*node(grv_bhIDMIN)) then
      node_accepted = .true.
    endif
  else
    ! with relative acceleration error control get the inverted value 
    ! of the acceleration at the point/block where the potential is calculated
    if (grv_bhUseRelAccErr) then
      if (point(IAXIS) == -1) then
        ! MAC evaluated for a block (point(IAXIS) == -1)
        gacc_err_i = grv_bhMeanBlockAccErrInv
      else
        gacc_err_i = solnData(ACEI_VAR, point(IAXIS), point(JAXIS), point(KAXIS))
      endif
    else
      gacc_err_i = grv_bhAccErrInv
    endif


    select case(grv_bhMACNum)
 
      case(grv_bhMAC_APE) ! "Gadget" MAC
        ! calculate limitting distance to the power MPDegree + 2
        ! and node size to the power MPDegree
        if (MOD(grv_bhMPDegree, 2) == 1) then
          dist_lim_mpdp2 = grv_bhNewton*node(grv_bhIM) &
          &              * sqrt(ndSize2)**grv_bhMPDegree * gacc_err_i
          dist_mpdp2 = sqrt(dr(MDIM+1))**grv_bhMPDegree_p2
        else
          dist_lim_mpdp2 = grv_bhNewton*node(grv_bhIM) &
          &              * ndSize2**grv_bhMPDegree_half * gacc_err_i
          dist_mpdp2 = dr(MDIM+1)**grv_bhMPDegree_halfp1
        endif
        !print *, "APE MAC: ", la, gacc_err_i, node(grv_bhIM), ndSize2,dr(MDIM+1), dist_lim_mpdp2, dist_mpdp2
 
        ! check if distance is larger than the limiting distance
        if (dist_mpdp2 > dist_lim_mpdp2) then
          node_accepted = .true.
        endif
      
      case(grv_bhMAC_MPE)
        ! Maximum Partial Error MAC can be used only with the minimum node 
        ! distance present in the tree in this version
        call Driver_abortFlash("Gravity_bhMAC: MPE misconfigured.")
 
      case(grv_bhMAC_SS)
        ! With SumSquare MAC, gr_bhMAC is never called.
        ! Error control is done using by calling gr_bhPartErr during 
        ! processing the priority queue
        call Driver_abortFlash("Gravity_bhMAC: SS misconfigured."//&
             " gr_bhMAC should not be called with SumSquare."//&
             " Set gr_bhPhysMACComm = .false.")
 
      case default
        ! shouldn't get here at all
        call Driver_abortFlash("Gravity_bhMAC: never should get here.")
    end select
  endif ! Use minimum node distance stored in the tree or not


  Gravity_bhMAC = node_accepted

  return
end function Gravity_bhMAC


