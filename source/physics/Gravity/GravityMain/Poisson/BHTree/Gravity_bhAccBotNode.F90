!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhAccBotNode
!!
!! NAME
!!
!!  Gravity_bhAccBotNode
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhAccBotNode(
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real(in)       :: locCoords(:,:,:),
!!                          real,pointer   :: solnData(:,:,:,:),
!!                          real(in)       :: botnode(:),
!!                          real(inout)    :: accnode(:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during tree build. Adds values of botnode into accnode.
!!
!! ARGUMENTS
!!
!!  blockno     : number of block into which the node belongs
!!  point       : indeces of the cell (botnode) in the block
!!  blkLimits   : limits of indeces in the block
!!  locCoords   : array with coordinates of allblocks on a given cpu
!!  solnData    : solution data from the grid
!!  botnode     : array of the bottom-most (leaf) node of the tree, i.e. a grid cell
!!  accnode     : array of the node into which botnode contribution is added
!!
!!
!!***

subroutine Gravity_bhAccBotNode(blockno, point, blkLimits, locCoords, solnData, botnode, accnode)
  use Gravity_data, ONLY : useGravity, grv_bhIB2, grv_bhNODE5, grv_bhN5_NONE, grv_bhIM
  use Grid_interface, ONLY : Grid_getDeltas
  implicit none
#include "constants.h"
#include "Flash.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, dimension(:,:,:), intent(IN) :: locCoords
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, dimension(:), intent(IN) :: botnode
  real, dimension(:), intent(INOUT) :: accnode
  real, dimension(MDIM) :: del
  real :: B2_bn, r2bn !, B3_bn, rbn, r3bn

  if (.not. useGravity) return

  if (grv_bhNode5 /= grv_bhN5_NONE) then

    call Grid_getDeltas(blockno, del)
    r2bn = (locCoords(point(IAXIS)-blkLimits(LOW,IAXIS)+1, IAXIS, blockno)**2 &
    &    +  locCoords(point(JAXIS)-blkLimits(LOW,JAXIS)+1, JAXIS, blockno)**2 &
    &    +  locCoords(point(KAXIS)-blkLimits(LOW,KAXIS)+1, KAXIS, blockno)**2)
    !rbn  = sqrt(r2bn)
    !r3bn = r2bn*rbn

    ! calculate second order moment B_2 to be used for formula 9 in S&W94
    ! B_2 = sum B_2_sn + sum m_sn*(r_sn - r_an)**2 = 
    ! = sum B_2_sn + sum m_sn*r_sn**2 - m_an*r_an**2
    ! (since sum m_sn*r_sn = m_an*r_an)
    ! in this moment we accumulate B_2_sn moments of subnodes plus
    ! a sum m_sn * r_sn**2
    ! the last term is subtracted in Gravity_bhNormalizeNode
    B2_bn = (1.0/12.0) * botnode(grv_bhIM) &
    & * (del(IAXIS)**2 + del(JAXIS)**2 + del(KAXIS)**2)
    accnode(grv_bhIB2) = accnode(grv_bhIB2) + B2_bn + botnode(grv_bhIM) * r2bn

    ! calculate second third moment B_3 to be used for formula 9 in S&W94
    ! B_3 = sum B_3_sn + 3*sum B_2_sn*(r_sn-r_an) + sum m_sn*(r_sn - r_an)**3 = 
    ! = sum B_3_sn + 3* sum B_2_sn*r_sn + sum m_sn*r_sn**3 
    ! - 3*r_an*(sum B_2_sn + sum m_sn*r_sn**2) + 2*m_an*r_an**3
    ! in this moment we accumulate first three terms
    ! remaining terms are added in Gravity_bhNormalizeNode

    ! B3 of a botnode assumes that grid cells are cubic!!!
    !B3_bn = grv_bhB3_CONST * botnode(grv_bhIM) * del(IAXIS)**3
    !accnode(grv_bhIB3) = accnode(grv_bhIB3) + B3_bn + 3 * B2_bn * rbn &
    !& + botnode(grv_bhIM) * r3bn

  endif

  return
end subroutine Gravity_bhAccBotNode
