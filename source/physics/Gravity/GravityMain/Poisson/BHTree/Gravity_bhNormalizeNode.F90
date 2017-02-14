!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhNormalizeNode
!!
!! NAME
!!
!!  Gravity_bhNormalizeNode
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhNormalizeNode(
!!                          real(in)       :: smr(MDIM),
!!                          real(inout)    :: node(:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during tree build after all subnodes are added to this node.
!!  Finish calculation of the B2 coefficient, if present in the node.
!!
!! ARGUMENTS
!!
!!  smr   : sum of mass times position of subnodes
!!  node  : array of the node
!!
!!
!!***

subroutine Gravity_bhNormalizeNode(smr, node)
  use Gravity_data, ONLY : useGravity, grv_bhIB2, grv_bhIB3, grv_bhNODE5, &
    grv_bhN5_NONE, grv_bhIM, grv_bhIX, grv_bhIY, grv_bhIZ
  implicit none
#include "constants.h"
#include "Flash.h"
  real, dimension(MDIM), intent(IN) :: smr
  real, dimension(:), intent(INOUT) :: node
  real :: rn, r2n, r3n

  if (.not. useGravity) return

  ! finish calculation of moments B_2 and B_3
  if (grv_bhNODE5 /= grv_bhN5_NONE) then
    !print *, "NN1: ", node
    r2n = node(grv_bhIX)**2 + node(grv_bhIY)**2 + node(grv_bhIZ)**2
    !rn  = sqrt(r2n)
    !r3n = r2n*rn

    ! First: subtract 3*r_n*(sum B2_sn + sum m_sn*r_sn**2) from B_3
    ! the term in parentheses is stored in B_2
    !node(grv_bhIB3) = node(grv_bhIB3) - 3.0*rn*node(grv_bhIB2)

    ! Second: subtract m_n*r_n**2 from B_2
    node(grv_bhIB2) = node(grv_bhIB2) - node(grv_bhIM) * r2n 
!    & - 2.0 * (node(grv_bhIX)*smr(IAXIS) + node(grv_bhIY)*smr(JAXIS) &
!    & + node(grv_bhIZ)*smr(KAXIS))

    ! Third: add 2 * m_n*r_n**3 to B_3 
    !node(grv_bhIB3) = node(grv_bhIB3) + 2.0*r3n*node(grv_bhIM)

    ! Fourth take absolute value
    !node(grv_bhIB3) = abs(node(grv_bhIB3))

    ! approximate value
    !node(grv_bhIB3) = sqrt(node(grv_bhIB2)**3/node(grv_bhIM))

    !print *, "NN2: ", node
  endif

  return
end subroutine Gravity_bhNormalizeNode

