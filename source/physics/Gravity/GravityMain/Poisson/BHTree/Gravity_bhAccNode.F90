!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhAccNode
!!
!! NAME
!!
!!  Gravity_bhAccNode
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhAccNode(
!!                          real(in)       :: subnode(:),
!!                          real(inout)    :: accnode(:)
!!        )
!!
!! DESCRIPTION
!!
!!  Called during tree build. Adds values of subnode into accnode.
!!
!! ARGUMENTS
!!
!!  subnode     : array of the node of the tree which is added to accnode
!!  accnode     : array of the node into which subnode contribution is added
!!
!!
!!
!!***

subroutine Gravity_bhAccNode(subnode, accnode)
  use Gravity_data, ONLY : useGravity, grv_bhIB2, grv_bhIB3, grv_bhNODE5, &
    grv_bhN5_NONE, grv_bhIM, grv_bhIX, grv_bhIY, grv_bhIZ
  implicit none
  real, dimension(:), intent(IN)  :: subnode
  real, dimension(:), intent(INOUT) :: accnode
  real :: r2sn !, rsn, r3sn

  if (.not. useGravity) return

  if (grv_bhNODE5 /= grv_bhN5_NONE) then
    r2sn = (subnode(grv_bhIX)**2 + subnode(grv_bhIY)**2 &
    &    + subnode(grv_bhIZ)**2)
    !rsn  = sqrt(r2sn)
    !r3sn = rsn*r2sn

    ! calculate second order moment B_2 to be used for formula 9 in S&W94
    ! B_2 = sum B_2_sn + sum m_sn*(r_sn - r_an)**2 = 
    ! = sum B_2_sn + sum m_sn*r_sn**2 - m_an*r_an**2
    ! (since sum m_sn*r_sn = m_an*r_an)
    ! in this moment we accumulate B_2_sn moments of subnodes plus
    ! a sum m_sn * r_sn**2
    ! the last term is subtracted in Gravity_bhNormalizeNode
    accnode(grv_bhIB2) = accnode(grv_bhIB2) + subnode(grv_bhIB2) &
    & + subnode(grv_bhIM) * r2sn

    ! calculate second third moment B_3 to be used for formula 9 in S&W94
    ! B_3 = sum B_3_sn + 3*sum B_2_sn*(r_sn-r_an) + sum m_sn*(r_sn - r_an)**3 = 
    ! = sum B_3_sn + 3* sum B_2_sn*r_sn + sum m_sn*r_sn**3 
    ! - 3*r_an*(sum B_2_sn + sum m_sn*r_sn**2) + 2*m_an*r_an**3
    ! in this moment we accumulate first three terms
    ! remaining terms are added in Gravity_bhNormalizeNode

    !accnode(grv_bhIB3) = accnode(grv_bhIB3) + subnode(grv_bhIB3) &
    !& + 3 * subnode(grv_bhIB2) * rsn + subnode(GR_TREE_IM) * r3sn

    !print *, "AN: ", accnode
  endif

  return
end subroutine Gravity_bhAccNode
