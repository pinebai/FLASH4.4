!!****if* source/physics/RadTrans/RadTransMain/MGD/RadTrans_mgdSetBc
!!
!! NAME
!!
!!  RadTrans_mgdSetBc
!!
!! SYNOPSIS
!!
!!  call RadTrans_mgdSetBc(integer(in) :: ig,
!!                         integer, optional(in) :: bctypes,
!!                         real, optional(in) :: bcvalues,
!!                         integer, optional(in) :: f,
!!                         integer, optional(in) :: bctype,
!!                         real, optional(in) :: bcvalue)
!!
!! DESCRIPTION
!!
!! Set the boundary condition for a specific group. Can be invoked in
!! one of two modes:
!! 
!! 1. bctypes/bcvalues present: This allows you to easily set the BC on all
!!    of the faces of the domain
!! 2. f/bctype/bcvalue present: This allows you to easily set the BC on a
!!    specific face of the domain
!!
!!
!! ARGUMENTS
!!
!!   ig : energy group number
!!   bctypes : boundary condition type on each boundary
!!   bcvalues : boundary condition value (for dirichlet) on each boundary
!!   f : boundary number
!!   bctype : type of bondary for face f
!!   bcvalue : dirichlet value for face f
!!
!!***

subroutine RadTrans_mgdSetBc(ig, bcTypes, bcValues, f, bcType, bcValue)
  use rt_data, only: rt_mgdDomainBC, rt_mgdBcVals
  implicit none
  integer, intent(in) :: ig

  integer, optional, intent(in) :: bcTypes(6)
  real, optional, intent(in) :: bcValues(6)

  integer, optional, intent(in) :: f
  integer, optional, intent(in) :: bcType
  real, optional, intent(in) :: bcValue

  if(present(bcType) .and. present(f)) then
     rt_mgdDomainBC(ig,f) = bcType
  end if

  if(present(bcValue) .and. present(f)) then
     rt_mgdBcVals(ig,f) = bcValue
  end if

  if(present(bcTypes)) then
     rt_mgdDomainBC(ig,:) = bcTypes
  end if

  if(present(bcValues)) then
     rt_mgdBcVals(ig,:) = bcValues
  end if

end subroutine RadTrans_mgdSetBc
