!!****if* source/Grid/localAPI/gr_mpoleCenterOfExpansion
!!
!! NAME
!!
!!  gr_mpoleCenterOfExpansion
!!
!! SYNOPSIS
!!
!!  gr_mpoleCenterOfExpansion (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of expansion for the multipoles.
!!  This routine is just the wrapper to call the appropriate routine according
!!  to the geometry present.
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!!***

subroutine gr_mpoleCenterOfExpansion (idensvar)

  implicit none
  
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleCenterOfExpansion
