!!****if* source/Grid/localAPI/gr_mpoleMoments
!!
!! NAME
!!
!!  gr_mpoleMoments
!!
!! SYNOPSIS
!!
!!  gr_mpoleMoments (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution,
!!  assuming the center of mass and the total mass have first been
!!  computed.  On output, the Moment_R () and Moment_I () arrays
!!  contain the mass moments over the regular and irregular
!!  solid harmonics. This routine calls the appropriate subroutines
!!  according to the geometry specified.
!!
!! ARGUMENTS
!!
!!  idensvar -- the index of the density variable
!!
!!***

subroutine gr_mpoleMoments (idensvar)

  implicit none
    
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleMoments
