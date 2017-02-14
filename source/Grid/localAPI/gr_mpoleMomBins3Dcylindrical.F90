!!****if* source/Grid/localAPI/gr_mpoleMomBins3Dcylindrical
!!
!! NAME
!!
!!  gr_mpoleMomBins3Dcylindrical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMomBins3Dcylindrical (integer (in) :: maxQtype)
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution in 3D cylindrical
!!  geometry. On output, the gr_mpoleMomentR and gr_mpoleMomentI arrays contain
!!  the mass moments over the regular and irregular solid harmonics. The moments
!!  are evaluated by an outer (threaded) loop over all radial bin types, ensuring
!!  separate evaluation of different radial bin types for all threads.
!!
!! ARGUMENTS
!!
!!  maxQtype : the total number of different radial bin types
!!
!!***

subroutine gr_mpoleMomBins3Dcylindrical (maxQtype)

  implicit none
  
  integer, intent (in) :: maxQtype
  
  return
end subroutine gr_mpoleMomBins3Dcylindrical
