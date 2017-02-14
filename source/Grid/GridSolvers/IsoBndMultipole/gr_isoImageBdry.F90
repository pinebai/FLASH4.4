!!****if* source/Grid/GridSolvers/IsoBndMultipole/gr_isoImageBdry
!!
!! NAME
!! 
!!  gr_isoImageBdry
!!
!! SYNOPSIS
!!
!!  gr_isoImageBdry(integer(IN) :: iiden,
!!                           integer(IN) :: iipot,
!!                           real(IN)    :: poisfact)
!!
!! DESCRIPTION
!!   Compute the boundary values of the image mass potential, this is
!!   necessary to apply isolated boundary conditions when using
!!   the multigrid solver
!!  
!! ARGUMENTS
!!    iiden : The variable index for density
!!    iipot : variable index for potential
!!    poisfact : factor to be used in calculation
!!
!!***


subroutine gr_isoImageBdry (iiden, iipot, poisfact)


  implicit none

  integer,intent(IN)       :: iiden, iipot
  real,intent(IN)          :: poisfact

  !==========================================================================


  call gr_isoFindMassCenter (iiden)
  call gr_isoMpoleMoments (iiden)
  call gr_isoMpolePotential (iipot, poisfact)

  !========================================================================

  return
end subroutine gr_isoImageBdry
