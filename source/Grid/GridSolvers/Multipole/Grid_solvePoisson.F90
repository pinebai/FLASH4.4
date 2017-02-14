!!****if* source/Grid/GridSolvers/Multipole/Grid_solvePoisson
!!
!! NAME
!!
!!  Grid_solvePoisson
!!
!!
!! SYNOPSIS
!!
!!   Grid_solvePoisson(integer(IN) :: iSoln,
!!                     integer(IN) :: iSrc, 
!!                  integer(6)(IN) :: bcTypes,
!!                   real(2,6)(IN) :: bcValues,
!!                     real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This module implements the multipole
!!   summation method for isolated problems.  Periodic problems are
!!   not supported; the multipole method should be used only for
!!   approximately spherical matter distributions.
!!
!!
!! ARGUMENTS
!!
!!  iSoln -  index to variable containing potential
!!  iSrc - index to variable containing density
!!  bcTypes - boundary types along various faces,
!!             only used in verifying that they are isolated
!!  bcValues - the values to boundary conditions, currently not used
!!  poisfact -  factor to be used in calculation
!!
!!***



subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)


  use Grid_interface, ONLY: GRID_PDE_BND_ISOLATED
  use Driver_interface, ONLY: Driver_abortFlash

  implicit none

#include "Flash.h"

  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes(6)
  real, intent(in)       :: bcValues(2,6)
  real, intent(inout)    :: poisfact
  
  !--------------------------------------------------------------------------
  
  if (ANY(bcTypes(1:2*NDIM) /= GRID_PDE_BND_ISOLATED)) call Driver_abortFlash &
       ("FATAL: Multipole Poisson solver requires isolated boundaries")

  call gr_mpoleCenterOfMass(iSrc)
  call gr_mpoleMoments(iSrc)
  call gr_mpolePotential(iSrc, iSoln, poisfact)

  
  return
end subroutine Grid_solvePoisson
