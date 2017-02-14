!!****f* source/Grid/Grid_dfsvCreateSystem
!!
!! NAME
!!  Grid_dfsvCreateSystem
!!
!! SYNOPSIS
!!
!!  call Grid_dfsvCreateSystem(
!!                            integer(IN), dimension(VARDESC_SIZE) :: baseVarDesc,
!!                            integer(IN)  :: ntotVars,
!!                        OPTIONAL,integer(IN) :: pass)
!!  DESCRIPTION 
!!
!!      This routine advances a generalized diffusion operator of the form
!!
!!         A*(df/dt) + C*f = div(B*grad(f)) + D ,
!!
!!      where
!!         f = f(x,t) is the  Variable to be diffused (x=1D..3D position);
!!
!!      Presently it is used to do heat conduction and multigroup diffusion.
!!
!! ARGUMENTS
!!   baseVarDesc    : A variable on which the diffusion operation is performed, usually
!!                    an extra variable, used to identify the system
!!   ntotVars       : total number or variables to solve for / equations to solve
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!  This implementation
!!    * supports: 
!!           1D, 2D PARAMESH (with local refinement)
!!           1D, 2D, 3D in UG (3D untested).
!!           1D Spherical in PARAMESH/UG
!!           1D, 2D Cylindrical in PARAMESH/UG (R-Z)
!!
!!    *  uses HYPRE library to solve AX = B
!!
!!
!!***

subroutine Grid_dfsvCreateSystem (baseVarDesc, ntotVars)
  

  implicit none

#include "constants.h"
  
  integer, dimension(VARDESC_SIZE), intent(IN):: baseVarDesc
  integer, intent(IN) :: ntotVars
  
end subroutine Grid_dfsvCreateSystem
