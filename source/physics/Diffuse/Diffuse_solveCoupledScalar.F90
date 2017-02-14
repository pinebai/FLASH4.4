!!****f* source/physics/Diffuse/Diffuse_solveCoupledScalar
!!
!!  NAME 
!!
!!  Diffuse_solveCoupledScalar
!!
!!  SYNOPSIS
!!
!!  call Diffuse_solveCoupledScalar (
!!                            integer(IN), dimension(VARDESC_SIZE) :: unkVarsDesc,
!!                            integer(IN), dimension(VARDESC_SIZE), OPTIONAL :: xtraVarDesc,
!!                            integer(IN), dimension(VARDESC_SIZE) :: diffCoeffDesc,
!!                            integer(IN), dimension(VARDESC_SIZE), OPTIONAL :: absorpCoeffDesc,
!!                            integer(IN)  :: bcTypes(6),
!!                            real(IN)     :: bcValues(2,6),
!!                            real(IN)     :: dt,
!!                            real(IN)     :: scaleFact,
!!                            real(IN)     :: theta,
!!                            integer(IN), OPTIONAL :: pass,
!!                            integer(IN), dimension(VARDESC_SIZE), OPTIONAL :: iFactorD,
!!                            integer(IN), dimension(VARDESC_SIZE) :: iFactorA)
!!
!!
!!  DESCRIPTION 
!!
!!    Advance a diffusion problem for multiple coupled scalar variables.
!!
!!      This routine advances a generalized diffusion operator of the form
!!
!!         A*(df/dt) + C*f = div(B*grad(f)) + D ,
!!
!!      where
!!         f = f(x,t) is the  Variable to be diffused (x=1D..3D position);
!!
!!         A,B,C,D are optional given scalar factors/terms that may depend
!!         on position; they are either physcially constant in time, or at
!!         least considered time-independent for the purpose of the operation
!!         implemented here (typically by computing their values from the
!!         solution state reached by the previous time step).
!!
!!      Presently it is used to do multigroup diffusion coupled to electron energy
!!      via absorption and emission.
!!
!! ARGUMENTS
!! 
!!   unkVarsDesc    : describes the main unknown variables on which the diffusion operation is performed.
!!                    NumGroups variables.
!!   xtraVarDesc    : describes an additional unknown variable coupled to the main unknown variables.
!!                    One variable.
!!   iFactorA       : unused
!!   diffCoeffDesc  : Sometimes referred to as "Factor B".
!!                    Represents the diffusion coefficient in MGD.
!!                    NumGroups variables.
!!   absorpCoeffDesc: Represents absorption in MGD.
!!                    NumGroups variables.
!!   iFactorD       : Represents emission in MGD.
!!                    Currently ignored, absorpCoeffDesc is used for emission as well as absorption.
!!                    NumGroups variables.
!!   bcTypes        : Presently OUTFLOW, VACUUM is supported, DIRICHLET is untested.
!!   bcValues       : Values of unkVars,diffCoeffDesc on boundary (DIRICHLET).                        
!!   dt             : The time step.
!!   scaleFact      : Factor by which the end solution is scaled (not used).
!!   theta          : varies scheme (0-> Explicit, 1-> backward Euler, 0.5 -> Crank Nicholson
!!   pass           : Ignored in unsplit solver.
!!                    pass=1 order of directional sweep X-Y-Z, 
!!                    pass=2 order of directional sweep Z-Y-X.
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!  Stub implementation.              
!!
!!  VARDESC_SIZE is defined in constants.h.  
!!
!!***

#include "constants.h"

subroutine Diffuse_solveCoupledScalar (unkVarsDesc, xtraVarDesc, diffCoeffDesc, absorpCoeffDesc, bcTypes, bcValues, &
     dt, scaleFact, theta, pass, iFactorD, iFactorA)
  
  implicit none
  
  integer, dimension(VARDESC_SIZE), intent(IN):: unkVarsDesc
  integer, dimension(VARDESC_SIZE), intent(IN):: xtraVarDesc
  integer, dimension(VARDESC_SIZE), intent(IN):: diffCoeffDesc
  integer, dimension(VARDESC_SIZE), OPTIONAL,intent(IN):: absorpCoeffDesc
  integer, intent(IN):: bcTypes(6)
  real, intent(IN),OPTIONAL:: bcValues(2,6)
  real, intent(IN):: dt
  real, intent(IN):: scaleFact
  real, intent(IN):: theta
  integer, OPTIONAL,intent(IN):: pass
  integer, dimension(VARDESC_SIZE), OPTIONAL,intent(IN):: iFactorD
  integer, dimension(VARDESC_SIZE), OPTIONAL, intent(IN):: iFactorA
  
  return
  
end subroutine Diffuse_solveCoupledScalar
