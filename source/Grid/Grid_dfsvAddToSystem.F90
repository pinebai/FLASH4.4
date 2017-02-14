!!****f* source/Grid/Grid_dfsvAddToSystem
!!
!! NAME
!!  Grid_dfsvAddToSystem
!!
!! SYNOPSIS
!!
!!  call Grid_dfsvAddToSystem(
!!                            integer(IN), dimension(VARDESC_SIZE) :: baseVarDesc,
!!                              integer(IN)  :: iVar,
!!                              integer(IN)         :: iSrc,
!!                                integer(IN)  :: iFactorA,
!!                                integer(IN)  :: bcTypes(6),
!!                                real(IN)     :: bcValues(2,6),
!!                                real(IN)     :: dt,
!!                                real(IN)     :: chi,
!!                                real(IN)     :: scaleFact,
!!                                real(IN)     :: theta,
!!                                logical(IN)  :: solnIsDelta,
!!                        OPTIONAL,integer(IN) :: iFactorC,
!!                        OPTIONAL,integer(IN) :: iFactorD
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
!!         A,B,C,D are optional given scalar factors/terms that may depend
!!         on position; they are either physcially constant in time, or at
!!         least considered time-independent for the purpose of the operation
!!         implemented here (typically by computing their values from the
!!         solution state reached by the previous time step).
!!
!!      Presently it is used to do heat conduction and multigroup diffusion.
!!
!! ARGUMENTS
!!   iVar           : Variable on which the diffusion operatorion is performed (e.g., TEMP_VAR)
!!   iFactorA       :| Are factors in the equation with spatial variation.
!!   factorB        :| Factor C,D are optional and are generally used
!!   iFactorC       :| to represent emission/absorption in MGD.
!!   iFactorD       :| iFactorA is needed only for conduction.
!!                   | For this FcB variant of the interface, factorB is passed in
!!                   | allocated scratch arrays, i.e., the implementation will have
!!                   | to call Grid_ascGetBlkPtr to get at it.
!!   bcTypes        : Presently OUTFLOW, VACUUM is supported, DIRICHLET is untested.
!!   bcValues       : Values of iVar,factorB on boundary (DIRICHLET).                        
!!   dt             : The time step.
!!   scaleFact      : Factor by which the end solution is scaled (not used).
!!   chi            : useful for constant diffusion problems (not used).
!!   theta          : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   pass           : Ignored in unsplit solver.
!!                    pass=1 order of directional sweep X-Y-Z, 
!!                    pass=2 order of directional sweep Z-Y-X.
!!   iSrc           : Ignored.
!!   solnIsDelta    : Is the solution only a delta that the caller has to apply to the
!!                    temperature, rather than temperature itself (ignored).
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

subroutine Grid_dfsvAddToSystem(baseVarDesc,unkVarsDesc,diffCoeffDesc,absorpCoeffDesc,pass)
  implicit none

#include "constants.h"
  
  integer,intent(in),dimension(VARDESC_SIZE) :: unkVarsDesc, baseVarDesc
  integer,intent(in),dimension(VARDESC_SIZE) :: diffCoeffDesc, absorpCoeffDesc
  integer, intent(IN), OPTIONAL :: pass
  
end subroutine Grid_dfsvAddToSystem
