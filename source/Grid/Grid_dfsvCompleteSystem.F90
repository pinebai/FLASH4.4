!!****f* source/Grid/Grid_dfsvCompleteSystem
!!
!! NAME
!!  Grid_dfsvCompleteSystem
!!
!! SYNOPSIS
!!
!!  call Grid_dfsvCompleteSystem(
!!                            integer(IN), dimension(VARDESC_SIZE) :: baseVarDesc,
!!                        OPTIONAL,integer(IN) :: iFactorC,
!!                        OPTIONAL,integer(IN) :: iFactorD,
!!                        OPTIONAL,integer(IN) :: bcTypes(6),
!!                        OPTIONAL,real(IN)    :: bcValues(2,6),
!!                        OPTIONAL,real(IN)    :: dtNow,
!!                        OPTIONAL,real(IN)    :: thetaNow,
!!                        OPTIONAL,integer(IN) :: pass)
!!  DESCRIPTION 
!!
!!      Completes setting matrix and vector data for a diffusion problem.
!!
!!      The solver advances a generalized diffusion operator of the form
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
!!   baseVarDesc    : A variable on which the diffusion operation is performed, usually
!!                    an extra variable, used to identify the system
!!   iFactorA       :| Are factors in the equation with spatial variation.
!!   factorB        :| Factor C,D are optional and are generally used
!!   iFactorC       :| to represent emission/absorption in MGD.
!!   iFactorD       :|
!!                   | For this FcB variant of the interface, factorB is passed in
!!                   | allocated scratch arrays, i.e., the implementation will have
!!                   | to call Grid_ascGetBlkPtr to get at it.
!!   bcTypes        : Presently OUTFLOW, VACUUM is supported, DIRICHLET is untested.
!!   bcValues       : Values of iVar,factorB on boundary (DIRICHLET).                        
!!   dtNow          : The time step.
!!   thetaNow       : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   pass           : Ignored in unsplit solver.
!!                    pass=1 order of directional sweep X-Y-Z, 
!!                    pass=2 order of directional sweep Z-Y-X.
!!
!!
!! EXAMPLE  
!!
!!  ueleVarDesc = (/MGDC_VAR,1,CENTER,VD_DUR_PERM/)
!!  call Grid_dfsvCreateSystem(baseVarDesc=ueleVarDesc,ntotVars=numVarsForHypre, &
!!       dt=dt, &
!!       theta=mgdTheta,thetaC=thetaC,thetaD=thetaD, maxChunkSize=groupChunkSize)
!!
!!  call Grid_dfsvBeginSystem(baseVarDesc=ueleVarDesc)
!!
!!     diffCoeffDesc = (/COND_VAR,numGroups,CENTER,VD_DUR_PERM/)
!!     call Grid_dfsvAddToSystem(unkVarsDesc=groupVarsDesc,baseVarDesc=ueleVarDesc, &
!!                               firstHypreVar=ch, &
!!                               diffCoeffDesc=diffCoeffDesc, &
!!                               absorpCoeffDesc=absorpCoeffDesc, &
!!                               emissCoeffDesc=emissCoeffDesc, &
!!                               emissTermDesc=emissTermDesc, &
!!                               bcTypes=bcTypes, bcValues=bcValues)
!!
!!  call Grid_dfsvCompleteSystem(baseVarDesc=ueleVarDesc)
!!
!!  call Grid_advanceMultiDiffusion(groupVarsDesc(VARDESC_VAR), &
!!       iSrc=-1,  &
!!       iFactorA=-1, &
!!       bcTypes=bcTypes, &
!!       dtNow=dt, chi=1.0, scaleFact=scaleFact, theta=mgdtheta, solnIsDelta=.FALSE., &
!!       iFactorC=absorpCoeffDesc(VARDESC_VAR), &
!!       iFactorD=-1, &
!!       pass=pass)
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!  This implementation
!!    * supports: 
!!           1D, 2D, 3D Cartesian PARAMESH (with local refinement)
!!           1D, 2D, 3D Cartesian in UG (3D untested).
!!           1D Spherical in PARAMESH/UG
!!           1D, 2D Cylindrical in PARAMESH/UG (R-Z)
!!
!!    *  uses HYPRE library to solve AX = B
!!
!! SEE ALSO
!!
!!  Grid_dfsvCreateSystem
!!  Grid_dfsvBeginSystem
!!  Grid_dfsvAddToSystem
!!  Grid_advanceMultiDiffusion
!!  Diffuse_solveCoupledScalar
!!
!!***

!!REORDER(4): solnVec

subroutine Grid_dfsvCompleteSystem (baseVarDesc, factorADesc, unkVarsDesc,diffCoeffDesc,absorpCoeffDesc, &
     bcTypes, bcValues, &
     dtNow, thetaNow, pass)
  

  implicit none

#include "constants.h"
  
  integer, dimension(VARDESC_SIZE), intent(IN):: baseVarDesc
  integer, dimension(VARDESC_SIZE), intent(IN), OPTIONAL :: factorADesc
  integer,intent(in),dimension(VARDESC_SIZE),OPTIONAL :: unkVarsDesc
  integer,intent(in),dimension(VARDESC_SIZE),OPTIONAL :: diffCoeffDesc, absorpCoeffDesc
  integer, dimension(6),  intent(IN),OPTIONAL :: bcTypes
  real   , dimension(2,6),intent(IN),OPTIONAL :: bcValues
  real,    intent(IN), OPTIONAL :: dtNow
  real,    intent(IN), OPTIONAL :: thetaNow
  integer, intent(IN), OPTIONAL :: pass
  
  
end subroutine Grid_dfsvCompleteSystem
