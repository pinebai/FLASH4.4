!!****f* source/physics/Diffuse/Diffuse_solveScalar
!!
!!  NAME 
!!
!!  Diffuse_solveScalar
!!
!!  SYNOPSIS
!!
!!  call Diffuse_solveScalar (integer(IN) :: iVar,
!!                            integer(IN) :: iFactorB,
!!                            integer(IN) :: iFactorA,
!!                            integer(IN) :: bcTypes(6),
!!                            real(IN)    :: bcValues(2,6),
!!                            real(IN)    :: dt,
!!                            real(IN)    :: scaleFact,
!!                            real(IN)    :: chi,
!!                            real(IN)    :: theta,
!!                            integer(IN), OPTIONAL :: pass,
!!                            integer(IN) :: blockCount,
!!                            integer(IN) :: blockList(blockCount),
!!                            integer(IN), OPTIONAL :: iFactorC,
!!                            integer(IN), OPTIONAL :: iFactorD)
!!
!!
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
!! 
!!   iVar           : Variable on which the diffusion operation is performed (e.g TEMP_VAR)
!!   iFactorA       :| Are factors in the equation with spatial variation.
!!   iFactorB       :| Factor C,D are optional and are generally used
!!   iFactorC       :| to represent emission/absorption in MGD.
!!   iFactorD       :| iFactorA is not needed (and thus set to constant 1.0) for
!!                   | radiation diffusion.
!!                   | If iFactorB is negative, the factor B (whose role is that of
!!                   | diffusion coefficient) is passed to the GridSolver implementation
!!                   | in face-centered allocatable scratch buffers (gds type FACEX,
!!                   ! FACEY, FACEZ) instead of UNK. [For HYPRE GridSolver only]
!!   bcTypes        : Presently OUTFLOW, VACUUM, DIRICHLET are supported, with additional
!!                    limited support for OUTSTREAM.
!!   bcValues       : Values of iVar,iFactorB on boundary (DIRICHLET).                        
!!   dt             : The time step.
!!   scaleFact      : Factor by which the end solution is scaled (not used).
!!   chi            : useful for constant diffusion problems (not used).
!!   theta          : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   pass           : Ignored in unsplit solver.
!!                    pass=1 order of directional sweep X-Y-Z, 
!!                    pass=2 order of directional sweep Z-Y-X.
!!   blockCount     : The number of blocks in the list.   
!!   blockList      : The list of blocks on which the solution must be updated.
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!  Stub implementation.              
!!
!! SEE ALSO
!!
!!  Grid_ascGetBlkPtr
!!
!!***

subroutine Diffuse_solveScalar (iVar, iFactorB, iFactorA, bcTypes, bcValues, &
     dt, scaleFact, chi, theta, pass, blockCount, blockList, iFactorC, iFactorD)
  
  implicit none
  
  integer, intent(IN):: iVar
  integer, intent(IN):: iFactorB
  integer, intent(IN):: iFactorA
  integer, intent(IN):: bcTypes(6)
  integer, intent(IN):: blockCount
  integer, dimension(blockCount),intent(IN):: blockList
  real, intent(IN):: bcValues(2,6)
  real, intent(IN):: dt
  real, intent(IN):: scaleFact
  real, intent(IN):: chi
  real, intent(IN):: theta
  integer, OPTIONAL,intent(IN):: pass
  integer, OPTIONAL,intent(IN):: iFactorC
  integer, OPTIONAL,intent(IN):: iFactorD
  
  return
  
end subroutine Diffuse_solveScalar
