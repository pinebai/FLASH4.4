!!****if* source/physics/Diffuse/DiffuseMain/Split/Diffuse_solveScalar
!!
!!  NAME 
!!
!!  Diffuse_solveScalar
!!
!!  SYNOPSIS
!!
!!  call Diffuse_solveScalar (integer, intent(IN) :: iVar,
!!                            integer, intent(IN) :: iFactorB,
!!                            integer, intent(IN) :: iFactorA,
!!                            integer, intent(IN) :: bcTypes(6),
!!                            real,    intent(IN) :: bcValues(2,6),
!!                            real,    intent(IN) :: dt,
!!                            real,    intent(IN) :: scaleFact,
!!                            real,    intent(IN) :: chi,
!!                            real,    intent(IN) :: theta,
!!                            integer, OPTIONAL, intent(IN) :: pass,
!!                            integer, intent(IN) :: blockCount,
!!                            integer,dimension(blockCount),intent(IN) :: blockList,
!!                            integer, intent(IN), OPTIONAL :: iFactorC,
!!                            integer, intent(IN), OPTIONAL :: iFactorD)
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
!!   iVar           :  Variable on which the diffusion operatorion is performed (e.g., TEMP_VAR)
!!   iFactorA       :| Are factors in the equation with spatial variation.
!!   iFactorB       :| Factor C,D are optional and are generally used
!!   iFactorC       :| to represent emission/absorption in MGD.
!!   iFactorD       :| iFactorA is needed only for conduction.
!!   bcTypes        :  Presently OUTFLOW, VACUUM is supported, DIRICHLET is untested.
!!   bcValues       :  Values of iVar,iFactorB on boundary (DIRICHLET).                        
!!   dt             :  The time step.
!!   scaleFact      :  Factor by which the end solution is scaled (not used).
!!   chi            :  useful for constant diffusion problems (not used).
!!   theta          :  varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   pass           :  Ignored in unsplit solver.
!!                     pass=1 order of directional sweep X-Y-Z, 
!!                     pass=2 order of directional sweep Z-Y-X.
!!   blockCount     :  The number of blocks in the list.   
!!   blockList      :  The list of blocks on which the solution must be updated.                    
!!
!!
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!  
!!
!!***


subroutine Diffuse_solveScalar(iVar, iFactorB, iFactorA, bcTypes, &
     bcValues, dt, scaleFact, chi, theta, pass, &
     blockCount, blockList, iFactorC, iFactorD)

 use diff_saData,    ONLY : diff_scaleFactThermSaTime
 use Grid_interface, ONLY : Grid_advanceDiffusion

 implicit none

 integer,                      intent(IN) :: iVar
 integer,                      intent(IN) :: iFactorB
 integer,                      intent(IN) :: iFactorA
 integer, OPTIONAL,            intent(IN) :: iFactorC
 integer, OPTIONAL,            intent(IN) :: iFactorD
 integer,                      intent(IN) :: bcTypes(6)
 real,                         intent(IN) :: bcValues(2,6)
 real,                         intent(IN) :: dt
 real,                         intent(IN) :: scaleFact
 real,                         intent(IN) :: chi
 real,                         intent(IN) :: theta
 integer, OPTIONAL,            intent(IN) :: pass
 integer,                      intent(IN) :: blockCount
 integer,dimension(blockCount),intent(IN) :: blockList

 !========================================================================= 
 logical :: solnIsDelta = .FALSE.
 integer :: iSrc = -1 
 
 
 call Grid_advanceDiffusion (iVar, iVar, iFactorB, iFactorA, bcTypes, bcValues,  &
      dt, chi, scaleFact,theta, solnIsDelta, iFactorC, iFactorD, pass)
   
  return

end subroutine Diffuse_solveScalar
