!!****if* source/physics/Diffuse/DiffuseMain/Split/Diffuse_advance1D
!!
!! NAME
!!
!!  Diffuse_advance1D
!!
!! SYNOPSIS
!!
!!  call Diffuse_advance1D(integer(IN)            :: nx,
!!                         integer(IN)            :: ny,
!!                         integer(IN)            :: nz,
!!                         integer(IN)            :: idirection,
!!                         integer(IN)            :: bctype(2),
!!                         real(IN)               :: del,
!!                         real(IN)               :: dt,
!!                         real(IN)               :: theta,
!!                         real(IN)               :: chi,
!!                         real(IN)               :: cond(nx,ny,nz),
!!                         real(INOUT)            :: soln(nx,ny,nz),
!!                         logical(IN)            :: useCFArray,
!!                         real,OPTIONAL, pointer :: ConstFact(:,:,:))
!!
!! DESCRIPTION
!!
!!  Subroutine to advance Diffusion in one dimension.
!!
!!  This routine advances a diffusive operator (think heat conduction) in time.
!!
!!  Generic form of the equation a*dV/dt = d/dx(b*dV/dx) + d/dy(b*dV/dy) + d/dz(b*dV/dz)
!!                               i.e for Heat conduction, a=rho*Cv, b=conductivity and V=Temp
!!
!!
!! ARGUMENTS
!!
!!   nx         : Total grid points along x-direction.
!!
!!   ny         : Total grid points along y-direction.
!!
!!   nz         : Total grid points along z-direction.
!!
!!   idirection : Direction of solve 
!!
!!   bctype     : Left and Right boundary condition type.
!!                Supports PERIODIC, OUTFLOW (tested).
!!                         DIRICHLET (untested).
!!
!!   del        : dx,dy or dz based on direction of solve
!!
!!   dt         : timestep.
!!
!!   theta      : 0.5 CN Scheme, 1.0 Fully Implicit, 0.0 Explicit
!!
!!   chi        : Diffusion coefficient, not used anymore.
!!
!!   cond       : Spatial Conductivity variation (K), "b" in the general diffusion equation.
!!
!!   soln       : discretised  variation  of variable to be diffused (input). On Output contains the diffused solution.
!!
!!   useCFArray : Flag to determine if we have ConstFact(:,:,:)
!!
!!   ConstFact  : The index of coefficient of temporal derivative (think A)
!!                rho*cv in Heat conduction.
!!
!! SEE ALSO
!!
!!  Grid_advanceDiffusion
!!  Diffuse
!!
!! NOTES
!!
!!  This file includes a small subroutine that implements the Thomas Algorithm 
!!  for a tridiagonal matrix.
!!
!!  This public API subroutine is so far used only from
!!  the pencil grid implementation of Grid_advanceDiffusion.
!!  The pencil grid implementation of Grid_advanceDiffusion, in turn,
!!  is currently only used by the Split implementation of the DiffuseMain
!!  subunit.
!!***

#include "constants.h"
#include "Flash.h"

subroutine Diffuse_advance1D (nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray, &
                              ConstFact,useMatCArray, MatC,useMatDArray,MatD)

  implicit none

  integer,                intent(IN)    :: nx, ny, nz
  integer,                intent(IN)    :: iDirection
  integer, dimension (2), intent(IN)    :: bcType
  real   , dimension (2), intent(IN)    :: bcValues
  real   ,                intent(IN)    :: DEL
  real   ,                intent(IN)    :: dt
  real   ,                intent(IN)    :: theta
  real   ,                intent(IN)    :: chi
  real   ,                intent(IN)    :: cond(nx,ny,nz)
  real   ,                intent(INOUT) :: soln(nx,ny,nz)
  logical,                intent(IN)    :: useCFArray
  real   ,OPTIONAL,       pointer       :: ConstFact(:,:,:) !ConstFact(nx,ny,nz) 
  logical,                intent(IN)    :: useMatCArray    
  real   ,OPTIONAL,       pointer       :: MatC(:,:,:)   !ConstFact(nx,ny,nz) 
  logical,                intent(IN)    :: useMatDArray
  real   ,OPTIONAL,       pointer       :: MatD(:,:,:)   !ConstFact(nx,ny,nz) 
    

  if (.NOT. associated(ConstFact)) then
     if (.NOT. associated(matC)) then
        if (.NOT. associated(matD)) then
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray,cond, &
                                   useMatCArray, cond, useMatDArray, cond)
        else
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray,cond, &
                                   useMatCArray, cond, useMatDArray, matD)
        end if
     else
        if (.NOT. associated(matD)) then
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray,cond, &
                                   useMatCArray, matC, useMatDArray, cond)
        else
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray,cond, &
                                   useMatCArray, matC, useMatDArray, matD)
        end if
     end if
  else
     if (.NOT. associated(matC)) then
        if (.NOT. associated(matD)) then
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray, &
                                   constFact, useMatCArray, cond, useMatDArray, cond)
        else
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray, &
                                   ConstFact, useMatCArray, cond, useMatDArray, matD)
        end if
     else
        if (.NOT. associated(matD)) then
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray, &
                                   constFact, useMatCArray, matC, useMatDArray, cond)
        else
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray, &
                                   ConstFact, useMatCArray, MatC, useMatDArray, MatD)
        end if
        
     end if
  end if
end subroutine Diffuse_advance1D

subroutine Diffuse_advance1D_1DArr (nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray, &
                                    ConstFact,useMatCArray, MatC, useMatDArray, MatD)

  implicit none

  integer,                intent(IN)    :: nx, ny, nz
  integer,                intent(IN)    :: iDirection
  integer, dimension (2), intent(IN)    :: bcType
  real   , dimension (2), intent(IN)    :: bcValues
  real   ,                intent(IN)    :: DEL
  real   ,                intent(IN)    :: dt
  real   ,                intent(IN)    :: theta
  real   ,                intent(IN)    :: chi
  real   ,                intent(IN)    :: cond(nx*ny*nz)
  real   ,                intent(INOUT) :: soln(nx*ny*nz)
  logical,                intent(IN)    :: useCFArray
  real   ,                pointer       :: ConstFact(:) !ConstFact(nx*ny*nz) 
  logical,                intent(IN)    :: useMatCArray
  real   ,                pointer       :: MatC(:)   !matC(nx,ny,nz) 
  logical,                intent(IN)    :: useMatDArray
  real   ,                pointer       :: MatD(:)   !matD(nx,ny,nz) 

!!$  if (.NOT. associated(ConstFact)) then
!!$     useCFArray = .FALSE.
!!$  end if
!!$  if (.NOT. associated(matC)) then
!!$     useMatCArray = .FALSE.
!!$  end if
!!$  if (.NOT. associated(matD)) then
!!$     useMatDArray = .FALSE.
!!$  end if
  if (.NOT. associated(ConstFact)) then
     if (.NOT. associated(matC)) then
        if (.NOT. associated(matD)) then
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray,cond, &
                                   useMatCArray, cond, useMatDArray, cond)
        else
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray,cond, &
                                   useMatCArray, cond, useMatDArray, matD)
        end if
     else
        if (.NOT. associated(matD)) then
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray,cond, &
                                   useMatCArray, matC, useMatDArray, cond)
        else
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray,cond, &
                                   useMatCArray, matC, useMatDArray, matD)
        end if
     end if
  else
     if (.NOT. associated(matC)) then
        if (.NOT. associated(matD)) then
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray, &
                                   constFact, useMatCArray, cond, useMatDArray, cond)
        else
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray, &
                                   ConstFact, useMatCArray, cond, useMatDArray, matD)
        end if
     else
        if (.NOT. associated(matD)) then
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray, &
                                   constFact, useMatCArray, matC, useMatDArray, cond)
        else
           call diff_advance1DImpl(nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray, &
                                   ConstFact, useMatCArray, MatC, useMatDArray, MatD)
        end if
        
     end if
  end if

end subroutine Diffuse_advance1D_1DArr


! Do not make an explicit interface for the following subroutine visible
! to the wrappers above - falling back to Fortran77-like level of knowledge
! is the whole point here! - KW
subroutine diff_advance1DImpl (nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray,ConstFact, &
     useMatCArray, matC, useMatDArray, matD)

  use Diffuse_data, ONLY : diff_asol, diff_geometry, diff_speedlt

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC, &
                               GRID_PDE_BND_NEUMANN,  &
                               GRID_PDE_BND_DIRICHLET


  implicit none

  integer,                intent(IN)    :: nx, ny, nz
  integer,                intent(IN)    :: iDirection
  integer, dimension (2), intent(IN)    :: bcType
  real   , dimension (2), intent(IN)    :: bcValues
  real   ,                intent(IN)    :: DEL
  real   ,                intent(IN)    :: dt
  real   ,                intent(IN)    :: theta
  real   ,                intent(IN)    :: chi
  real   ,                intent(IN)    :: cond(nx,ny,nz)
  real   ,                intent(INOUT) :: soln(nx,ny,nz)

  logical,                intent(IN)    :: useCFArray
  real   ,                intent(IN)    :: ConstFact(nx,ny,nz) 

  logical,                intent(IN)    :: useMatCArray
  real   ,                intent(IN)    :: matC(nx,ny,nz)

  logical,                intent(IN)    :: useMatDArray
  real   ,                intent(IN)    :: matD(nx,ny,nz)
  


  !**** Local variables ****
  real                           :: T_0, T_N
  real                           :: DELSQ
  real, allocatable,dimension(:) :: AA,BB,DD,S, U
  integer                        :: NL, LL, ML
  integer                        :: IL, JL, KL
  real                           :: CondL, CondR, cDiv
  real                           :: bcL, FInc, bcLE, bcLI

  real                           :: r, dr, rphalf, rmhalf
  real                           :: coeff
 
  
  ! =====================================================================

  DELSQ = DEL**2
!!$  if (iDirection == SWEEP_X) then
!!$     NL = nx
!!$     LL = ny
!!$     ML = nz
!!$  elseif (iDirection == SWEEP_Y) then
!!$     NL = ny
!!$     LL = nx
!!$     ML = nz
!!$  elseif (iDirection == SWEEP_Z) then
!!$     NL = nz
!!$     LL = ny
!!$     ML = nx
!!$  else
!!$     call Driver_abortFlash("[Diffuse_advance1D] invalid sweep direction: must be SWEEP_X, SWEEP_Y or SWEEP_Z")
!!$  endif

  NL = nx
  LL = ny
  ML = nz
  
  allocate (AA(NL))  ! Coefficient ahead of the diagonal.
  allocate (DD(NL))  ! Coefficient of the diagonal.
  allocate (BB(NL))  ! Coefficient after the diagonal.
  allocate ( S(NL))  ! RHS: Source Term.
  allocate ( U(NL))
 ! ************************************* X SWEEP ***************************************
 cDiv = 1.0                     ! set beelow to 1/ConstFact(...) if ConstFact is provided
 do KL = 1, ML
    do JL = 1, LL   

       if  ((diff_geometry == CYLINDRICAL) .and. iDirection == SWEEP_X) then 
          do IL=2, NL-1 
             !! psuedo code to compute AA,BB,DD for Radial diffusion.       
             ! Since it's uniform grid
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = r*dr**2              
          
             if (useCFArray) cDiv = 1.0/ConstFact(IL,JL,KL)           
             CondL = 0.5 * (cond(IL-1,JL,KL) + cond(IL,JL,KL))*cDiv ! (IL,JL,KL)
             CondR = 0.5 * (cond(IL+1,JL,KL) + cond(IL,JL,KL))*cDiv ! (IL,JL,KL)           
             !! The matrix is modified.
             AA(IL) = -theta*CondL*rmhalf*dt/DELSQ
             BB(IL) = -theta*CondR*rphalf*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*rphalf + CondL*rmhalf)*dt/DELSQ
             
             S(IL)  = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ) *(CondL*(rmhalf)*soln(IL-1,JL,KL) &
                  -(CondL*(rmhalf)+CondR*(rphalf))*soln(IL,JL,KL) &
                  + CondR*(rphalf)*soln(IL+1,JL,KL))                                      
   
          end do

       else if  ((diff_geometry == SPHERICAL) .and. iDirection == SWEEP_X) then         

          do IL=2, NL-1 
             
             if (useCFArray) cDiv = 1.0/ConstFact(IL,JL,KL)           
             CondL = 0.5 * (cond(IL-1,JL,KL) + cond(IL,JL,KL))*cDiv ! (IL,JL,KL)
             CondR = 0.5 * (cond(IL+1,JL,KL) + cond(IL,JL,KL))*cDiv ! (IL,JL,KL)     

             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = dr*((rphalf**3) -(rmhalf**3))/3.0
             
             !! The matrix is modified.
             AA(IL) = -theta*CondL*(rmhalf**2)*dt/DELSQ
             BB(IL) = -theta*CondR*(rphalf**2)*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*(rphalf**2) + CondL*(rmhalf**2))*dt/DELSQ
             
             
             S(IL) = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ)*(CondL*(rmhalf**2)*soln(IL-1,JL,KL) & 
                  -(CondL*(rmhalf**2)+CondR*(rphalf**2))*soln(IL,JL,KL) &
                  +CondR*(rphalf**2)* soln(IL+1,JL,KL))      
          end do
       else if ((diff_geometry == CARTESIAN) .or. &
            (diff_geometry == CYLINDRICAL .and. iDirection == SWEEP_Y)) then           
          do IL=2, NL-1      
             if (useCFArray) cDiv = 1.0/ConstFact(IL,JL,KL)
             
             CondL = 0.5 * (cond(IL-1,JL,KL) + cond(IL,JL,KL))*cDiv ! (IL,JL,KL)
             CondR = 0.5 * (cond(IL+1,JL,KL) + cond(IL,JL,KL))*cDiv ! (IL,JL,KL)

             DELSQ = DEL**2
             
             AA(IL) = -CondL*theta*dt/(DELSQ)
             BB(IL) = -CondR*theta*dt/(DELSQ)      
             
             DD(IL) = (1.0 + theta*(CondL+CondR)*dt/DELSQ)
             
             S(IL) = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ)* &
                     (CondL*(soln(IL-1,JL,KL)-soln(IL,JL,KL))-CondR*(soln(IL,JL,KL)-soln(IL+1,JL,KL)))
             
          end do
          
       end if
       
       select case (bcType(1))           
          
       case (GRID_PDE_BND_PERIODIC)
          T_0 = soln(NL,JL,KL)

          if (useCFArray) cDiv = 1.0/ConstFact(1,JL,KL)
          CondR =  0.5 * (cond( 1,JL,KL) + cond(2,JL,KL))*cDiv ! (1,JL,KL)
          CondL =  0.5 * (cond(NL,JL,KL) + cond(1,JL,KL))*cDiv ! (1,JL,KL)
          
          
          if  ((diff_geometry == CYLINDRICAL) .and. iDirection == SWEEP_X) then
             IL = 1
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = r*dr**2                        
             
             !! The matrix is modified.
             AA(IL) = -theta*CondL*rmhalf*dt/DELSQ
             BB(IL) = -theta*CondR*rphalf*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*rphalf + CondL*rmhalf)*dt/DELSQ - AA(IL)

             S(IL)  = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ) *(CondL*(rmhalf)*T_0 &
                  -(CondL*(rmhalf)+CondR*(rphalf))*soln(IL,JL,KL) &
                  + CondR*(rphalf)*soln(IL+1,JL,KL))    
             
             
          else        
             DELSQ = DEL**2
             AA(1) = -CondL*theta*dt/(DELSQ) 
             BB(1) = -CondR*theta*dt/(DELSQ) !not used
             DD(1) = (1.0 + theta*(CondL+CondR)*dt/DELSQ) - AA(1)

             S (1) = soln(1,JL,KL)  + dt*(1.0-theta)* &
                     (CondL*T_0 - (CondL+CondR)*soln(1,JL,KL) + CondR*soln(2,JL,KL))/DELSQ         

          end if
          
          
       case (GRID_PDE_BND_NEUMANN)
          T_0 = soln(1,JL,KL)
          
          if (useCFArray) cDiv = 1.0/ConstFact(1,JL,KL)
          CondR =  0.5 * (cond(1,JL,KL) + cond(2,JL,KL))*cDiv ! (1,JL,KL)
          CondL =  0.5 * (cond(1,JL,KL) + cond(1,JL,KL))*cDiv ! (1,JL,KL)
          
          if  ((diff_geometry == CYLINDRICAL) .and. iDirection == SWEEP_X) then
             IL = 1
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = r*dr**2    
             
             AA(IL) = -theta*CondL*rmhalf*dt/DELSQ
             BB(IL) = -theta*CondR*rphalf*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*rphalf + CondL*rmhalf)*dt/DELSQ + AA(IL)

             S(IL)  = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ) *(CondL*(rmhalf)*T_0 &
                  -(CondL*(rmhalf)+CondR*(rphalf))*soln(IL,JL,KL) &
                  + CondR*(rphalf)*soln(IL+1,JL,KL))

          else if  ((diff_geometry == SPHERICAL) .and. iDirection == SWEEP_X) then
             IL = 1
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = r*dr**2    
             DELSQ = dr*((rphalf**3) -(rmhalf**3))/3.0
             
             !! The matrix is modified.
             AA(IL) = -theta*CondL*(rmhalf**2)*dt/DELSQ
             BB(IL) = -theta*CondR*(rphalf**2)*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*(rphalf**2) + CondL*(rmhalf**2))*dt/DELSQ + AA(IL)
             
          
             S(IL) = soln(IL,JL,KL) &
                     + (1.0-theta)*(dt/DELSQ)*(CondL*(rmhalf**2)*T_0 & 
                     - (CondL*(rmhalf**2)+CondR*(rphalf**2))*soln(IL,JL,KL) &
                     + CondR*(rphalf**2)* soln(IL+1,JL,KL))   

          else                    
             DELSQ = DEL**2
             AA(1) = -CondL*theta*dt/(DELSQ)
             BB(1) = -CondR*theta*dt/(DELSQ) !not used
             DD(1) = (1.0 + theta*(CondL+CondR)*dt/DELSQ) + AA(1)
             S (1) = soln(1,JL,KL) &
                     + dt*(1.0-theta)*(CondL*T_0 - (CondL+CondR)*soln(1,JL,KL) + CondR*soln(2,JL,KL))/DELSQ      
          end if
          
          
       case(GRID_PDE_BND_DIRICHLET)
!!$          T_0 =(2.0*2000.0-soln(1,JL,KL)) !bcValues(1) ! 0.0
          
          T_0 = bcValues(1)
          
          if (useCFArray) cDiv = 1.0/ConstFact(1,JL,KL)
          CondR =  0.5 * (cond(1,JL,KL) + cond(2,JL,KL))*cDiv ! (1,JL,KL)
          CondL =  0.5 * (cond(1,JL,KL) + cond(1,JL,KL))*cDiv ! for now - should use cond(0,JL,KL) if available

          if  ((diff_geometry == CYLINDRICAL) .and. iDirection == SWEEP_X) then

             IL = 1
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = r*dr**2    
             
             AA(IL) = -theta*CondL*rmhalf*dt/DELSQ
             BB(IL) = -theta*CondR*rphalf*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*rphalf + CondL*rmhalf)*dt/DELSQ
             
             S(IL)  = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ) *(CondL*(rmhalf)*T_0 &
                  -(CondL*(rmhalf)+CondR*(rphalf))*soln(IL,JL,KL) &
                  + CondR*(rphalf)*soln(IL+1,JL,KL)) - AA(IL)*T_0
             
             
          else if  ((diff_geometry == SPHERICAL) .and. iDirection == SWEEP_X) then
             
             IL = 1
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = r*dr**2    
             DELSQ = dr*((rphalf**3) -(rmhalf**3))/3.0

             
             !! The matrix is modified.
             AA(IL) = -theta*CondL*(rmhalf**2)*dt/DELSQ
             BB(IL) = -theta*CondR*(rphalf**2)*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*(rphalf**2) + CondL*(rmhalf**2))*dt/DELSQ
             
             
             S(IL) = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ)*(CondL*(rmhalf**2)*T_0 & 
                  -(CondL*(rmhalf**2)+CondR*(rphalf**2))*soln(IL,JL,KL) &
                  +CondR*(rphalf**2)* soln(IL+1,JL,KL)) - AA(IL)*T_0
             

          else 
             AA(1) = -CondL*theta*dt/(DELSQ)
             BB(1) = -CondR*theta*dt/(DELSQ) !not used
             DD(1) = (1.0 + theta*(CondL+CondR)*dt/DELSQ)
             S (1) = soln(1,JL,KL) &
                     + dt*(1.0-theta)*(CondL*T_0 - (CondL+CondR)*soln(1,JL,KL) + CondR*soln(2,JL,KL))/DELSQ - AA(1)*T_0
          end if

          
       case (MARSHAK)

          
!!$          ! For Su-Olson (1996) problem
!!$          
!!$          FInc =diff_asol*(1.0E3*1.16045221E4) ** 4 ! TRAD=1.0E3 eV (BC), such that FInc = aT**4
!!$          
!!$          ! 0.5*(BC_L + soln(1,JL,KL)) - (2/3*K)*(BC_L - soln(1,JL,KL))/DEL = FInc - Implemented Boundary condition. K=1
!!$          bcL = 0.0  ! BC at left  side of domain.
!!$          
!!$          ! BC_L has an implicit part &  explicit part.
!!$          bcL = (FInc -soln(1,JL,KL)*(0.5 - (2.0/3.0)*(1.0/(1.0*DEL)))) / (0.5 + (2.0/3.0)*(1.0/(1.0*DEL)))
!!$
!!$          bcLE = FInc  / (0.5 + (2.0/3.0)*(1.0/(1.0*DEL)))
!!$          bcLI = - ((0.5 - (2.0/3.0)*(1.0/(1.0*DEL)))/(0.5 + (2.0/3.0)*(1.0/(1.0*DEL))))
!!$
!!$          condL = chi
!!$          condR = chi
!!$
!!$          AA(1) = -CondL*theta*dt/(DELSQ)
!!$          BB(1) = -CondR*theta*dt/(DELSQ)
!!$          DD(1) = (1.0 + theta*(CondL+CondR)*dt/DELSQ) + AA(1)*bcLI
!!$          S (1) = soln(1,JL,KL)  + dt*(1.0-theta)*(CondL*bcL - (CondL+CondR)*soln(1,JL,KL) + CondR*soln(2,JL,KL))/DELSQ - AA(1)*bcLE

       case (VACUUM)                
          

          if (useCFArray) cDiv = 1.0/ConstFact(1,JL,KL)

          CondR =  0.5 * (cond(1,JL,KL) + cond(2,JL,KL))*cDiv ! (NL,JL,KL)
          CondL =  0.5 * (cond(1,JL,KL) + cond(1,JL,KL))*cDiv ! (NL,JL,KL)

          coeff = (2.0*cond(1,JL,KL)/(diff_speedlt*DEL))
          T_0 = soln(1,JL,KL)*((coeff - 0.5)/(coeff + 0.5))
          
          if  ((diff_geometry == CYLINDRICAL) .and. iDirection == SWEEP_X) then

             
             IL = 1
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = r*dr**2              

             
             !! The matrix is modified.
             AA(IL) = -theta*CondL*rmhalf*dt/DELSQ
             BB(IL) = -theta*CondR*rphalf*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*rphalf + CondL*rmhalf)*dt/DELSQ  + AA(IL)*((coeff - 0.5)/(coeff + 0.5))
             
             S(IL)  = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ) *(CondL*(rmhalf)*T_0 &
                  -(CondL*(rmhalf)+CondR*(rphalf))*soln(IL,JL,KL) &
                  + CondR*(rphalf)*soln(IL+1,JL,KL))    

             AA(IL) = 0.
          else                              
             
             AA(1) = -CondL*theta*dt/(DELSQ)
             BB(1) = -CondR*theta*dt/(DELSQ)
             
             !! Modify diagnol to handle implicit part
             DD(1) = (1.0 + theta*(CondL+CondR)*dt/DELSQ) + AA(1)*((coeff - 0.5)/(coeff + 0.5))

             
             
             !! use BC for RHS
             S (1) = soln(1,JL,KL) &
                     + dt*(1.0-theta)*(CondL*T_0 - (CondL+CondR)*soln(1,JL,KL) + CondR*soln(2,JL,KL))/DELSQ   
          end if
          


       end select
  
       select case (bcType(2))

       case (GRID_PDE_BND_PERIODIC)

          T_N = soln(1,JL,KL)
          if (useCFArray) cDiv = 1.0/ConstFact(NL,JL,KL)
          CondR =  0.5 * (cond(NL,JL,KL) + cond(1,JL,KL))*cDiv ! (NL,JL,KL)
          CondL =  0.5 * (cond(NL,JL,KL) + cond(NL-1,JL,KL))*cDiv ! (NL,JL,KL)

          AA(NL) = -CondL*theta*dt/(DELSQ)
          BB(NL) = -CondR*theta*dt/(DELSQ)
          DD(NL) = (1.0 + theta*(CondL+CondR)*dt/DELSQ) - BB(NL)
          S (NL) = soln(NL,JL,KL) &
                   + dt*(1.0-theta)*(CondL*soln(NL-1,JL,KL)- (CondL+CondR)*soln(NL,JL,KL) + CondR*T_N)/DELSQ

       case (GRID_PDE_BND_NEUMANN)
          T_N = soln(NL,JL,KL)
          if (useCFArray) cDiv = 1.0/ConstFact(NL,JL,KL)
          CondR =  0.5 * (cond(NL,JL,KL) + cond(NL,JL,KL))*cDiv ! (NL,JL,KL)
          CondL =  0.5 * (cond(NL,JL,KL) + cond(NL-1,JL,KL))*cDiv ! (NL,JL,KL)

          if  ((diff_geometry == CYLINDRICAL) .and. iDirection == SWEEP_X) then
             IL = NL
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = r*dr**2    
             
             AA(IL) = -theta*CondL*rmhalf*dt/DELSQ
             BB(IL) = -theta*CondR*rphalf*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*rphalf + CondL*rmhalf)*dt/DELSQ + BB(NL)
             
             S(IL)  = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ) *(CondL*(rmhalf)*soln(NL-1,JL,KL) &
                  -(CondL*(rmhalf)+CondR*(rphalf))*soln(IL,JL,KL) &
                  + CondR*(rphalf)*T_N)

          else if  ((diff_geometry == SPHERICAL) .and. iDirection == SWEEP_X) then
             IL = NL
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = dr*((rphalf**3) -(rmhalf**3))/3.0 

             !! The matrix is modified.
             AA(IL) = -theta*CondL*(rmhalf**2)*dt/DELSQ
             BB(IL) = -theta*CondR*(rphalf**2)*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*(rphalf**2) + CondL*(rmhalf**2))*dt/DELSQ + BB(NL)
             
          
             S(IL) = soln(IL,JL,KL) &
                     + (1.0-theta)*(dt/DELSQ)*(CondL*(rmhalf**2)*soln(NL-1,JL,KL) & 
                     - (CondL*(rmhalf**2)+CondR*(rphalf**2))*soln(IL,JL,KL) &
                     + CondR*(rphalf**2)* T_N) 
          else
             DELSQ = DEL**2
             AA(NL) = -CondL*theta*dt/(DELSQ)
             BB(NL) = -CondR*theta*dt/(DELSQ)
             DD(NL) = (1.0 + theta*(CondL+CondR)*dt/DELSQ) + BB(NL)
             S (NL) = soln(NL,JL,KL) &
                      + dt*(1.0-theta)*(CondL*soln(NL-1,JL,KL)- (CondL+CondR)*soln(NL,JL,KL) + CondR*T_N)/DELSQ 
          end if
           
       case (GRID_PDE_BND_DIRICHLET)     !DEV: Even more untested than other parts! - KW         

          T_N = 1000.0 !bcValues(2)
          
          if (useCFArray) cDiv = 1.0/ConstFact(NL,JL,KL)
          CondR =  0.5 * (cond(NL,JL,KL) + cond(NL,JL,KL))*cDiv ! for now - should use cond(NL+1,JL,KL) if available
          CondL =  0.5 * (cond(NL,JL,KL) + cond(NL-1,JL,KL))*cDiv ! (NL,JL,KL)

          if  ((diff_geometry == CYLINDRICAL) .and. iDirection == SWEEP_X) then

             IL = NL
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = r*dr**2    
             
             AA(IL) = -theta*CondL*rmhalf*dt/DELSQ
             BB(IL) = -theta*CondR*rphalf*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*rphalf + CondL*rmhalf)*dt/DELSQ
             
             S(IL)  = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ) *(CondL*soln(IL-1,JL,KL) &
                  -(CondL*(rmhalf)+CondR*(rphalf))*soln(IL,JL,KL) &
                  + CondR*(rphalf)*soln(IL+1,JL,KL)) - BB(NL)*T_N
             
             

          else if  ((diff_geometry == SPHERICAL) .and. iDirection == SWEEP_X) then
             IL = NL
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = dr*((rphalf**3) -(rmhalf**3))/3.0
             !! The matrix is modified.
             AA(IL) = -theta*CondL*(rmhalf**2)*dt/DELSQ
             BB(IL) = -theta*CondR*(rphalf**2)*dt/DELSQ
             DD(IL) = 1.0 + theta*(CondR*(rphalf**2) + CondL*(rmhalf**2))*dt/DELSQ
             
          
             S(IL) = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ)*(CondL*(rmhalf**2)*soln(IL-1,JL,KL) & 
                  -(CondL*(rmhalf**2)+CondR*(rphalf**2))*soln(IL,JL,KL) &
                  +CondR*(rphalf**2)* T_N) - BB(NL)*T_N

          else 
             AA(NL) = -CondL*theta*dt/(DELSQ) !not used
             BB(NL) = -CondR*theta*dt/(DELSQ)
             DD(NL) = (1.0 + theta*(CondL+CondR)*dt/DELSQ)
             S (NL) = soln(NL,JL,KL) + dt*(1.0-theta)* &
                      (CondL*soln(NL-1,JL,KL)- (CondL+CondR)*soln(NL,JL,KL) + CondR*T_N)/DELSQ &
                      - BB(NL)*T_N
          end if


       case (VACUUM)       

          if (useCFArray) cDiv = 1.0/ConstFact(NL,JL,KL)

          CondR =  0.5 * (cond(NL,JL,KL) + cond(NL,  JL,KL))*cDiv ! (NL,JL,KL)
          CondL =  0.5 * (cond(NL,JL,KL) + cond(NL-1,JL,KL))*cDiv ! (NL,JL,KL)
          
          coeff = (2.0*cond(NL,JL,KL)/(diff_speedlt*DEL))
          T_N = soln(NL,JL,KL)*(coeff-0.5)/(coeff+0.5)
        
          
          if  ((diff_geometry == CYLINDRICAL) .and. iDirection == SWEEP_X) then
             IL = NL
             r  = (IL-0.5)*DEL
             dr = DEL      
             rphalf = r + 0.5*dr
             rmhalf = r - 0.5*dr
             DELSQ = r*dr**2              
             
             !! The matrix is modified.
             AA(NL) = -theta*CondL*rmhalf*dt/DELSQ
             BB(NL) = -theta*CondR*rphalf*dt/DELSQ
             DD(NL) = (1.0 + theta*(CondR*rphalf + CondL*rmhalf)*dt/DELSQ)  + BB(NL)*((coeff - 0.5)/(coeff + 0.5))

             BB(NL) = 0.0

             
             S(IL)  = soln(IL,JL,KL) + (1.0-theta)*(dt/DELSQ) *(CondL*(rmhalf)*soln(NL-1,JL,KL) &
                  -(CondL*(rmhalf)+CondR*(rphalf))*soln(IL,JL,KL) &
                  + CondR*(rphalf)*T_N)
          else                              
             
             AA(NL) = -CondL*theta*dt/(DELSQ)
             BB(NL) = -CondR*theta*dt/(DELSQ)
             
             !! Modify diagnol to handle implicit part
             DD(NL) = (1.0 + theta*(CondL+CondR)*dt/DELSQ) + (coeff-0.5)/(coeff+0.5)*BB(NL)

             
             !! use BC for RHS
             S (NL) = soln(NL,JL,KL) &
                      + dt*(1.0-theta)*(CondL*soln(NL-1,JL,KL) - (CondL+CondR)*soln(NL,JL,KL) + CondR*T_N)/DELSQ
          end if
                    
       end select
       
       if(useMatDArray) then  
          S(1:NL) = S(1:NL) + dt*matD(1:NL,JL,KL)/NDIM
       end if
       
       if(useMatCArray) then          
          DD(1:NL) = DD(1:NL) + dt*MatC(1:NL,JL,KL)/NDIM     
!!$      DD(:) = DD(:) + theta*MatC(:,JL,KL)/NDIM
!!$      S(:) =  S(:) - (1-theta)*matC(:,JL,KL)*soln(:,JL,KL)/NDIM
       end if
       


       !*** To handle periodic BC within the framework of Thomas Algorithm ***
       U(:)  = 0.0
       
       if (bcType(1) == PERIODIC .or. bcType(2) == PERIODIC) then          
          if (bcType(1) == PERIODIC) then        
             U(1)  = AA(1)
          endif

          if (bcType(2) == PERIODIC) then
             U(NL)  = BB(NL)
          endif
          
          call Thomas3D(AA,DD,BB,U,NL)
       end if

       call Thomas3D(AA,DD,BB,S,NL)
       soln(:,JL,KL)  = S(:) - ((S(1) + S(NL))/(1.0 + U(1) + U(NL)))*U(:)

    end do
 end do

 deallocate (AA)
 deallocate (DD)
 deallocate (BB)
 deallocate (S)
 deallocate (U)

contains 
  ! *** Thomas Algorithm - for tridiag matrix ***
  subroutine Thomas3D (BB,DDL,AA,S,n)

    implicit none

    integer, intent (IN)    :: n
    real   , intent (IN) :: BB(n), DDL(n), AA(n)
    real   , intent (INOUT) :: S(n)

    integer i,j
    real    R
    real    DD(n)

    DD (:) = DDL(:)
    ! Elimination.
    do  i  = 2, n
       R     = BB(i)/DD(i-1)
       DD(i) = DD(i)-R*AA(i-1)
       S (i) = S(i)-R*S(i-1)
    end do

    ! Back Substitution.
    S (n) = S(n)/DD(n)
    do i    = 2, n
       j    = n - i + 1
       S(j) = (S(j) - AA(j)*S(j+1))/DD(j)
    end do

  end subroutine Thomas3D

end subroutine diff_advance1DImpl
