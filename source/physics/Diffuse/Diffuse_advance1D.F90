!!****f* source/physics/Diffuse/Diffuse_advance1D
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
!!                         DIRICHLET,PNEUMAN (untested).
!!
!!   del        : dx,dy or dz based on direction of solve
!!
!!   dt         : timestep.
!!
!!   theta      : 0.5 CN Scheme, 1.0 Fully Implicit, 0.0 Explicit
!!
!!   chi        : Diffusion coefficient, not used anymore.
!!
!!   cond       : Spatial Conductivity variation (K)
!!
!!   soln       : discretised  variation  of variable to be diffused (input). On Output contains the diffused solution.
!!
!!   useCFArray : Flag to determine if we have ConstFact(:,:,:)
!!
!!   ConstFact  : Spatial distributio of a factor (think rho*cv in heat conduction). !!
!!
!! SEE ALSO
!!
!!  Grid_advanceDiffusion
!!  Diffuse
!!
!! NOTES
!!
!!
!!  This public API subroutine is so far used only from
!!  the pencil grid implementation of Grid_advanceDiffusion.
!!  The pencil grid implementation of Grid_advanceDiffusion, in turn,
!!  is currently only used by the Split implementation of the DiffuseMain
!!  subunit.
!!***


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
  real   ,OPTIONAL,       pointer       :: ConstFact(:,:,:)
  logical,                intent(IN)    :: useMatCArray    
  real   ,OPTIONAL,       pointer       :: MatC(:,:,:)
  logical,                intent(IN)    :: useMatDArray
  real   ,OPTIONAL,       pointer       :: MatD(:,:,:)
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
end subroutine Diffuse_advance1D_1DArr
