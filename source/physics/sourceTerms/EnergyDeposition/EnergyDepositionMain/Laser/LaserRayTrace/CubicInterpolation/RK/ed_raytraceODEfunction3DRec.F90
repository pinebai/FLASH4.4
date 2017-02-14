!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/RK/ed_raytraceODEfunction3DRec
!!
!! NAME
!!
!!  ed_raytraceODEfunction3DRec
!!
!! SYNOPSIS
!!
!!  ed_raytraceODEfunction3DRec (real, intent (in) :: t,
!!                               real, intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Defines the ray trace ODE function in 3D rectangular coordinates to be used for Runge Kutta
!!  ray tracing.
!!
!! ARGUMENTS
!!
!!  t    : independent variable (time)
!!  y    : dependent variable array (rx,ry,rz,vx,vy,vz,P)
!!
!!***

function ed_raytraceODEfunction3DRec (t,y)

  use EnergyDeposition_data,      ONLY : ed_Boltzmann,      &
                                         ed_cellCubicNele,  &
                                         ed_cellCubicTele,  &
                                         ed_electronMass,   &
                                         ed_electronCharge

  use ed_raytraceODEfunctionData, ONLY : accX, accY, accZ,                         &
                                         accFactorX, accFactorY, accFactorZ,       &
                                         cellEdgeInvX, cellEdgeInvY, cellEdgeInvZ, &
                                         cellZbar,                                 &
                                         i, j, k,                                  &
                                         Nele,                                     &
                                         rayCritDens,                              &
                                         saveComputations,                         &
                                         x01, y01, z01,                            &
                                         xminCell, yminCell, zminCell

  use Driver_interface,           ONLY : Driver_abortFlash

  use Interpolate_interface,      ONLY : Interpolate_cubic3DF,  &
                                         Interpolate_cubic3DFd1

  use ed_interface,               ONLY : ed_CoulombFactor,             &
                                         ed_inverseBremsstrahlungRate

  implicit none

  real, intent (in) :: t                             ! it is absolutely mandatory to
  real, intent (in) :: y (:)                         ! declare the variables and
                                                     ! array function in the way shown.
  real :: ed_raytraceODEfunction3DRec (1:size (y))   ! (compatible to Runge Kutta interfaces)

  real :: lnLambda
  real :: nu
  real :: Tele
!
!
!     ...Calculate the rescaled coordinates and acceleration (if needed) and the power rate
!        loss (- nu * P) at the current ray position.
!
!
  if (.not. saveComputations) then

       x01 = (y (1) - xminCell) * cellEdgeInvX     ! rescaled [0,1] ray x coordinate
       y01 = (y (2) - yminCell) * cellEdgeInvY     ! rescaled [0,1] ray y coordinate
       z01 = (y (3) - zminCell) * cellEdgeInvZ     ! rescaled [0,1] ray z coordinate

       Nele (1:4) = Interpolate_cubic3DFd1 (ed_cellCubicNele (1:64,i,j,k), x01,y01,z01)

       if (Nele (1) <= 0.0) then
           call Driver_abortFlash ('[ed_raytraceODEfunction3DRec] ERROR: Nele <= 0 for a cell')
       end if

       accX = accFactorX * Nele (2)
       accY = accFactorY * Nele (3)
       accZ = accFactorZ * Nele (4)

  end if

  Tele = Interpolate_cubic3DF (ed_cellCubicTele (1:64,i,j,k), x01,y01,z01)

  if (Tele <= 0.0) then
      call Driver_abortFlash ('[ed_raytraceODEfunction3DRec] ERROR: Tele <= 0 for a cell')
  end if

  lnLambda = ed_CoulombFactor (cellZbar,          &
                               ed_electronCharge, &
                               ed_Boltzmann,      &
                               Tele,              &
                               Nele (1)           )

  nu = ed_inverseBremsstrahlungRate (cellZbar,          &
                                     ed_electronCharge, &
                                     ed_electronMass,   &
                                     ed_Boltzmann,      &
                                     Tele,              &
                                     Nele (1),          &
                                     rayCritDens,       &
                                     lnLambda           )
!
!
!     ...Calculate the new ODE vector.
!
!
  ed_raytraceODEfunction3DRec (1) = y (4)         ! the old velX
  ed_raytraceODEfunction3DRec (2) = y (5)         ! the old velY
  ed_raytraceODEfunction3DRec (3) = y (6)         ! the old velZ
  ed_raytraceODEfunction3DRec (4) = accX
  ed_raytraceODEfunction3DRec (5) = accY
  ed_raytraceODEfunction3DRec (6) = accZ
  ed_raytraceODEfunction3DRec (7) = - nu * y (7)  ! the old power is in y (7)
  if (size(y) .GE. 8) &
       ed_raytraceODEfunction3DRec (8) =        y (7) 
!
!
!     ...reset the saving computations keyword.
!
!
  saveComputations = .false.
!
!
!     ...Ready!
!
!
  return
end function ed_raytraceODEfunction3DRec
