!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/RK/ed_raytraceODEfunction2DRec
!!
!! NAME
!!
!!  ed_raytraceODEfunction2DRec
!!
!! SYNOPSIS
!!
!!  ed_raytraceODEfunction2DRec (real, intent (in) :: t,
!!                               real, intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Defines the ray trace ODE function in 2D rectangular coordinates to be used for Runge Kutta
!!  ray tracing.
!!
!! ARGUMENTS
!!
!!  t    : independent variable (time)
!!  y    : dependent variable array (rx,ry,vx,vy,P)
!!
!!***

function ed_raytraceODEfunction2DRec (t,y)

  use EnergyDeposition_data,      ONLY : ed_Boltzmann,      &
                                         ed_cellCubicNele,  &
                                         ed_cellCubicTele,  &
                                         ed_electronMass,   &
                                         ed_electronCharge

  use ed_raytraceODEfunctionData, ONLY : accX, accY,                 &
                                         accFactorX, accFactorY,     &
                                         cellEdgeInvX, cellEdgeInvY, &
                                         cellZbar,                   &
                                         i, j,                       &
                                         Nele,                       &
                                         rayCritDens,                &
                                         saveComputations,           &
                                         x01, y01,                   &
                                         xminCell, yminCell

  use Driver_interface,           ONLY : Driver_abortFlash

  use Interpolate_interface,      ONLY : Interpolate_cubic2DF,  &
                                         Interpolate_cubic2DFd1

  use ed_interface,               ONLY : ed_CoulombFactor,             &
                                         ed_inverseBremsstrahlungRate

  implicit none

  real, intent (in) :: t                             ! it is absolutely mandatory to
  real, intent (in) :: y (:)                         ! declare the variables and
                                                     ! array function in the way shown.
  real :: ed_raytraceODEfunction2DRec (1:size (y))   ! (compatible to Runge Kutta interfaces)

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

       Nele (1:3) = Interpolate_cubic2DFd1 (ed_cellCubicNele (1:16,i,j,1), x01,y01)

       if (Nele (1) <= 0.0) then
           call Driver_abortFlash ('[ed_raytraceODEfunction2DRec] ERROR: Nele <= 0 for a cell')
       end if

       accX = accFactorX * Nele (2)
       accY = accFactorY * Nele (3)

  end if

  Tele = Interpolate_cubic2DF (ed_cellCubicTele (1:16,i,j,1), x01,y01)

  if (Tele <= 0.0) then
      call Driver_abortFlash ('[ed_raytraceODEfunction2DRec] ERROR: Tele <= 0 for a cell')
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
  ed_raytraceODEfunction2DRec (1) = y (3)         ! the old velX
  ed_raytraceODEfunction2DRec (2) = y (4)         ! the old velY
  ed_raytraceODEfunction2DRec (3) = accX
  ed_raytraceODEfunction2DRec (4) = accY
  ed_raytraceODEfunction2DRec (5) = - nu * y (5)  ! the old power is in y (5)
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
end function ed_raytraceODEfunction2DRec
