!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/RK/ed_raytraceODEfunction1DRec
!!
!! NAME
!!
!!  ed_raytraceODEfunction1DRec
!!
!! SYNOPSIS
!!
!!  ed_raytraceODEfunction1DRec (real, intent (in) :: t,
!!                               real, intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Defines the ray trace ODE function in 1D rectangular coordinates to be used for Runge Kutta
!!  ray tracing.
!!
!! ARGUMENTS
!!
!!  t    : independent variable (time)
!!  y    : dependent variable array (rx,vx,P)
!!
!!***

function ed_raytraceODEfunction1DRec (t,y)

  use EnergyDeposition_data,      ONLY : ed_Boltzmann,      &
                                         ed_cellCubicNele,  &
                                         ed_cellCubicTele,  &
                                         ed_electronMass,   &
                                         ed_electronCharge

  use ed_raytraceODEfunctionData, ONLY : accX,             &
                                         accFactorX,       &
                                         cellEdgeInvX,     &
                                         cellZbar,         &
                                         i,                &
                                         Nele,             &
                                         rayCritDens,      &
                                         saveComputations, &
                                         x01,              &
                                         xminCell

  use Driver_interface,           ONLY : Driver_abortFlash

  use Interpolate_interface,      ONLY : Interpolate_cubic1DF,  &
                                         Interpolate_cubic1DFd1

  use ed_interface,               ONLY : ed_CoulombFactor,             &
                                         ed_inverseBremsstrahlungRate

  implicit none

  real, intent (in) :: t                             ! it is absolutely mandatory to
  real, intent (in) :: y (:)                         ! declare the variables and
                                                     ! array function in the way shown.
  real :: ed_raytraceODEfunction1DRec (1:size (y))   ! (compatible to Runge Kutta interfaces)

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

       Nele (1:2) = Interpolate_cubic1DFd1 (ed_cellCubicNele (1:4,i,1,1), x01)

       if (Nele (1) <= 0.0) then
           call Driver_abortFlash ('[ed_raytraceODEfunction1DRec] ERROR: Nele <= 0 for a cell')
       end if

       accX = accFactorX * Nele (2)

  end if

  Tele = Interpolate_cubic1DF (ed_cellCubicTele (1:4,i,1,1), x01)

  if (Tele <= 0.0) then
      call Driver_abortFlash ('[ed_raytraceODEfunction1DRec] ERROR: Tele <= 0 for a cell')
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
  ed_raytraceODEfunction1DRec (1) = y (2)         ! the old velX
  ed_raytraceODEfunction1DRec (2) = accX
  ed_raytraceODEfunction1DRec (3) = - nu * y (3)  ! the old power is in y (3)
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
end function ed_raytraceODEfunction1DRec
