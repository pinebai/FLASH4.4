!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/RK/ed_raytraceODEfunction2DCyl3D
!!
!! NAME
!!
!!  ed_raytraceODEfunction2DCyl3D
!!
!! SYNOPSIS
!!
!!  ed_raytraceODEfunction2DCyl3D (real, intent (in) :: t,
!!                                 real, intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Defines the ray trace ODE function for the Runge Kutta unit, to perform 3D ray tracing in
!!  2D cylindrical coordinates.
!!
!! ARGUMENTS
!!
!!  t    : independent variable (time)
!!  y    : dependent variable array (rx,ry,rz,vx,vy,vz,P)
!!
!!***

function ed_raytraceODEfunction2DCyl3D (t,y)

  use EnergyDeposition_data,      ONLY : ed_Boltzmann,      &
                                         ed_cellCubicNele,  &
                                         ed_cellCubicTele,  &
                                         ed_electronMass,   &
                                         ed_electronCharge

  use ed_raytraceODEfunctionData, ONLY : accX, accZ,                   &
                                         accFactorX, accFactorZ,       &
                                         wedgeEdgeInvX, wedgeEdgeInvZ, &
                                         wedgeZbar,                    &
                                         i, j,                         &
                                         Nele,                         &
                                         rayCritDens,                  &
                                         saveComputations,             &
                                         x01, z01,                     &
                                         xminWedge, zminWedge

  use Driver_interface,           ONLY : Driver_abortFlash

  use Interpolate_interface,      ONLY : Interpolate_cubic2DF,  &
                                         Interpolate_cubic2DFd1

  use ed_interface,               ONLY : ed_CoulombFactor,             &
                                         ed_inverseBremsstrahlungRate

  implicit none

  real, intent (in) :: t                             ! it is absolutely mandatory to
  real, intent (in) :: y (:)                         ! declare the variables and
                                                     ! array function in the way shown.
  real :: ed_raytraceODEfunction2DCyl3D (1:size (y)) ! (compatible to Runge Kutta interfaces)

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

       x01 = (y (1) - xminWedge) * wedgeEdgeInvX     ! rescaled [0,1] ray x coordinate
       z01 = (y (3) - zminWedge) * wedgeEdgeInvZ     ! rescaled [0,1] ray z coordinate

       Nele (1:3) = Interpolate_cubic2DFd1 (ed_cellCubicNele (1:16,i,j,1), x01,z01)

       if (Nele (1) <= 0.0) then
           call Driver_abortFlash ('[ed_raytraceODEfunction2DCyl3D] ERROR: Nele <= 0 for a cell')
       end if

       accX = accFactorX * Nele (2)
       accZ = accFactorZ * Nele (3)

  end if

  Tele = Interpolate_cubic2DF (ed_cellCubicTele (1:16,i,j,1), x01,z01)

  if (Tele <= 0.0) then
      call Driver_abortFlash ('[ed_raytraceODEfunction2DCyl3D] ERROR: Tele <= 0 for a cell')
  end if

  lnLambda = ed_CoulombFactor (wedgeZbar,         &
                               ed_electronCharge, &
                               ed_Boltzmann,      &
                               Tele,              &
                               Nele (1)           )

  nu = ed_inverseBremsstrahlungRate (wedgeZbar,         &
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
  ed_raytraceODEfunction2DCyl3D (1) = y (4)         ! the old velX
  ed_raytraceODEfunction2DCyl3D (2) = y (5)         ! the old velY
  ed_raytraceODEfunction2DCyl3D (3) = y (6)         ! the old velZ
  ed_raytraceODEfunction2DCyl3D (4) = accX
  ed_raytraceODEfunction2DCyl3D (5) = 0.0           ! no acceleration along the angular part
  ed_raytraceODEfunction2DCyl3D (6) = accZ
  ed_raytraceODEfunction2DCyl3D (7) = - nu * y (7)  ! the old power is in y (7)
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
end function ed_raytraceODEfunction2DCyl3D
