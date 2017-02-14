!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonsTrace/pi_traceODEfunction3DRec
!!
!! NAME
!!
!!  pi_traceODEfunction3DRec
!!
!! SYNOPSIS
!!
!!  pi_traceODEfunction3DRec (real, intent (in) :: t,
!!                            real, intent (in) :: y (:))
!!
!! DESCRIPTION
!!
!!  Defines the proton trace ODE function in 3D rectangular coordinates to be used for Runge Kutta
!!  proton tracing.
!!
!! ARGUMENTS
!!
!!  t    : independent variable (time)
!!  y    : dependent variable array (px,py,pz,vx,vy,vz,Jv,Kx,Ky,Kz   position,velocity,diagnostics)
!!
!!***

function pi_traceODEfunction3DRec (t,y)

  use ProtonImaging_data,      ONLY : pi_screenProtonDiagnostics

  use pi_traceODEfunctionData, ONLY : Bx, By, Bz,             &
                                      BxInvC, ByInvC, BzInvC, &
                                      CurlBx, CurlBy, CurlBz, &
                                      Ex, Ey, Ez,             &
                                      Qm

  implicit none

  real, intent (in) :: t                             ! it is absolutely mandatory to
  real, intent (in) :: y (:)                         ! declare the variables and
                                                     ! array function in the way shown.
  real :: pi_traceODEfunction3DRec (1:size (y))      ! (compatible to Runge Kutta interfaces)

  real :: vx, vy, vz, vMag
!
!
!     ...Calculate the new ODE vector. The acceleration components are calculated
!        using the Lorentz force.
!
!
  vx = y (4)         ! the old x-component of the velocity
  vy = y (5)         ! the old y-component of the velocity
  vz = y (6)         ! the old z-component of the velocity

  pi_traceODEfunction3DRec (1) = vx
  pi_traceODEfunction3DRec (2) = vy
  pi_traceODEfunction3DRec (3) = vz
  pi_traceODEfunction3DRec (4) = Qm * (Ex + vy * BzInvC - vz * ByInvC)
  pi_traceODEfunction3DRec (5) = Qm * (Ey + vz * BxInvC - vx * BzInvC)
  pi_traceODEfunction3DRec (6) = Qm * (Ez + vx * ByInvC - vy * BxInvC)
!
!
!     ...If requested, calculate additional diagnostic values.
!
!
  if (pi_screenProtonDiagnostics) then

      vMag = sqrt (vx * vx + vy * vy + vz * vz)

      pi_traceODEfunction3DRec (7)  = vx * CurlBx + vy * CurlBy + vz * CurlBz
      pi_traceODEfunction3DRec (8)  = vMag * Bx
      pi_traceODEfunction3DRec (9)  = vMag * By
      pi_traceODEfunction3DRec (10) = vMag * Bz

  end if
!
!
!     ...Ready!
!
!
  return
end function pi_traceODEfunction3DRec
