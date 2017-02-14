!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/2Dellipse/sim_RungeKuttaTest
!!
!! NAME
!!
!!  sim_RungeKuttaTest
!! 
!! SYNOPSIS
!!
!!  call sim_RungeKuttaTest
!!
!! DESCRIPTION
!!
!!  This function tests the Runge Kutta stepper of the Runge Kutta utility unit.
!!  A ODE function is defined for this test and the stepper is called many times
!!  recording advancement of the function + error.
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_RungeKuttaTest ()

  use  Simulation_data,      ONLY: sim_ellipseAspectRatio,   &
                                   sim_ellipseMajorSemiAxis, &
                                   sim_ellipseMinorSemiAxis, &
                                   sim_ellipseRotationAngle, &
                                   sim_ellipseCenterX,       &
                                   sim_ellipseCenterY,       &
                                   sim_stepSize,             &
                                   sim_errorFraction,        &
                                   sim_numberOfEllipses,     &
                                   sim_RungeKuttaMethod,     &
                                   sim_x0,                   &
                                   sim_y0


  use  RungeKutta_interface, ONLY: RungeKutta_Step,             &
                                   RungeKutta_StepSizeEstimate

  use  sim_interface,        ONLY: sim_distancePoint2Ellipse2Dxy, &
                                   sim_ODEfunction

  implicit none

  character (len=14), parameter :: rTitle   = '    radius    '
  character (len=14), parameter :: xTitle   = '      x       '
  character (len=14), parameter :: yTitle   = '      y       '
  character (len=14), parameter :: xtTitle  = '     x(t)     '
  character (len=14), parameter :: ytTitle  = '     y(t)     '
  character (len=14), parameter :: exTitle  = '     err x    '
  character (len=14), parameter :: eyTitle  = '     err y    '
  character (len=14), parameter :: erTitle  = '     err r    '
  character (len=14), parameter :: epTitle  = '     err p    '
  character (len=14), parameter :: hTitle   = '    h used    '
  character (len=14), parameter :: tTitle   = '       t      '
  character (len=14), parameter :: nHETitle = ' Half-Ellipse '

  integer :: n
  integer :: nStep
  integer :: nHalfEllipses
  integer :: nStepPerEllipse

  real    :: A
  real    :: ex, ey, er, ep, ermax
  real    :: erAvgPerEllipse, erSumPerEllipse
  real    :: htry, hused, hnext, hmin, hmax
  real    :: k1, k2, k3
  real    :: phi, sinPhi, cosPhi
  real    :: t
  real    :: x, y, r
  real    :: xt, yt
  real    :: yprev

  real    :: d     (1:2)
  real    :: p     (1:2)
  real    :: pnew  (1:2)
  real    :: error (1:2)
  real    :: eBase (1:2)
!
!
!   ...Set error base for all dependent variables.
!
!  
  eBase (1) = 1.0
  eBase (2) = 1.0
!
!
!   ...Initialize dependent variables.
!
!  
  p (1) = sim_x0
  p (2) = sim_y0

  x = p (1)
  y = p (2)
  r = sqrt (x * x + y * y)

  xt = sim_x0
  yt = sim_y0

  A  = sim_ellipseAspectRatio
  k1 = (A + A) / (A * A - 1.0)         ! this is sqrt (k^2 - 1)
  k2 = (A * A - 1.0) / (A + A)         ! this is 1/ sqrt (k^2 - 1)
  k3 = (A * A + 1.0) / (A + A)         ! this is k/ sqrt (k^2 - 1)

  ex = 0.0
  ey = 0.0
  er = 0.0
  ep = 0.0

  write (*,'(12a14)') rTitle, xTitle, yTitle, xtTitle, ytTitle, &
                      exTitle, eyTitle, erTitle, epTitle,       &
                      hTitle, tTitle, nHETitle
  write (*,'(150a1)') ('-',n = 1,150)

  write (*,'(f10.4,4x,4(f12.6,2x),6(es10.3,4x),i9,5x)') r,x,y,xt,yt,ex,ey,er,ep
!  write (*,'(f9.4,2x,f9.4,2x,es10.3,2x,es10.3,2x,es10.3)') x,y
!
!
!   ...Call the Runge Kutta stepper.
!
!  
  t = 0.0
  nStep = 0

  htry = RungeKutta_StepSizeEstimate (sim_RungeKuttaMethod, &
                                      sim_ODEfunction,      &
                                      t,                    &
                                      p,                    &
                                      sim_errorFraction,    &
                                      eBase,                &
                                      hmax = 1000.0         )

  write (*,*) ' h estimate = ',htry

  hmin = htry
  hmax = htry
  ermax = 0.0

  nHalfEllipses   = 0
  nStepPerEllipse = 0
  erSumPerEllipse = 0.0

  do while (nHalfEllipses < 2 * sim_numberOfEllipses)

     call RungeKutta_Step (sim_RungeKuttaMethod,   &
                           sim_ODEfunction,        &
                           t,                      &
                           p,                      &
                           sim_errorFraction,      &
                           eBase,                  &
                                            htry,  &
                                            hused, &
                                            hnext, &
                                            pnew,  &
                                            error  )


     nStep = nStep + 1
!     nStepPerEllipse = nStepPerEllipse + 1

     p (:) = pnew (:)

     yprev = y

     x = p (1)
     y = p (2)
     r = sqrt (x * x + y * y)

     t = t + hused
     hmin = min (hmin, hused)
     hmax = max (hmax, hused)

     phi = t * k1
     sinPhi = sin (phi)
     cosPhi = cos (phi)

     xt = sim_x0 * (cosPhi + sinPhi * k2) + sim_y0 * sinPhi * k3
     yt = sim_y0 * (cosPhi - sinPhi * k2) - sim_x0 * sinPhi * k3

     d (1:2) = sim_distancePoint2Ellipse2Dxy (sim_ellipseMajorSemiAxis, &
                                              sim_ellipseMinorSemiAxis, &
                                              sim_ellipseRotationAngle, &
                                              sim_ellipseCenterX,       &
                                              sim_ellipseCenterY,       &
                                              x, y                      )

     ex = abs (error (1))
     ey = abs (error (2))
     er = d (1)
     ermax = max (ermax, er)
     ep = sqrt ((x - xt) * (x - xt) + (y - yt) * (y - yt))
!
!     erSumPerEllipse = erSumPerEllipse + er

     if (abs (yprev - y) > max (abs (yprev - sim_y0), abs (y - sim_y0))) then    ! this records sign changes when
         nHalfEllipses = nHalfEllipses + 1                                       ! crossing the y = sim_y0 line

!         if (mod (nHalfEllipses,2) == 0) then
!             erAvgPerEllipse = erSumPerEllipse / real (nStepPerEllipse)
!             erSumPerEllipse = er
!             nStepPerEllipse = 1
!             write (*,'(i8,2x,es10.3)') nHalfEllipses/2 , erAvgPerEllipse
!         end if

     end if

     write (*,'(f10.4,4x,4(f12.6,2x),6(es10.3,4x),i9,5x)') r,x,y,xt,yt,ex,ey,er,ep,hused,t,nHalfEllipses
!     write (*,'(f9.4,2x,f9.4,2x,es10.3,2x,es10.3,2x,es10.3)') x,y
!
!     write (*,'(a,f6.2,a,f6.2,a,4(f8.3,a),4(es8.1,a))') ' $',r, '$ & $',t, '$ & $', &
!                                                             x, '$ & $',y, '$ & $',xt,'$ & $',yt,'$ & $', &
!                                                             ex,'$ & $',ey,'$ & $',ep,'$ & $',er,'$ \\'

     htry = hnext

  end do

  write (*,*) ' Number of steps = ', nStep
  write (*,*) ' Avg stepsize    = ', t / real (nStep)
  write (*,*) ' Min stepsize    = ', hmin
  write (*,*) ' Max stepsize    = ', hmax
  write (*,*) ' Max er          = ', ermax

!
!
!   ...Ready!
!
!
  return
end subroutine sim_RungeKuttaTest
