!!****if* source/numericalTools/RungeKutta/RungeKuttaMain/RungeKutta_stepSizeEstimate
!!
!! NAME
!!
!!  RungeKutta_stepSizeEstimate
!!
!! SYNOPSIS
!!
!!  RungeKutta_stepSizeEstimate (character (len=*), intent (in)    :: method,
!!                               function,          intent (in)    :: f,
!!                               real,              intent (in)    :: x,
!!                               real,              intent (in)    :: y     (:),
!!                               real,              intent (in)    :: eFrac,
!!                               real,              intent (in)    :: eBase (:),
!!                               real,              intent (in)    :: hmax)
!!
!! DESCRIPTION
!!
!!  This function provides the user with an initial estimate in step size, based on the
!!  error vector supplied and the ODE function 'f' to be solved. The estimate is based on
!!  all local information at 'x' and is obtained as follows. Consider the finite Taylor
!!  expansion of y (x+h) around 'x':
!!
!!           y(x+h) = y(x) + hy'(x) + (h^2/2)y''(x) + ... + (h^p/p!)y^p(x) + O[h^(p+1)]
!!
!!  The truncation error O[h^(p+1)] can be expressed exactly as:
!!
!!           O[h^(p+1)] = h^(p+1)/(p+1)! * y^(p+1)(x*)   x* between x and x+h
!!
!!  where the neglected higher order terms of the Taylor expansion have been accounted for
!!  by the definition of x*. Taking only the leading error term E(p+!) in O[h^(p+1)] is
!!  equivalent of setting x*=x in the above equation:
!!
!!              E(p+1) = h^(p+1)/(p+1)! * y^(p+1)(x)
!!
!!  A Runge Kutta step of order p is like a Taylor expansion of order p and leads also to
!!  an error term O[h^(p+1)]. Hence to estimate the step size of a Runge Kutta of order p
!!  one could simply use the expression for E(p+1) and equate it with the error goal of the
!!  Runge Kutta step:
!!
!!               Runge Kutta error goal <= eFrac * eBase
!!
!!  and solve for h in the equation:
!!
!!                      eFrac * eBase >= E(p+1)
!!                      eFrac * eBase >= h^(p+1)/(p+1)! * y^(p+1)(x)
!!
!!  Since eBase and y^(p+1)(x) are vectors, the equation is actually a set of equations,
!!  one for each of the corresponding vector elements. This results in a set of different
!!  h values, of which we will take the smallest:
!!
!!                      h <= [eFrac * (p+1)! * min |eBase/y^(p+1)(x)|]^(1/(p+1))
!!  or
!!                      h <= [eFrac * (p+1)! / max |y^(p+1)(x)/eBase|]^(1/(p+1))
!!
!!  One of the difficulties with the above equation is that for high orders it becomes difficult
!!  to estimate the higher order derivative vector y^(p+1)(x). To circumvent this problem and
!!  to use lower order derivatives, it is posulated that the eBase-rescaled higher order Taylor
!!  expansion errors decay in form of a power law:
!!
!!                  E(p+1)/eBase = [E(m)/eBase]^((p+1)/m)       p+1 >= m
!!
!!  The eBase vector contains (or should contain) implicitly information about the size of
!!  the ODE's dependent variables and should never contain zero values. eFrac relates to
!!  the relative accuracy wanted for each dependent variable along the solution path of the
!!  ODE. Hence eFrac should always be < 1. The step size estimate for a Runge Kutta of order
!!  p >= m starts from the equation:
!!
!!                      eFrac * eBase >= E(p+1)
!!
!!  and uses the above postulate for E(p+1) to get:
!!
!!                h <= eFrac^(1/(p+1)) * [m! / max |y^(m)(x)/eBase|]^(1/m)
!!
!!  For calculating our step size estimate for a given Runge Kutta order p, we use the m = 1
!!  and m = 2 error estimate formulas. The following steps are taken:
!!
!!                1) Is max |y''(x)| > 0? If yes, use:
!!
!!                       h = eFrac^(1/(p+1)) * sqrt (2) / sqrt (max |y''(x)/eBase|)
!!                       h = min (h, hmax)
!!
!!                   for orders p >= 1. If no, go to step 2).
!!
!!                2) Is max |y'(x)| > 0? If yes, use:
!!
!!                       h = eFrac^(1/(p+1)) / max |y'(x)/eBase|
!!                       h = min (h, hmax)
!!
!!                   for orders p >= 1. If no, goto step 3).
!!
!!                3) Use:
!!
!!                       h = hmax
!!
!!                   for all orders.
!!
!! ARGUMENTS
!!
!!  method       : the type of RK method for which the estimate is
!!  f            : the ODE function containing details of the ODE system to be solved
!!  x            : the initial value of the independent variable
!!  y            : the initial values of the dependent variable(s)
!!  eFrac        : the fractional value for each of the error base values
!!  eBase        : the error base values for each of the dependent variables
!!  hmax         : the safeguard maximum allowed step size
!!
!! NOTES
!!
!!  Contains an interfaced function in the declarations (see below).
!!
!!***

real function RungeKutta_stepSizeEstimate (method,f,x,y,eFrac,eBase,hmax)

  use rk_interface,      ONLY : rk_orderRKmethod

  use Driver_interface,  ONLY : Driver_abortFlash

  use RungeKutta_data,   ONLY : rk_cubeRootMacheps

  implicit none

  interface
    function f (x,y)
      real, intent (in) :: x
      real, intent (in) :: y (:)
      real              :: f (1:size (y))
    end function f
  end interface

  character (len=*), intent (in)  :: method
  real,              intent (in)  :: x
  real,              intent (in)  :: y     (:)
  real,              intent (in)  :: eFrac
  real,              intent (in)  :: eBase (:)
  real,              intent (in)  :: hmax

  integer :: p
  real    :: d,h
!
!
!   ...Check, if eFrac is < 1 and all eBase values are different from zero.
!
!
  if (eFrac >= 1.0) then
      call Driver_abortFlash ('[RungeKutta_stepSizeEstimate] ERROR: eFrac >= 1!')
  end if

  if (any (abs (eBase) == 0.0)) then
      call Driver_abortFlash ('[RungeKutta_stepSizeEstimate] ERROR: eBase has zero values!')
  end if
!
!
!   ...Check, if y (:) and eBase (:) have the same number of elements.
!
!
  if (size (y) /= size (eBase)) then
      call Driver_abortFlash ('[RungeKutta_stepSizeEstimate] ERROR: y/eBase size mismatch!')
  end if
!
!
!   ...Calculate the step size estimate.
!
!
  p = rk_orderRKmethod (method)

  h = max (abs (x), 1.0) * rk_cubeRootMacheps    ! estimate optimum numerical derivative step

  d = maxval (abs ( ((f (x + h, y (:) + h * f (x, y (:))) &
                    - f (x - h, y (:) - h * f (x, y (:)))) / (h + h)) / eBase (:)))

  if (d /= 0.0) then

      h = sqrt (2.0 / d) * eFrac ** (1.0 / real (p + 1))     ! 2nd derivative estimate
      RungeKutta_stepSizeEstimate = min (h, hmax)

  else
      d = maxval (abs (f (x, y (:)) / eBase (:)))

      if (d /= 0.0) then
          h = (1.0 / d) * eFrac ** (1.0 / real (p + 1))      ! 1st derivative estimate
          RungeKutta_stepSizeEstimate = min (h, hmax)
      else
          RungeKutta_stepSizeEstimate = hmax
      end if

  end if
!
!
!   ...Ready! 
!
!
  return
end function RungeKutta_stepSizeEstimate
