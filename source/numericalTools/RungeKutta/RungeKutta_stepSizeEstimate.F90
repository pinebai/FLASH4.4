!!****f* source/numericalTools/RungeKutta/RungeKutta_stepSizeEstimate
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
!!                               real, optional,    intent (in)    :: hmax)
!!
!! DESCRIPTION
!!
!!  This function provides the user with an initial estimate in step size, based on the
!!  error vector supplied and the ODE function 'f' to be solved. The estimate is based on
!!  all local information at 'x' and is obtained as follows. Consider the Taylor expansion
!!  of y (x+h) around 'x':
!!
!!           y(x+h) = y(x) + hy'(x) + (h^2/2)y''(x) + ... + (h^m/m!)y^m(x)
!!
!!  Its truncation error can be expressed as:
!!
!!             E(m) = h^(m+1)/(m+1)! * y^(m+1)(x*)      x* between x and x+h
!!
!!  For smooth functions we can approximate the (m+1)-th derivative at x* by the
!!  (m+1)-th derivative at x, and the error becomes:
!!
!!             E(m) = h^(m+1)/(m+1)! * y^(m+1)(x)
!!
!!  We now transform the error based on the 'y' vector to an error based on eFrac. The
!!  reason for doing so is to be able to apply a postulate (see below) linking the
!!  the error E(m) to a postulated error E(p) for a Runge Kutta step of order p >= m.
!!  We have, from the definition of the error goal for the Runge Kutta step:
!!
!!                          y(x+h) - y(x) = eFrac * eBase
!!
!!  eFrac is a single number that, when combined with the error base vector eBase, gives
!!  an advancement error limit for each 'y' component. Hence we can rewrite the above
!!  Taylor expansion as an expansion for eFrac:
!!
!!          eFrac = hy'(x)/eBase + h^2y''(x)/(2*eBase) + ... + h^my^m(x)/(m!*eBase)
!!
!!  and the eFrac truncation error is:
!!
!!           E(m) = h^(m+1) * y^(m+1)(x) / ((m+1)! * eBase)
!!
!!  We next state a posutale about the behavior of E(m) when estimating a Runge Kutta step
!!  error of order p >= m:
!!
!!           E(p) = E(m)^((p+1)/(m+1))
!!
!!  This assumes a power law decrease of the error as the Runge Kutta order p increases,
!!  provided the Taylor series eFrac error E(m) is < 1. If eFrac is < 1 and the Taylor series
!!  is continuously converging, then E(p) <= E(m).
!!
!!  As a first application of these ideas, we use m = 0 for the initial step size estimate.
!!  This means, that we use the first (linear) term in 'h' of the eFrac Taylor expansion series:
!!
!!           eFrac = [hy'(x)/eBase]^(p+1)
!!
!!  which gives for h:
!!
!!               h = eFrac^(1/(p+1)) * |eBase / y'(x)|
!!
!!  Note, that this last equation is actually a set of equations for each i-th component
!!  of both vectors eBase and y'(x):
!!
!!              h(i) = eFrac^(1/(p+1)) * |eBase (i) / [y'(x)](i)|
!!
!!  and the step size estimate must be such that it leads to acceptable errors for all
!!  individual components:
!!
!!                 h = min { h(i) }
!!
!!  One problem can arise, if a particular y'(x) component is equal to zero or very small,
!!  as this would lead to a very large h(i). If some of the other h(i) are well behaved, this
!!  situation would not be a problem due to taking the minimum, however if all y'(x) are very
!!  small, the situation leads to an overly large and unacceptable step size. For this reason
!!  an (optional) maximum allowed step size estimate can be provided by the user, safeguarding
!!  against any unreasonably large step size estimate.
!!
!! ARGUMENTS
!!
!!  method       : the type of RK method for which the estimate is
!!  f            : the ODE function containing details of the ODE system to be solved
!!  x            : the initial value of the independent variable
!!  y            : the initial values of the dependent variable(s)
!!  eFrac        : the fractional value for each of the error base values
!!  eBase        : the error base values for each of the dependent variables
!!  hmax         : the safeguard maximum allowed step size (optional)
!!
!! NOTES
!!
!!  Contains an interfaced function in the declarations (see below).
!!
!!***

real function RungeKutta_stepSizeEstimate (method,f,x,y,eFrac,eBase,hmax)

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
  real, optional,    intent (in)  :: hmax

  RungeKutta_stepSizeEstimate = 0.0

  return
end function RungeKutta_stepSizeEstimate
