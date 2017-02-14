!!****f* source/numericalTools/RungeKutta/RungeKutta_step
!!
!! NAME
!!
!!  RungeKutta_step
!!
!! SYNOPSIS
!!
!!  call RungeKutta_step (character (len=*), intent (in)  :: method,
!!                        function,          intent (in)  :: f,
!!                        real,              intent (in)  :: x,
!!                        real,              intent (in)  :: y     (:),
!!                        real,              intent (in)  :: eFrac,
!!                        real,              intent (in)  :: eBase (:),
!!                        real,              intent (in)  :: htry,
!!                        real,              intent (out) :: hused,
!!                        real,              intent (out) :: hnext,
!!                        real,              intent (out) :: yout  (:),
!!                        real,              intent (out) :: eout  (:))
!!
!! DESCRIPTION
!!
!!  This routine performs a single adaptive Runge Kutta step. It copies the appropriate
!!  Butcher tableau data from the repository and calls the corresponding RK stepper.
!!  The details of the ODE function are contained in the passed array function 'f',
!!  which constitutes the link of the general Runge Kutta module to the individual
!!  applications. Thus 'f' has to be provided by an application outside of the Runge
!!  Kutta unit.
!!
!!  The error vector for the dependent variables is composed of two pieces (as suggested
!!  by Numerical Recipies):
!!
!!                            e (i) = eFrac * eBase (i)
!!
!!  In this error vector definition, eFrac denotes the fractional value each of the
!!  individual eBase (i) values that define the total allowed maximum error for each
!!  of the i-th dependent variables.
!!
!! ARGUMENTS
!!
!!  method       : the type of RK method to be applied
!!  f            : the ODE function containing details of the ODE system to be solved
!!  x            : the initial value of the independent variable
!!  y            : the initial values of the dependent variable(s)
!!  eFrac        : the fractional value for each of the error base values
!!  eBase        : the error base values for each of the dependent variables
!!  htry         : the initial step size on the independent variable
!!  hused        : the final step size applied on the independent variable
!!  hnext        : the suggested next step size on the independent variable
!!  yout         : the values of the dependent variable(s) after the accepted RK45 step
!!  eout         : the returned error of all dependent variable(s)
!!
!! NOTES
!!
!!  Contains an interfaced function in the declarations (see below).
!!
!!***

subroutine RungeKutta_step (method,                      &
                            f, x, y,                     &
                            eFrac, eBase,                &
                            htry,                        &
                                           hused, hnext, &
                                           yout, eout    )

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
  real,              intent (in)  :: htry
  real,              intent (out) :: hused
  real,              intent (out) :: hnext
  real,              intent (out) :: yout  (:)
  real,              intent (out) :: eout  (:)

  hused = 0.0
  hnext = 0.0

  yout (:) = 0.0
  eout (:) = 0.0

  return
end subroutine RungeKutta_step
