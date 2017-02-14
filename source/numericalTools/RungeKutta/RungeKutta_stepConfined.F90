!!****f* source/numericalTools/RungeKutta/RungeKutta_stepConfined
!!
!! NAME
!!
!!  RungeKutta_stepConfined
!!
!! SYNOPSIS
!!
!!  call RungeKutta_stepConfined (character (len=*), intent (in)  :: method,
!!                                function,          intent (in)  :: f,
!!                                integer,           intent (in)  :: nc,
!!                                real,              intent (in)  :: x,
!!                                real,              intent (in)  :: y     (:),
!!                                function,          intent (in)  :: ymin,
!!                                function,          intent (in)  :: ymax,
!!                                real,              intent (in)  :: eFrac,
!!                                real,              intent (in)  :: eBase (:),
!!                                real,              intent (in)  :: htry,
!!                                real,              intent (out) :: hused,
!!                                real,              intent (out) :: hnext,
!!                                real,              intent (out) :: yout  (:),
!!                                real,              intent (out) :: eout  (:))
!!
!! DESCRIPTION
!!
!!  This routine performs a single confined Runge Kutta step. It copies the appropriate
!!  Butcher tableau data from the repository and calls the corresponding RK stepper. The details
!!  of the ODE function are contained in the passed array function 'f', which constitutes
!!  the link of the general Runge Kutta module to the individual applications. Thus 'f' has
!!  to be provided by an application outside of the Runge Kutta unit.
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
!!  The code performs a confined RK step. The confinement boundaries are (potentially)
!!  dependent on the current value of the dependent variables at the intermediate and final
!!  RK points. This boundary dependency has to be transmitted via via 2 functions 'ymin' and
!!  'ymax', which accept the current dependent variables array as argument input. The user
!!  has to provide the external definition of these 2 functions.
!!
!!  Note that the layout of the dependent variable array can always be done such that the
!!  confined variables come first. The user must therefore be careful to design its ODE
!!  function in such a way that those confined dependent variables are at the beginning
!!  of the dependent variable array.
!!
!! ARGUMENTS
!!
!!  method       : the type of RK method to be applied
!!  f            : the ODE function containing details of the ODE system to be solved
!!  nc           : the number of confined dependent variables (must be <= # dep variables)
!!  x            : the initial value of the independent variable
!!  y            : the initial values of the dependent variable(s)
!!  ymin         : the lower bound confinement function (as function of dependent variables)
!!  ymax         : the upper bound confinement function (as function of dependent variables)
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
!!  1) Contains an interfaced function in the declarations (see below).
!!
!!  2) The sizes of the arrays eBase, yout and eout can be larger than needed, i.e.,
!!     they must not match exactly the number of dependent variables. Likewise with
!!     the sizes of the boundary arrays ymin and ymax, which can be greater than
!!     the number of confined dependent variables.
!!
!!***

subroutine RungeKutta_stepConfined (method,                      &
                                    f, nc,                       &
                                    x, y, ymin, ymax,            &
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

  interface
    function ymin (nc,y)
      integer, intent (in) :: nc
      real,    intent (in) :: y (:)
      real                 :: ymin (1:nc)
    end function ymin
  end interface

  interface
    function ymax (nc,y)
      integer, intent (in) :: nc
      real,    intent (in) :: y (:)
      real                 :: ymax (1:nc)
    end function ymax
  end interface

  character (len=*), intent (in)  :: method
  integer,           intent (in)  :: nc
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
end subroutine RungeKutta_stepConfined
