!!****if* source/numericalTools/RungeKutta/RungeKuttaMain/RungeKutta_stepConfined
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

  use RungeKutta_data,   ONLY : rk_a,                           &
                                rk_b,                           &
                                rk_c,                           &
                                rk_aTableauCashKarp45,          &
                                rk_bTableauCashKarp45,          &
                                rk_cTableauCashKarp45,          &
                                rk_aTableauBogShamp23,          &
                                rk_bTableauBogShamp23,          &
                                rk_cTableauBogShamp23,          &
                                rk_aTableauFehlberg34,          &
                                rk_bTableauFehlberg34,          &
                                rk_cTableauFehlberg34,          &
                                rk_aTableauFehlberg45,          &
                                rk_bTableauFehlberg45,          &
                                rk_cTableauFehlberg45,          &
                                rk_aTableauEulerHeu12,          &
                                rk_bTableauEulerHeu12,          &
                                rk_cTableauEulerHeu12,          &
                                rk_maxNumberDependentVariables 

  use rk_interface,      ONLY : rk_stepEAC

  use Driver_interface,  ONLY : Driver_abortFlash

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

  integer :: m,n
  real    :: eMin, ePower
!
!
!   ...Check:
!
!         1) if the stepper can handle the transmitted number of dependent variables.
!         2) if the arrays for the new dependent variables and their errors are of
!            appropriate size.
!         3) if the number of confined dependent variables is less or equal the total
!            number of dependent variables.
!
!
  m = min (size (eBase), size (yout), size (eout), rk_maxNumberDependentVariables)
  n = size (y)

  if (n > m) then
      call Driver_abortFlash ('[RungeKutta_stepConfined] ERROR: too many dependent variables!')
  end if

  if (nc > n) then
      call Driver_abortFlash ('[RungeKutta_stepConfined] ERROR: too many confined variables!')
  end if
!
!
!   ...Select the appropriate Runge Kutta method. eMin and ePower are defined as
!      follows (where p is the order of the lowest order RK step):
!
!         eMin = 1/10^(p+1), ensures that step size increase is never larger than x10.
!       ePower = - 1/(p+1), the exponent that governs decrease/increase of step size.
!
!
  select case (method)

  case ('EulerHeu12')

    rk_c  => rk_cTableauEulerHeu12
    rk_b  => rk_bTableauEulerHeu12
    rk_a  => rk_aTableauEulerHeu12
    eMin   =  1.e-2  ! p = 1
    ePower = -0.5    !

    call rk_stepEAC (f,n,nc,2,x,y,ymin,ymax,eMin,ePower,eFrac,eBase,htry,   hused,hnext,yout,eout)

  case ('BogShamp23')

    rk_c  => rk_cTableauBogShamp23
    rk_b  => rk_bTableauBogShamp23
    rk_a  => rk_aTableauBogShamp23
    eMin   =  1.e-3   ! p = 2
    ePower = -(1./3.) !

    call rk_stepEAC (f,n,nc,4,x,y,ymin,ymax,eMin,ePower,eFrac,eBase,htry,   hused,hnext,yout,eout)

  case ('Fehlberg34')

    rk_c  => rk_cTableauFehlberg34
    rk_b  => rk_bTableauFehlberg34
    rk_a  => rk_aTableauFehlberg34
    eMin   =  1.e-4  ! p = 3
    ePower = -0.25   !

    call rk_stepEAC (f,n,nc,5,x,y,ymin,ymax,eMin,ePower,eFrac,eBase,htry,   hused,hnext,yout,eout)

  case ('Fehlberg45')

    rk_c  => rk_cTableauFehlberg45
    rk_b  => rk_bTableauFehlberg45
    rk_a  => rk_aTableauFehlberg45
    eMin   =  1.e-5  ! p = 4
    ePower = -0.2    !

    call rk_stepEAC (f,n,nc,6,x,y,ymin,ymax,eMin,ePower,eFrac,eBase,htry,   hused,hnext,yout,eout)

  case ('CashKarp45')

    rk_c  => rk_cTableauCashKarp45
    rk_b  => rk_bTableauCashKarp45
    rk_a  => rk_aTableauCashKarp45
    eMin   =  1.e-5  ! p = 4
    ePower = -0.2    !

    call rk_stepEAC (f,n,nc,6,x,y,ymin,ymax,eMin,ePower,eFrac,eBase,htry,   hused,hnext,yout,eout)

  case default
        call Driver_abortFlash ('[RungeKutta_stepConfined] ERROR: unknown RK method')
  end select
!
!
!   ...Ready! 
!
!
  return
end subroutine RungeKutta_stepConfined
