!!****if* source/numericalTools/RungeKutta/localAPI/rk_stepEA
!!
!! NAME
!!
!!  rk_stepEA
!!
!! SYNOPSIS
!!
!!  call rk_stepEA (function, intent (in)  :: f,
!!                  integer,  intent (in)  :: n,
!!                  integer,  intent (in)  :: s,
!!                  real,     intent (in)  :: x,
!!                  real,     intent (in)  :: y     (:),
!!                  real,     intent (in)  :: eMin,
!!                  real,     intent (in)  :: ePower,
!!                  real,     intent (in)  :: eFrac,
!!                  real,     intent (in)  :: eBase (:),
!!                  real,     intent (in)  :: htry,
!!                  real,     intent (out) :: hused,
!!                  real,     intent (out) :: hnext,
!!                  real,     intent (out) :: yout  (:),
!!                  real,     intent (out) :: eout  (:))
!!
!! DESCRIPTION
!!
!!  This routine performs a single (E)mbedded (A)daptive Runge Kuta (RK) step.
!!  The type of embedded adaptive step performed is determined by the previously copied
!!  Butcher tableau a,b,c values. The Butcher tableau copied contains implicitly the order,
!!  its size and the kind of adapted step to be performed.
!!
!!  The details of the ODE function are contained in the passed array function 'f', which
!!  constitutes the link of the general Runge Kutta module to the individual applications.
!!  Thus 'f' has to be provided by an application outside of the Runge Kutta unit.
!!
!! ARGUMENTS
!!
!!  f            : the ODE function containing details of the ODE system to be solved
!!  n            : the number of dependent variables
!!  s            : the Butcher tableau size
!!  x            : the initial value of the independent variable
!!  y            : the initial values of the dependent variable(s)
!!  eMin         : the lowest allowed relative error (limits increase of step size)
!!  ePower       : the error rescaling power (equal to -1/(p+1), p = order of lowest order RK step)
!!  eFrac        : the fractional value for each of the error base values
!!  eBase        : the error base values for each of the dependent variables
!!  htry         : the initial step size on the independent variable
!!  hused        : the final step size applied on the independent variable
!!  hnext        : the suggested next step size on the independent variable
!!  yout         : the values of the dependent variable(s) after the accepted RK step
!!  eout         : the returned error of all dependent variable(s)
!!
!! NOTES
!!
!!  1) Contains an interfaced function in the declarations (see below).
!!
!!  2) The routine should be threadsave: the intermediate arrays and pointers used
!!     (rk_a, rk_b, rk_c, rk_initialFunction, rk_kVectors) have all been declared
!!     as 'threadprivate'.
!!
!!  3) For performance reasons, no local arrays are declared here. This is to avoid
!!     repeated allocation/deallocation of memory in a routine which can possibly
!!     be called many times (millions!) during an application.
!!
!!***

subroutine rk_stepEA (f, n, s,                      &
                      x, y,                         &
                      eMin, ePower, eFrac, eBase,   &
                      htry,                         &
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

  integer, intent (in)  :: n
  integer, intent (in)  :: s
  real,    intent (in)  :: x
  real,    intent (in)  :: y     (:)
  real,    intent (in)  :: eMin
  real,    intent (in)  :: ePower
  real,    intent (in)  :: eFrac
  real,    intent (in)  :: eBase (:)
  real,    intent (in)  :: htry
  real,    intent (out) :: hused
  real,    intent (out) :: hnext
  real,    intent (out) :: yout  (:)
  real,    intent (out) :: eout  (:)

  hused = 0.0
  hnext = 0.0

  yout (:) = 0.0
  eout (:) = 0.0

  return
end subroutine rk_stepEA
