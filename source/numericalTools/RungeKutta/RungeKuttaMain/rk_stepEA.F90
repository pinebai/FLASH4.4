!!****if* source/numericalTools/RungeKutta/RungeKuttaMain/rk_stepEA
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

  use RungeKutta_data,   ONLY : rk_a,                    &
                                rk_b,                    &
                                rk_c,                    &
                                rk_initialFunction,      &
                                rk_kVectors,             &
                                rk_stepSizeSafetyFactor

  use Driver_interface,  ONLY : Driver_abortFlash

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

  logical :: badStep
  logical :: stagnation

  integer :: i,j

  real    :: eMax
  real    :: h
!
!
!   ...Do the embedded adaptive Runge Kutta step.
!
!
  h = htry
  badStep = .true.

  rk_initialFunction (1:n) = f (x, y (1:n))      ! save initial function (no recompute, if h changes)

  stepSizeLoop: do while (badStep)

     rk_kVectors (1:n,1) = h * rk_initialFunction (1:n)

     do i = 2,s
        yout (1:n) = y (1:n)                                              ! temporarily use the 'yout'
        do j = 1,i-1                                                      ! vectors for accumulating
           yout (1:n) = yout (1:n) + rk_a (i,j) * rk_kVectors (1:n,j)     ! the k vectors for error
        end do                                                            ! analysis
        rk_kVectors (1:n,i) = h * f (x + rk_c (i) * h, yout (1:n))
     end do

     eout (1:n) = (rk_b (1,1) - rk_b (1,2)) * rk_kVectors (1:n,1)

     do i = 2,s
        eout (1:n) = eout (1:n) + (rk_b (i,1) - rk_b (i,2)) * rk_kVectors (1:n,i)
     end do

     hused   = h
     eMax    = maxval (abs (eout (1:n)) / eBase (1:n)) / eFrac           ! must 1st take abs of eout values
     badStep = (eMax > 1.0)

     if (badStep) then
         h = h * rk_stepSizeSafetyFactor * eMax ** ePower                ! decrease step size
     else
         yout (1:n) = y (1:n) + rk_b (1,1) * rk_kVectors (1:n,1)         ! add advancement vector

         do i = 2,s
            yout (1:n) = yout (1:n) + rk_b (i,1) * rk_kVectors (1:n,i)   ! add advancement vector
         end do

         h = h * rk_stepSizeSafetyFactor * max (eMin,eMax) ** ePower     ! increase step size

     end if

     hnext = h

  end do stepSizeLoop
!
!
!   ...Alert, if complete non-advancement of dependent variables. Note, that this check
!      has to be done AFTER you evaluated the new dependent variables. There is no
!      point in checking, separately, if the advancement vector for y (1:n) is equal
!      to zero, as the advancement vector can have tiny elements which get lost during
!      addition to y (1:n) due to finite working precision. 
!
!
  stagnation = (maxval (abs (y (1:n) - yout (1:n))) == 0.0)

  if (stagnation) then
      call Driver_abortFlash ('[rk_stepEA] ERROR: no advancement of dependent variables!')
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine rk_stepEA
