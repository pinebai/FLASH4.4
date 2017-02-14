!!****if* source/numericalTools/RungeKutta/RungeKuttaMain/rk_stepEAC
!!
!! NAME
!!
!!  rk_stepEAC
!!
!! SYNOPSIS
!!
!!  call rk_stepEAC (function, intent (in)  :: f,
!!                   integer,  intent (in)  :: n,
!!                   integer,  intent (in)  :: nc,
!!                   integer,  intent (in)  :: s,
!!                   real,     intent (in)  :: x,
!!                   real,     intent (in)  :: y     (:),
!!                   function, intent (in)  :: ymin,
!!                   function, intent (in)  :: ymax,
!!                   real,     intent (in)  :: eMin,
!!                   real,     intent (in)  :: ePower,
!!                   real,     intent (in)  :: eFrac,
!!                   real,     intent (in)  :: eBase (:),
!!                   real,     intent (in)  :: htry,
!!                   real,     intent (out) :: hused,
!!                   real,     intent (out) :: hnext,
!!                   real,     intent (out) :: yout  (:),
!!                   real,     intent (out) :: eout  (:))
!!
!! DESCRIPTION
!!
!!  This routine performs a single (E)mbedded (A)daptive (C)onfined Runge Kuta (RK)
!!  step. The type of embedded adaptive step performed is determined by the previously copied
!!  Butcher tableau a,b,c values. The Butcher tableau copied contains implicitly the order,
!!  its size and the kind of adapted step to be performed.
!!
!!  The details of the ODE function are contained in the passed array function 'f', which
!!  constitutes the link of the general Runge Kutta module to the individual applications.
!!  Thus 'f' has to be provided by an application outside of the Runge Kutta unit.
!!
!!  The code performs a confined RK step. The confinement boundaries are (potentially) dependent
!!  on the current value of the dependent variables at the intermediate and final RK points. This
!!  boundary dependency has to be transmitted via via 2 functions 'ymin' and 'ymax', which accept
!!  the current dependent variables array as argument input. The user has to provide the external
!!  definition of these 2 functions.
!!
!!  The code checks if at any stage during the RK step the intermediate RK points (i.e. dependent
!!  variables) are out of bounds for the confined dependent variables. If the case, the step
!!  size must be reduced (currently by a fixed factor < 1, but this factor might be passed in
!!  the future as argument as well).
!!
!!  Note that the layout of the dependent variable array can always be done such that the
!!  confined variables come first. The user must therefore be careful to design its ODE
!!  function in such a way that those confined dependent variables are at the beginning
!!  of the dependent variable array.
!!
!! ARGUMENTS
!!
!!  f            : the ODE function containing details of the ODE system to be solved
!!  n            : the number of dependent variables
!!  nc           : the number of confined dependent variables (nc <= n, not checked here!)
!!  s            : the Butcher tableau size
!!  x            : the initial value of the independent variable
!!  y            : the initial values of the dependent variable(s)
!!  ymin         : the lower bound confinement function (as function of dependent variables)
!!  ymax         : the upper bound confinement function (as function of dependent variables)
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
!!  1) Contains interfaced functions in the declarations (see below).
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

subroutine rk_stepEAC (f, n, nc, s,                  &
                       x, y, ymin, ymax,             &
                       eMin, ePower, eFrac, eBase,   &
                       htry,                         &
                                       hused, hnext, &
                                       yout, eout    )

  use RungeKutta_data,   ONLY : rk_a,                         &
                                rk_b,                         &
                                rk_c,                         &
                                rk_initialFunction,           &
                                rk_kVectors,                  &
                                rk_stepSizeConfinementFactor, &
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

  integer, intent (in)  :: n
  integer, intent (in)  :: nc
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
!   ...Do the embedded adaptive confined Runge Kutta step.
!
!
  h = htry
  badStep = .true.

  rk_initialFunction (1:n) = f (x, y (1:n))     ! save initial function (no recompute, if h changes)

  stepSizeLoop: do while (badStep)

     rk_kVectors (1:n,1) = h * rk_initialFunction (1:n)

     do i = 2,s
        yout (1:n) = y (1:n)                                              ! temporarily use the 'yout'
        do j = 1,i-1                                                      ! vectors for accumulating
           yout (1:n) = yout (1:n) + rk_a (i,j) * rk_kVectors (1:n,j)     ! the k vectors for error
        end do                                                            ! analysis

        if (any (yout (1:nc) < ymin (nc, yout (1:n)) .or. &               ! check, if intermediate
                 yout (1:nc) > ymax (nc, yout (1:n)))     ) then          ! RK points step outside
            h = h * rk_stepSizeConfinementFactor                          ! the confinement barriers
            cycle stepSizeLoop
        end if

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

         if (any (yout (1:nc) < ymin (nc, yout (1:n)) .or. &             ! check, if RK endpoints
                  yout (1:nc) > ymax (nc, yout (1:n)))     ) then        ! step outside the confinement
             h = h * rk_stepSizeConfinementFactor                        ! barriers
             badStep = .true.
             cycle stepSizeLoop
         end if

         h = h * rk_stepSizeSafetyFactor * max (eMin,eMax) ** ePower     ! increase step size

     end if

     hnext = h

  end do stepSizeLoop
!
!
!   ...Alert, if complete non-advancement of dependent variables. Note, that this check
!      has to be done AFTER we evaluated the new dependent variables. There is no
!      point in checking, separately, if the advancement vector for y (1:n) is equal
!      to zero, as the advancement vector can have tiny elements which get lost during
!      addition to y (1:n) due to finite working precision. 
!
!
  stagnation = (maxval (abs (y (1:n) - yout (1:n))) == 0.0)

  if (stagnation) then
      call Driver_abortFlash ('[rk_stepEAC] ERROR: no advancement of dependent variables!')
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine rk_stepEAC
