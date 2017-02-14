!!****if* source/numericalTools/Roots/RootsMain/Roots_x2Polynomial
!!
!! NAME
!!
!!  Roots_x2Polynomial
!!
!! SYNOPSIS
!!
!!  call Roots_x2Polynomial (real,    intent (in)  :: q1,
!!                           real,    intent (in)  :: q0,
!!                           integer, intent (out) :: nReal,
!!                           real,    intent (out) :: root (1:2,1:2))
!!
!! DESCRIPTION
!!
!!  Calculates all real + complex roots of the quadratic polynomial:
!!
!!                 x^2 + q1 * x + q0
!!
!!  The code checks internally, if rescaling of the coefficients is needed to
!!  avoid overflow.
!!
!!  The order of the roots is as follows:
!!
!!        1) For real roots, the order is according to their algebraic value
!!           on the number scale (largest positive first, largest negative last).
!!
!!        2) Since there can be only one complex conjugate pair root, no order
!!           is necessary.
!!
!! ARGUMENTS
!!
!!  q1         : coefficient of x term
!!  q0         : independent coefficient
!!  nReal      : number of real roots found
!!  root (n,1) : real part of n-th root
!!  root (n,2) : imaginary part of n-th root
!!
!! NOTES
!!
!!***

subroutine Roots_x2Polynomial (q1, q0,         &
                                        nReal, &
                                        root   )

  use Roots_data, ONLY: rt_LPN,     &
                        rt_sqrtLPN

  implicit none

  real   , intent (in)  :: q1, q0
  integer, intent (out) :: nReal
  real   , intent (out) :: root (1:2,1:2)

  logical :: rescale

  real :: a0, a1
  real :: k, x, y, z

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2
!
!
!     ...Handle special cases.
!
!
  if (q0 == 0.0 .and. q1 == 0.0) then

      nReal = 2

      root (:,Re) = 0.0
      root (:,Im) = 0.0

  else if (q0 == 0.0) then

      nReal = 2

      root (1,Re) = max (0.0, - q1)
      root (2,Re) = min (0.0, - q1)
      root (:,Im) = 0.0

  else if (q1 == 0.0) then

      x = sqrt (abs (q0))

      if (q0 < 0.0) then

          nReal = 2

          root (1,Re) = x
          root (2,Re) = - x
          root (:,Im) = 0.0

      else

          nReal = 0

          root (:,Re) = 0.0
          root (1,Im) = x
          root (2,Im) = - x

      end if

  else
!
!
!     ...The general case. Do rescaling, if either squaring of q1/2 or evaluation of
!        (q1/2)^2 - q0 will lead to overflow. This is better than to have the solver
!        crashed. Note, that rescaling might lead to loss of accuracy, so we only
!        invoke it when absolutely necessary.
!
!
      rescale = (q1 > rt_sqrtLPN + rt_sqrtLPN)    ! this detects overflow of (q1/2)^2

      if (.not.rescale) then
           x = q1 * 0.5                           ! we are sure here that x*x will not overflow
           rescale = (q0 < x * x - rt_LPN)        ! this detects overflow of (q1/2)^2 - q0
      end if

      if (rescale) then

          x = abs (q1)
          y = sqrt (abs (q0))

          if (x > y) then
              k  = x
              z  = 1.0 / x
              a1 = sign (1.0, q1)
              a0 = (q0 * z) * z
          else
              k  = y
              a1 = q1 / y
              a0 = sign (1.0, q0)
          end if

      else
          a1 = q1
          a0 = q0
      end if
!
!
!     ...Determine the roots of the quadratic. Note, that either a1 or a0 might
!        have become equal to zero due to underflow. But both cannot be zero.
!
!
      x = a1 * 0.5
      y = x * x - a0

      if (y >= 0.0) then

          y = sqrt (y)

          if (x > 0.0) then
              y = - x - y
          else
              y = - x + y
          end if

          if (rescale) then
              y = y * k                       ! very important to convert to original
              z = q0 / y                      ! root first, otherwise complete loss of
          else                                ! root due to possible a0 = 0 underflow
              z = a0 / y
          end if

          nReal = 2

          root (1,Re) = max (y,z)             ! 1st real root of x^2 + a1 * x + a0
          root (2,Re) = min (y,z)             ! 2nd real root of x^2 + a1 * x + a0
          root (:,Im) = 0.0

      else

          y = sqrt (- y)

          nReal = 0

          root (1,Re) = - x
          root (2,Re) = - x
          root (1,Im) = y                     ! complex conjugate pair of roots
          root (2,Im) = - y                   ! of x^2 + a1 * x + a0

          if (rescale) then
              root = root * k
          end if

      end if

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine Roots_x2Polynomial
