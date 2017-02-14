!!****if* source/numericalTools/Roots/RootsMain/Roots_x3Polynomial
!!
!! NAME
!!
!!  Roots_x3Polynomial
!!
!! SYNOPSIS
!!
!!  call Roots_x3Polynomial (real,              intent (in)  :: c2,
!!                           real,              intent (in)  :: c1,
!!                           real,              intent (in)  :: c0,
!!                           integer,           intent (out) :: nReal,
!!                           real,              intent (out) :: root (1:3,1:2),
!!                           logical, optional, intent (in)  :: printInfo,
!!                           integer, optional, intent (in)  :: printUnit)
!!
!! DESCRIPTION
!!
!!  Calculates all real + complex roots of the cubic polynomial:
!!
!!                 x^3 + c2 * x^2 + c1 * x + c0
!!
!!  The first real root (which always exists) is obtained using an optimized
!!  Newton-Raphson scheme. The other remaining roots are obtained through
!!  composite deflation to a quadratic.
!!
!!  The cubic root solver can handle any size of cubic coefficients and there is
!!  no danger of overflow due to proper rescaling of the cubic polynomial.
!!
!!  The order of the roots is as follows:
!!
!!        1) For real roots, the order is according to their algebraic value
!!           on the number scale (largest positive first, largest negative last).
!!
!!        2) Since there can be only one complex conjugate pair root, no order
!!           is necessary.
!!
!!        3) All real roots preceede the complex ones.
!!
!! ARGUMENTS
!!
!!  c2         : coefficient of x^2 term
!!  c1         : coefficient of x term
!!  c0         : independent coefficient
!!  nReal      : number of different real roots found
!!  root (n,1) : real part of n-th root
!!  root (n,2) : imaginary part of n-th root
!!  printInfo  : if given and true, detailed info will be printed about intermediate stages
!!  printUnit  : the unit ID, where the info will be printed
!!
!! NOTES
!!
!!  Only passing both printInfo AND printUnit will result in printing out the info.
!!  Giving only one of them results in no printing action.
!!
!!***

subroutine Roots_x3Polynomial (c2, c1, c0,         &
                                            nReal, &
                                            root,  &
                               printInfo,          &
                               printUnit           )

  use Roots_data,      ONLY: rt_macheps
  use Roots_interface, ONLY: Roots_x2Polynomial

  implicit none

  real   ,           intent (in)  :: c2, c1, c0
  integer,           intent (out) :: nReal
  real   ,           intent (out) :: root (1:3,1:2)
  logical, optional, intent (in)  :: printInfo
  integer, optional, intent (in)  :: printUnit

  character (len=*), parameter :: printFormat = "(a,es25.16)"

  logical :: bisection
  logical :: converged
  logical :: doPrint = .false.            ! the default

  integer :: cubicType
  integer :: deflateCase
  integer :: oscillate

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2

  real :: a0, a1, a2
  real :: a, b, c, k, s, t, u, x, y, z
  real :: xShift

  real, parameter :: one27th = 1.0 / 27.0
  real, parameter :: two27th = 2.0 / 27.0
  real, parameter :: third   = 1.0 /  3.0

  real, parameter :: p1 = 1.09574         !
  real, parameter :: q1 = 3.23900e-1      ! Newton-Raphson coeffs for class 1 and 2
  real, parameter :: r1 = 3.23900e-1      !
  real, parameter :: s1 = 9.57439e-2      !

  real, parameter :: p3 = 1.14413         !
  real, parameter :: q3 = 2.75509e-1      ! Newton-Raphson coeffs for class 3
  real, parameter :: r3 = 4.45578e-1      !
  real, parameter :: s3 = 2.59342e-2      !

  real, parameter :: q4 = 7.71845e-1      ! Newton-Raphson coeffs for class 4
  real, parameter :: s4 = 2.28155e-1      !

  real, parameter :: p51 = 8.78558e-1     !
  real, parameter :: p52 = 1.92823e-1     !
  real, parameter :: p53 = 1.19748        !
  real, parameter :: p54 = 3.45219e-1     !
  real, parameter :: q51 = 5.71888e-1     !
  real, parameter :: q52 = 5.66324e-1     !
  real, parameter :: q53 = 2.83772e-1     ! Newton-Raphson coeffs for class 5 and 6
  real, parameter :: q54 = 4.01231e-1     !
  real, parameter :: r51 = 7.11154e-1     !
  real, parameter :: r52 = 5.05734e-1     !
  real, parameter :: r53 = 8.37476e-1     !
  real, parameter :: r54 = 2.07216e-1     !
  real, parameter :: s51 = 3.22313e-1     !
  real, parameter :: s52 = 2.64881e-1     !
  real, parameter :: s53 = 3.56228e-1     !
  real, parameter :: s54 = 4.45532e-3     !
!
!
!     ...Start.
!
!
  if (present (printInfo)) then
      doPrint = printInfo .and. present (printUnit)
  end if

  if (doPrint) then
      write (printUnit,printFormat) ' initial cubic c2      = ',c2
      write (printUnit,printFormat) ' initial cubic c1      = ',c1
      write (printUnit,printFormat) ' initial cubic c0      = ',c0
      write (printUnit,printFormat) ' ------------------------------------------------'
  end if
!
!
!     ...Handle special cases.
!
!            1) all terms zero
!            2) only quadratic term is nonzero -> linear equation.
!            3) only independent term is zero -> quadratic equation.
!
!
  if (c0 == 0.0 .and. c1 == 0.0 .and. c2 == 0.0) then

      cubicType = 0

  else if (c0 == 0.0 .and. c1 == 0.0) then

      k  = 1.0
      a2 = c2

      cubicType = 1

  else if (c0 == 0.0) then

      k  = 1.0
      a2 = c2
      a1 = c1

      cubicType = 2

  else
!
!
!     ...The general case. Rescale cubic polynomial, such that largest absolute coefficient
!        is (exactly!) equal to 1. Honor the presence of a special cubic case that might have
!        been obtained (due to underflow in the coefficients).
!
!
      x = abs (c2)
      y = sqrt (abs (c1))
      z = abs (c0) ** third
      u = max (x,y,z)

      if (u == x) then

          k  = 1.0 / x
          a2 = sign (1.0, c2)
          a1 = (c1 * k) * k
          a0 = ((c0 * k) * k) * k

      else if (u == y) then

          k  = 1.0 / y
          a2 = c2 * k
          a1 = sign (1.0, c1)
          a0 = ((c0 * k) * k) * k

      else

          k  = 1.0 / z
          a2 = c2 * k
          a1 = (c1 * k) * k
          a0 = sign (1.0, c0)

      end if

      if (doPrint) then
          write (printUnit,printFormat) ' rescaling factor      = ',k
          write (printUnit,printFormat) ' ------------------------------------------------'
          write (printUnit,printFormat) ' rescaled cubic c2     = ',a2
          write (printUnit,printFormat) ' rescaled cubic c1     = ',a1
          write (printUnit,printFormat) ' rescaled cubic c0     = ',a0
          write (printUnit,printFormat) ' ------------------------------------------------'
      end if

      k = 1.0 / k

      if (a0 == 0.0 .and. a1 == 0.0 .and. a2 == 0.0) then
          cubicType = 0
      else if (a0 == 0.0 .and. a1 == 0.0) then
          cubicType = 1
      else if (a0 == 0.0) then
          cubicType = 2
      else
          cubicType = 3
      end if

  end if
!
!
!     ...Select the case.
!
!        1) Only zero roots.
!
!
  select case (cubicType)

    case (0)

      nReal = 3

      root (:,Re) = 0.0
      root (:,Im) = 0.0
!
!
!     ...2) The linear equation case -> additional 2 zeros.
!
!
    case (1)

      x = - a2 * k

      nReal = 3

      root (1,Re) = max (0.0, x)
      root (2,Re) = 0.0
      root (3,Re) = min (0.0, x)
      root (:,Im) = 0.0
!
!
!     ...3) The quadratic equation case -> additional 1 zero.
!
!
    case (2)

      call Roots_x2Polynomial (a2, a1, nReal, root (1:2,1:2))

      if (nReal == 2) then

          x = root (1,1) * k         ! real roots of quadratic are ordered x >= y
          y = root (2,1) * k

          nReal = 3

          root (1,Re) = max (x, 0.0)
          root (2,Re) = max (y, min (x, 0.0))
          root (3,Re) = min (y, 0.0)
          root (:,Im) = 0.0

      else

          nReal = 1

          root (3,Re) = root (2,Re) * k
          root (2,Re) = root (1,Re) * k
          root (1,Re) = 0.0
          root (3,Im) = root (2,Im) * k
          root (2,Im) = root (1,Im) * k
          root (1,Im) = 0.0

      end if
!
!
!     ...3) The general cubic case. Set the best Newton-Raphson root estimates for the cubic.
!           The easiest and most robust conditions are checked first. The most complicated
!           ones are last and only done when absolutely necessary.
!
!
    case (3)

      if (a0 == 1.0) then

          x = - p1 + q1 * a1 - a2 * (r1 - s1 * a1)

          a = a2
          b = a1
          c = a0
          xShift = 0.0

      else if (a0 == - 1.0) then

          x = p1 - q1 * a1 - a2 * (r1 - s1 * a1)

          a = a2
          b = a1
          c = a0
          xShift = 0.0

      else if (a1 == 1.0) then

          if (a0 > 0.0) then
              x = a0 * (- q4 - s4 * a2)
          else
              x = a0 * (- q4 + s4 * a2)
          end if

          a = a2
          b = a1
          c = a0
          xShift = 0.0

      else if (a1 == - 1.0) then

          y = - two27th
          y = y * a2
          y = y * a2 - third
          y = y * a2

          if (a0 < y) then
              x = + p3 - q3 * a0 - a2 * (r3 + s3 * a0)               ! + guess
          else
              x = - p3 - q3 * a0 - a2 * (r3 - s3 * a0)               ! - guess
          end if

          a = a2
          b = a1
          c = a0
          xShift = 0.0

      else if (a2 == 1.0) then

          b = a1 - third
          c = a0 - one27th

          if (abs (b) < rt_macheps .and. abs (c) < rt_macheps) then  ! triple -1/3 root

              x = - third * k

              nReal = 3

              root (:,Re) = x
              root (:,Im) = 0.0

              return

          else

              y = third * a1 - two27th

              if (a1 <= third) then
                  if (a0 > y) then
                      x = - p51 - q51 * a0 + a1 * (r51 - s51 * a0)   ! - guess
                  else
                      x = + p52 - q52 * a0 - a1 * (r52 + s52 * a0)   ! + guess
                  end if
              else
                  if (a0 > y) then
                      x = - p53 - q53 * a0 + a1 * (r53 - s53 * a0)   ! <-1/3 guess
                  else
                      x = + p54 - q54 * a0 - a1 * (r54 + s54 * a0)   ! >-1/3 guess
                  end if
              end if

              if (abs (b) < 1.e-2 .and. abs (c) < 1.e-2) then        ! use shifted root
                  c = - third * b + c
                  if (abs (c) < rt_macheps) c = 0.0                  ! prevent random noise
                  a = 0.0
                  xShift = third
                  x = x + xShift
              else
                  a = a2
                  b = a1
                  c = a0
                  xShift = 0.0
              end if

          end if

      else if (a2 == - 1.0) then

          b = a1 - third
          c = a0 + one27th

          if (abs (b) < rt_macheps .and. abs (c) < rt_macheps) then  ! triple 1/3 root

              x = third * k

              nReal = 3

              root (:,Re) = x
              root (:,Im) = 0.0

              return

          else

              y = two27th - third * a1

              if (a1 <= third) then
                  if (a0 < y) then
                      x = + p51 - q51 * a0 - a1 * (r51 + s51 * a0)   ! +1 guess
                  else
                      x = - p52 - q52 * a0 + a1 * (r52 - s52 * a0)   ! -1 guess
                  end if
              else
                  if (a0 < y) then
                      x = + p53 - q53 * a0 - a1 * (r53 + s53 * a0)   ! >1/3 guess
                  else
                      x = - p54 - q54 * a0 + a1 * (r54 - s54 * a0)   ! <1/3 guess
                  end if
              end if

              if (abs (b) < 1.e-2 .and. abs (c) < 1.e-2) then        ! use shifted root
                  c = third * b + c
                  if (abs (c) < rt_macheps) c = 0.0                  ! prevent random noise
                  a = 0.0
                  xShift = - third
                  x = x + xShift
              else
                  a = a2
                  b = a1
                  c = a0
                  xShift = 0.0
              end if

          end if

      end if
!
!
!     ...Perform Newton/Bisection iterations on x^3 + ax^2 + bx + c.
!
!
      z = x + a
      y = x + z
      z = z * x + b
      y = y * x + z       ! C'(x)
      z = z * x + c       ! C(x)
      t = z               ! save C(x) for sign comparison
      x = x - z / y       ! 1st improved root

      oscillate = 0
      bisection = .false.
      converged = .false.

      do while (.not.converged .and. .not.bisection)    ! Newton-Raphson iterates)

         z = x + a
         y = x + z
         z = z * x + b
         y = y * x + z
         z = z * x + c

         if (z * t < 0.0) then                          ! does Newton start oscillating ?
             if (z < 0.0) then
                 oscillate = oscillate + 1              ! increment oscillation counter
                 s = x                                  ! save lower bisection bound
             else
                 u = x                                  ! save upper bisection bound
             end if
             t = z                                      ! save current C(x)
         end if

         y = z / y                                      ! Newton correction
         x = x - y                                      ! new Newton root

         bisection = oscillate > 2                      ! activate bisection
         converged = abs (y) <= abs (x) * rt_macheps    ! Newton convergence indicator

         if (doPrint) write (printUnit,printFormat)     ' Newton root           = ',x

      end do

      if (bisection) then

          t = u - s                                     ! initial bisection interval
          do while (abs (t) > abs (x) * rt_macheps)     ! bisection iterates

             z = x + a                                  !
             z = z * x + b                              ! C (x)
             z = z * x + c                              !

             if (z < 0.0) then                          !
                 s = x                                  !
             else                                       ! keep bracket on root
                 u = x                                  !
             end if                                     !

             t = 0.5 * (u - s)                          ! new bisection interval
             x = s + t                                  ! new bisection root

             if (doPrint) write (printUnit,printFormat) ' Bisection root        = ',x

          end do
      end if

      if (doPrint) write (printUnit,printFormat) ' ------------------------------------------------'

      x = x - xShift                                    ! unshift root
!
!
!     ...Forward / backward deflate rescaled cubic (if needed) to check for other real roots.
!        The deflation analysis is performed on the rescaled cubic. The actual deflation must
!        be performed on the original cubic, not the rescaled one. Otherwise deflation errors
!        will be enhanced when undoing the rescaling on the extra roots.
!
!
      z = abs (x)
      s = abs (a2)
      t = abs (a1)
      u = abs (a0)

      y = z * max (s,z)           ! take maximum between |x^2|,|a2 * x|

      deflateCase = 1             ! up to now, the maximum is |x^3| or |a2 * x^2|

      if (y < t) then             ! check maximum between |x^2|,|a2 * x|,|a1|
          y = t * z               ! the maximum is |a1 * x|
          deflateCase = 2         ! up to now, the maximum is |a1 * x|
      else
          y = y * z               ! the maximum is |x^3| or |a2 * x^2|
      end if

      if (y < u) then             ! check maximum between |x^3|,|a2 * x^2|,|a1 * x|,|a0|
          deflateCase = 3         ! the maximum is |a0|
      end if

      y = x * k                   ! real root of original cubic

      select case (deflateCase)

      case (1)
        x = 1.0 / y
        t = - c0 * x              ! t -> backward deflation on unscaled cubic
        s = (t - c1) * x          ! s -> backward deflation on unscaled cubic
      case (2)
        s = c2 + y                ! s ->  forward deflation on unscaled cubic
        t = - c0 / y              ! t -> backward deflation on unscaled cubic
      case (3)
        s = c2 + y                ! s ->  forward deflation on unscaled cubic
        t = c1 + s * y            ! t ->  forward deflation on unscaled cubic
      end select

      if (doPrint) then
          write (printUnit,printFormat) ' Residual quadratic q1 = ',s
          write (printUnit,printFormat) ' Residual quadratic q0 = ',t
          write (printUnit,printFormat) ' ------------------------------------------------'
      end if

      call Roots_x2Polynomial (s, t, nReal, root (1:2,1:2))

      if (nReal == 2) then

          x = root (1,Re)         ! real roots of quadratic are ordered x >= z
          z = root (2,Re)         ! use 'z', because 'y' is original cubic real root

          nReal = 3

          root (1,Re) = max (x, y)
          root (2,Re) = max (z, min (x, y))
          root (3,Re) = min (z, y)
          root (:,Im) = 0.0

      else

          nReal = 1

          root (3,Re) = root (2,Re)
          root (2,Re) = root (1,Re)
          root (1,Re) = y
          root (3,Im) = root (2,Im)
          root (2,Im) = root (1,Im)
          root (1,Im) = 0.0

      end if

  end select
!
!
!     ...Ready!
!
!
  return
end subroutine Roots_x3Polynomial
