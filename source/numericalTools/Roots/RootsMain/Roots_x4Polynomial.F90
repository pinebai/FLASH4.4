!!****if* source/numericalTools/Roots/RootsMain/Roots_x4Polynomial
!!
!! NAME
!!
!!  Roots_x4Polynomial
!!
!! SYNOPSIS
!!
!!  call Roots_x4Polynomial (real,              intent (in)  :: q3,
!!                           real,              intent (in)  :: q2,
!!                           real,              intent (in)  :: q1,
!!                           real,              intent (in)  :: q0,
!!                           integer,           intent (out) :: nReal,
!!                           real,              intent (out) :: root (1:4,1:2),
!!                           logical, optional, intent (in)  :: printInfo,
!!                           integer, optional, intent (in)  :: printUnit)
!!
!! DESCRIPTION
!!
!!  Calculates all real + complex roots of the quartic polynomial:
!!
!!                 x^4 + q3 * x^3 + q2 * x^2 + q1 * x + q0
!!
!!  An option for printing a detailed info about the intermediate stages in solving
!!  the quartic is available. Since the code has not yet been extensively tested,
!!  this enables a detailed check in case something went wrong and the roots obtained
!!  are not proper.
!!
!!  The quartic root solver can handle any size of quartic coefficients and there is
!!  no danger of overflow, due to proper rescaling of the quartic polynomial.
!!
!!  The order of the roots is as follows:
!!
!!        1) For real roots, the order is according to their algebraic value
!!           on the number scale (largest positive first, largest negative last).
!!
!!        2) For complex conjugate pair roots, the order is according to the
!!           algebraic value of their real parts (largest positive first). If
!!           the real parts are equal, the order is according to the algebraic
!!           value of their imaginary parts (largest first).
!!
!!        3) All real roots preceede the complex ones.
!!
!! ARGUMENTS
!!
!!  q3         : coefficient of x^3 term
!!  q2         : coefficient of x^2 term
!!  q1         : coefficient of x term
!!  q0         : independent coefficient
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

subroutine Roots_x4Polynomial (q3, q2, q1, q0,        &
                                               nReal, &
                                               root,  &
                               printInfo,             &
                               printUnit              )
  
  use Roots_data,      ONLY: rt_macheps

  use Roots_interface, ONLY: Roots_x2Polynomial, &
                             Roots_x3Polynomial

  implicit none

  real   ,           intent (in)  :: q3, q2, q1, q0
  integer,           intent (out) :: nReal
  real   ,           intent (out) :: root (1:4,1:2)
  logical, optional, intent (in)  :: printInfo
  integer, optional, intent (in)  :: printUnit

  character (len=*), parameter :: printFormat = "(a,es25.16)"

  logical :: bisection
  logical :: converged
  logical :: iterate
  logical :: minimum
  logical :: notZero
  logical :: doPrint = .false.     ! the default

  integer :: deflateCase
  integer :: oscillate
  integer :: quarticType

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2

  real :: a0, a1, a2, a3
  real :: a, b, c, d, k, s, t, u, x, y, z

  real, parameter :: third = 1.0 / 3.0
!
!
!     ...Start.
!
!
  if (present (printInfo)) then
      doPrint = printInfo .and. present (printUnit)
  end if

  if (doPrint) then
      write (printUnit,printFormat) ' initial quartic q3    = ',q3
      write (printUnit,printFormat) ' initial quartic q2    = ',q2
      write (printUnit,printFormat) ' initial quartic q1    = ',q1
      write (printUnit,printFormat) ' initial quartic q0    = ',q0
      write (printUnit,printFormat) ' ------------------------------------------------'
  end if
!
!
!     ...Handle special cases. Since the cubic solver handles all its
!        special cases by itself, we need to check only for two cases:
!
!            1) independent term is zero -> solve cubic and include
!               the zero root
!
!            2) the biquadratic case.
!
!
  if (q0 == 0.0) then

      k  = 1.0
      a3 = q3
      a2 = q2
      a1 = q1

      quarticType = 3

  else if (q3 == 0.0 .and. q1 == 0.0) then

      k  = 1.0
      a2 = q2
      a0 = q0

      quarticType = 2

  else
!
!
!     ...The general case. Rescale quartic polynomial, such that largest absolute coefficient
!        is (exactly!) equal to 1. Honor the presence of a special quartic case that might have
!        been obtained (due to underflow in the coefficients).
!
!
      s = abs (q3)
      t = sqrt (abs (q2))
      u = abs (q1) ** third
      x = sqrt (sqrt (abs (q0)))
      y = max (s,t,u,x)

      if (y == s) then

          k  = 1.0 / s
          a3 = sign (1.0, q3)
          a2 = (q2 * k) * k
          a1 = ((q1 * k) * k) * k
          a0 = (((q0 * k) * k) * k) * k

      else if (y == t) then

          k  = 1.0 / t
          a3 = q3 * k
          a2 = sign (1.0, q2)
          a1 = ((q1 * k) * k) * k
          a0 = (((q0 * k) * k) * k) * k

      else if (y == u) then

          k  = 1.0 / u
          a3 = q3 * k
          a2 = (q2 * k) * k
          a1 = sign (1.0, q1)
          a0 = (((q0 * k) * k) * k) * k

      else

          k  = 1.0 / x
          a3 = q3 * k
          a2 = (q2 * k) * k
          a1 = ((q1 * k) * k) * k
          a0 = sign (1.0, q0)

      end if

      k = 1.0 / k

      if (doPrint) then
          write (printUnit,printFormat) ' rescaling factor      = ',k
          write (printUnit,printFormat) ' ------------------------------------------------'
          write (printUnit,printFormat) ' rescaled quartic q3   = ',a3
          write (printUnit,printFormat) ' rescaled quartic q2   = ',a2
          write (printUnit,printFormat) ' rescaled quartic q1   = ',a1
          write (printUnit,printFormat) ' rescaled quartic q0   = ',a0
          write (printUnit,printFormat) ' ------------------------------------------------'
      end if

      if (a0 == 0.0) then
          quarticType = 3
      else if (a3 == 0.0 .and. a1 == 0.0) then
          quarticType = 2
      else
          quarticType = 4
      end if

  end if
!
!
!     ...Select the case.
!
!        1) The quartic with independent term = 0 -> solve cubic and add a zero root.
!
!
  select case (quarticType)

    case (3)

      call Roots_x3Polynomial (a3, a2, a1, nReal, root (1:3,1:2), printInfo, printUnit)

      if (nReal == 3) then

          x = root (1,Re) * k       ! real roots of cubic are ordered x >= y >= z
          y = root (2,Re) * k
          z = root (3,Re) * k

          nReal = 4

          root (1,Re) = max (x, 0.0)
          root (2,Re) = max (y, min (x, 0.0))
          root (3,Re) = max (z, min (y, 0.0))
          root (4,Re) = min (z, 0.0)
          root (:,Im) = 0.0

      else                          ! there is only one real cubic root here

          x = root (1,Re) * k

          nReal = 2

          root (4,Re) = root (3,Re) * k
          root (3,Re) = root (2,Re) * k
          root (2,Re) = min (x, 0.0)
          root (1,Re) = max (x, 0.0)

          root (4,Im) = root (3,Im) * k
          root (3,Im) = root (2,Im) * k
          root (2,Im) = 0.0
          root (1,Im) = 0.0

      end if
!
!
!     ...2) The quartic with x^3 and x terms = 0 -> solve biquadratic.
!
!
    case (2)

      call Roots_x2Polynomial (q2, q0, nReal, root (1:2,1:2))

      if (nReal == 2) then

          x = root (1,Re)         ! real roots of quadratic are ordered x >= y
          y = root (2,Re)

          if (y >= 0.0) then

              x = sqrt (x) * k
              y = sqrt (y) * k

              nReal = 4

              root (1,Re) = x
              root (2,Re) = y
              root (3,Re) = - y
              root (4,Re) = - x
              root (:,Im) = 0.0

          else if (x >= 0.0 .and. y < 0.0) then

              x = sqrt (x)       * k
              y = sqrt (abs (y)) * k

              nReal = 2

              root (1,Re) = x
              root (2,Re) = - x
              root (3,Re) = 0.0
              root (4,Re) = 0.0
              root (1,Im) = 0.0
              root (2,Im) = 0.0
              root (3,Im) = y
              root (4,Im) = - y

          else if (x < 0.0) then

              x = sqrt (abs (x)) * k
              y = sqrt (abs (y)) * k

              nReal = 0

              root (:,Re) = 0.0
              root (1,Im) = y
              root (2,Im) = x
              root (3,Im) = - x
              root (4,Im) = - y

          end if

      else          ! complex conjugate pair biquadratic roots x +/- iy.
              
          x = root (1,Re) * 0.5
          y = root (1,Im) * 0.5
          z = sqrt (x * x + y * y)
          y = sqrt (z - x) * k
          x = sqrt (z + x) * k

          nReal = 0

          root (1,Re) = x
          root (2,Re) = x
          root (3,Re) = - x
          root (4,Re) = - x
          root (1,Im) = y
          root (2,Im) = - y
          root (3,Im) = y
          root (4,Im) = - y

      end if
!
!
!     ...3) The general quartic case. Search for stationary points. Set the first
!           derivative polynomial (cubic) equal to zero and find its roots.
!           Check, if any minimum point of Q(x) is below zero, in which case we
!           must have real roots for Q(x). Hunt down only the real root, which
!           will potentially converge fastest during Newton iterates. The remaining
!           roots will be determined by deflation Q(x) -> cubic.
!
!           The best roots for the Newton iterations are the two on the opposite
!           ends, i.e. those closest to the +2 and -2. Which of these two roots
!           to take, depends on the location of the Q(x) minima x = s and x = u,
!           with s > u. There are three cases:
!
!              1) both Q(s) and Q(u) < 0
!                 ----------------------
!
!                 The best root to target is the one that has the lowest derivative
!                 value at the initial root guess. Let x be the initial root guess
!                 for Q(s): +2 (or zero, if s < 0 and a0 > 0). Likewise, let y be
!                 the initial root guess for Q(u): -2 (or zero, if u > 0 and a0 > 0).
!                 We then calculate Q'(x) and Q'(y). Ideally we should have always
!                 Q'(x) > 0 and Q'(y) < 0. To be safe we operate with their absolute
!                 values and choose as follows: if |Q'(x)| =< |Q'(y)| we pick x as
!                 our initial Newton guess, if |Q'(y)| < |Q'(x)| we pick y as our
!                 initial Newton guess. The idea behind this is that the lower the
!                 first derivative, the faster we approach the root through Newton
!                 iterations.
!
!              2) only Q(s) < 0
!                 -------------
!
!                 With both sides +2 and -2 possible as a Newton starting point,
!                 we have to avoid the area in the Q(x) graph, where inflection
!                 points are present. Solving Q''(x) = 0, leads to solutions
!                 x = -a3/4 +/- discriminant, i.e. they are centered around -a3/4.
!                 Since both inflection points must be either on the r.h.s or l.h.s.
!                 from x = s, a simple test where s is in relation to -a3/4 allows
!                 us to avoid the inflection point area.
!
!              3) only Q(u) < 0
!                 -------------
!
!                 Same of what has been said under 2) but with x = u.
!
!           On rare occasions, the chosen Q(s) or Q(u) are < 0, but the Newton
!           iterates result in oscillations where Q(x) is always > 0. This happens
!           if Q(s) or Q(u) were < 0 just by chance due to rounding errors in their
!           evaluation. In this case we need to store the root corresponding to Q(s)
!           or Q(u) < 0 as a lower bound for the Bisection method that follows the
!           Newton procedure in case of Newton oscillations.
!
!
    case (4)

      x = 0.75 * a3
      y = 0.50 * a2
      z = 0.25 * a1

      if (doPrint) then
          write (printUnit,printFormat) ' dQ(x)/dx cubic c2     = ',x
          write (printUnit,printFormat) ' dQ(x)/dx cubic c1     = ',y
          write (printUnit,printFormat) ' dQ(x)/dx cubic c0     = ',z
          write (printUnit,printFormat) ' ------------------------------------------------'
      end if

      call Roots_x3Polynomial (x, y, z, nReal, root (1:3,1:2), printInfo, printUnit)

      s = root (1,Re)        ! Q'(x) root s (real for sure)
      x = s + a3
      x = x * s + a2
      x = x * s + a1
      x = x * s + a0         ! Q(s)

      y = 1.0                ! dual info: Q'(x) has more real roots, and if so, is Q(u) < 0 ? 

      if (nReal > 1) then
          u = root (3,Re)    ! Q'(x) root u
          y = u + a3
          y = y * u + a2
          y = y * u + a1
          y = y * u + a0     ! Q(u)
      end if

      if (doPrint) then
          write (printUnit,printFormat) ' dQ(x)/dx root s       = ',s
          write (printUnit,printFormat) ' Q(s)                  = ',x
          write (printUnit,printFormat) ' dQ(x)/dx root u       = ',u
          write (printUnit,printFormat) ' Q(u)                  = ',y
          write (printUnit,printFormat) ' ------------------------------------------------'
      end if

      if (x < 0.0 .and. y < 0.0) then

          if (s < 0.0 .and. a0 > 0.0) then
              x = 0.0
          else
              x = 2.0
          end if

          if (u > 0.0 .and. a0 > 0.0) then
              y = 0.0
          else
              y = - 2.0
          end if

          a = x + a3
          b = x + a
          a = a * x + a2
          b = b * x + a
          a = a * x + a1
          b = b * x + a     ! b = Q'(x)

          c = y + a3
          d = y + c
          c = c * y + a2
          d = d * y + c
          c = c * y + a1
          d = d * y + c     ! d = Q'(y)

          if (doPrint) then
              write (printUnit,printFormat) ' dQ(x)/dx at root x    = ',b
              write (printUnit,printFormat) ' dQ(x)/dx at root y    = ',d
              write (printUnit,printFormat) ' ------------------------------------------------'
          end if

          if (abs (b) > abs (d)) then   ! if Q'(y) < Q'(x),
              x = y                     ! take root u for Newton iterations
              s = u                     ! save for lower bisecion bound just in case
          end if

          nReal = 1

      else if (x < 0.0) then

          if (s < - a3 * 0.25) then
              if (s > 0.0 .and. a0 > 0.0) then
                  x = 0.0
              else
                  x = - 2.0
              end if
          else
              if (s < 0.0 .and. a0 > 0.0) then
                  x = 0.0
              else
                  x = 2.0
              end if
          end if

          nReal = 1

      else if (y < 0.0) then

          if (u < - a3 * 0.25) then
              if (u > 0.0 .and. a0 > 0.0) then
                  x = 0.0
              else
                  x = - 2.0
              end if
          else
              if (u < 0.0 .and. a0 > 0.0) then
                  x = 0.0
              else
                  x = 2.0
              end if
          end if
          s = u                 ! save for lower bisecion bound just in case

          nReal = 1
      else
          nReal = 0
      end if
!
!
!     ...Do all necessary Newton iterations. In case we have more than 2 oscillations,
!        exit the Newton iterations and switch to bisection. Note, that from the
!        definition of the Newton starting point, we always have Q(x) > 0 and Q'(x)
!        starts (-ve/+ve) for the (-2/+2) starting points and (increase/decrease) smoothly
!        and staying (< 0 / > 0). In practice, for extremely shallow Q(x) curves near the
!        root, the Newton procedure can overshoot slightly due to rounding errors when
!        approaching the root. The result are tiny oscillations around the root. If such
!        a situation happens, the Newton iterations are abandoned after 3 oscillations
!        and further location of the root is done using bisection starting with the
!        oscillation brackets.
!
!
      if (nReal > 0) then

          t = 32.0                                          ! maximum possible Q(x) value
          oscillate = 0
          bisection = .false.
          converged = .false.

          do while (.not.converged .and. .not.bisection)    ! Newton-Raphson iterates

             y = x + a3                                     !
             z = x + y                                      !
             y = y * x + a2                                 ! y = Q(x)
             z = z * x + y                                  !
             y = y * x + a1                                 ! z = Q'(x)
             z = z * x + y                                  !
             y = y * x + a0                                 !

             if (y > t) then                                ! does Newton start oscillating ?
                 oscillate = oscillate + 1                  ! increment oscillation counter
             end if

             if (y < 0.0) then
                 s = x                                      ! save lower bisection bound
             else
                 u = x                                      ! save upper bisection bound
             end if

             if (z == 0.0) then                             ! safeguard against accidental
                 bisection = .true.                         ! Q'(x) = 0 due to roundoff
                 exit                                       ! errors -> bisection takes
             end if                                         ! over

             t = y                                          ! save Q(x) for next Newton step
             y = y / z                                      ! Newton correction
             x = x - y                                      ! new Newton root

             bisection = oscillate > 2                      ! activate bisection
             converged = abs (y) <= abs (x) * rt_macheps    ! Newton convergence indicator

             if (doPrint) write (printUnit,printFormat)     ' Newton root           = ',x

          end do

          if (bisection) then

              t = u - s                                     ! initial bisection interval
              do while (abs (t) > abs (x) * rt_macheps)     ! bisection iterates

                 y = x + a3                                 !
                 y = y * x + a2                             ! y = Q(x)
                 y = y * x + a1                             !
                 y = y * x + a0                             !

                 if (y < 0.0) then                          !
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
!
!
!     ...Find remaining roots -> reduce to cubic. The reduction to a cubic polynomial
!        is done using composite deflation to minimize rounding errors. Also, while
!        the composite deflation analysis is done on the reduced quartic, the actual
!        deflation is being performed on the original quartic again to avoid enhanced
!        propagation of root errors.
!
!
          z = abs (x)            !
          a = abs (a3)           !
          b = abs (a2)           ! prepare for composite deflation
          c = abs (a1)           !
          d = abs (a0)           !

          y = z * max (a,z)      ! take maximum between |x^2|,|a3 * x|

          deflateCase = 1        ! up to now, the maximum is |x^4| or |a3 * x^3|

          if (y < b) then        ! check maximum between |x^2|,|a3 * x|,|a2|
              y = b * z          ! the maximum is |a2| -> form |a2 * x|
              deflateCase = 2    ! up to now, the maximum is |a2 * x^2|
          else
              y = y * z          ! the maximum is |x^3| or |a3 * x^2|
          end if

          if (y < c) then        ! check maximum between |x^3|,|a3 * x^2|,|a2 * x|,|a1|
              y = c * z          ! the maximum is |a1| -> form |a1 * x|
              deflateCase = 3    ! up to now, the maximum is |a1 * x|
          else
              y = y * z          ! the maximum is |x^4|,|a3 * x^3| or |a2 * x^2|
          end if

          if (y < d) then        ! check maximum between |x^4|,|a3 * x^3|,|a2 * x^2|,|a1 * x|,|a0|
              deflateCase = 4    ! the maximum is |a0|
          end if

          x = x * k              ! 1st real root of original Q(x)

          select case (deflateCase)

          case (1)
            z = 1.0 / x
            u = - q0 * z         ! u -> backward deflation on original Q(x)
            t = (u - q1) * z     ! t -> backward deflation on original Q(x)
            s = (t - q2) * z     ! s -> backward deflation on original Q(x)
          case (2)
            z = 1.0 / x
            u = - q0 * z         ! u -> backward deflation on original Q(x)
            t = (u - q1) * z     ! t -> backward deflation on original Q(x)
            s = q3 + x           ! s ->  forward deflation on original Q(x)
          case (3)
            s = q3 + x           ! s ->  forward deflation on original Q(x)
            t = q2 + s * x       ! t ->  forward deflation on original Q(x)
            u = - q0 / x         ! u -> backward deflation on original Q(x)
          case (4)
            s = q3 + x           ! s ->  forward deflation on original Q(x)
            t = q2 + s * x       ! t ->  forward deflation on original Q(x)
            u = q1 + t * x       ! u ->  forward deflation on original Q(x)
          end select

          if (doPrint) then
              write (printUnit,printFormat) ' Residual cubic c2     = ',s
              write (printUnit,printFormat) ' Residual cubic c1     = ',t
              write (printUnit,printFormat) ' Residual cubic c0     = ',u
              write (printUnit,printFormat) ' ------------------------------------------------'
          end if

          call Roots_x3Polynomial (s, t, u, nReal, root (1:3,1:2), printInfo, printUnit)

          if (nReal == 3) then

              s = root (1,Re)    !
              t = root (2,Re)    ! real roots of cubic are ordered s >= t >= u
              u = root (3,Re)    !

              root (1,Re) = max (s, x)
              root (2,Re) = max (t, min (s, x))
              root (3,Re) = max (u, min (t, x))
              root (4,Re) = min (u, x)
              root (:,Im) = 0.0

              nReal = 4

          else                   ! there is only one real cubic root here

              s = root (1,Re)

              root (4,Re) = root (3,Re)
              root (3,Re) = root (2,Re)
              root (2,Re) = min (s, x)
              root (1,Re) = max (s, x)
              root (4,Im) = root (3,Im)
              root (3,Im) = root (2,Im)
              root (2,Im) = 0.0
              root (1,Im) = 0.0

              nReal = 2

          end if

      else
!
!
!     ...If no real roots have been found by now, only complex roots are possible.
!        Find real parts of roots first, followed by imaginary components.
!
!
          s = a3 * 0.5
          t =  s * s - a2
          u =  s * t + a1                     ! value of Q'(-a3/4) at stationary point -a3/4

          notZero = (abs (u) >= rt_macheps)   ! H(-a3/4) is considered > 0 at stationary point

          if (doPrint) then
              write (printUnit,printFormat) ' dQ/dx (-a3/4) value   = ',u
              write (printUnit,printFormat) ' ------------------------------------------------'
          end if

          if (a3 /= 0.0) then
              s = a1 / a3
              minimum = (a0 > s * s)
          else
              minimum = (4 * a0 > a2 * a2)                      ! H''(-a3/4) > 0 -> minimum
          end if

          iterate = notZero .or. (.not.notZero .and. minimum)

          if (iterate) then

              x = sign (2.0, a3)                                ! initial root -> target = smaller mag root

              oscillate = 0
              bisection = .false.
              converged = .false.

              do while (.not.converged .and. .not.bisection)    ! Newton-Raphson iterates

                 a = x + a3                                     !
                 b = x + a                                      ! a = Q(x)
                 c = x + b                                      !
                 d = x + c                                      ! b = Q'(x)
                 a = a * x + a2                                 !
                 b = b * x + a                                  ! c = Q''(x) / 2
                 c = c * x + b                                  !
                 a = a * x + a1                                 ! d = Q'''(x) / 6
                 b = b * x + a                                  !
                 a = a * x + a0                                 !
                 y = a * d * d - b * c * d + b * b              ! y = H(x), usually < 0
                 z = 2 * d * (4 * a - b * d - c * c)            ! z = H'(x)

                 if (y > 0.0) then                              ! does Newton start oscillating ?
                     oscillate = oscillate + 1                  ! increment oscillation counter
                     s = x                                      ! save upper bisection bound
                 else
                     u = x                                      ! save lower bisection bound
                 end if

                 y = y / z                                      ! Newton correction
                 x = x - y                                      ! new Newton root

                 bisection = oscillate > 2                      ! activate bisection
                 converged = abs (y) <= abs (x) * rt_macheps    ! Newton convergence criterion

                 if (doPrint) write (printUnit,printFormat)     ' Newton H(x) root      = ',x

              end do

              if (bisection) then

                  t = u - s                                     ! initial bisection interval
                  do while (abs (t) > abs (x * rt_macheps))     ! bisection iterates

                     a = x + a3                                 !
                     b = x + a                                  ! a = Q(x)
                     c = x + b                                  !
                     d = x + c                                  ! b = Q'(x)
                     a = a * x + a2                             !
                     b = b * x + a                              ! c = Q''(x) / 2
                     c = c * x + b                              !
                     a = a * x + a1                             ! d = Q'''(x) / 6
                     b = b * x + a                              !
                     a = a * x + a0                             !
                     y = a * d * d - b * c * d + b * b          ! y = H(x)

                     if (y > 0.0) then                          !
                         s = x                                  !
                     else                                       ! keep bracket on root
                         u = x                                  !
                     end if                                     !

                     t = 0.5 * (u - s)                          ! new bisection interval
                     x = s + t                                  ! new bisection root

                     if (doPrint) write (printUnit,printFormat) ' Bisection H(x) root   = ',x

                  end do
              end if

              if (doPrint) write (printUnit,printFormat) ' ------------------------------------------------'

              a = x * k                                         ! 1st real component -> a
              b = - 0.5 * q3 - a                                ! 2nd real component -> b

              x = 4 * a + q3                                    ! Q'''(a)
              y = x + q3 + q3                                   !
              y = y * a + q2 + q2                               ! Q'(a)
              y = y * a + q1                                    !
              y = max (y / x, 0.0)                              ! ensure Q'(a) / Q'''(a) >= 0
              x = 4 * b + q3                                    ! Q'''(b)
              z = x + q3 + q3                                   !
              z = z * b + q2 + q2                               ! Q'(b)
              z = z * b + q1                                    !
              z = max (z / x, 0.0)                              ! ensure Q'(b) / Q'''(b) >= 0
              c = a * a                                         ! store a^2 for later
              d = b * b                                         ! store b^2 for later
              s = c + y                                         ! magnitude^2 of (a + iy) root
              t = d + z                                         ! magnitude^2 of (b + iz) root

              if (s > t) then                                   ! minimize imaginary error
                  c = sqrt (y)                                  ! 1st imaginary component -> c
                  d = sqrt (q0 / s - d)                         ! 2nd imaginary component -> d
              else
                  c = sqrt (q0 / t - c)                         ! 1st imaginary component -> c
                  d = sqrt (z)                                  ! 2nd imaginary component -> d
              end if

          else                                                  ! no bisection -> real components equal

              a = - 0.25 * q3                                   ! 1st real component -> a
              b = a                                             ! 2nd real component -> b = a

              x = a + q3                                        !
              x = x * a + q2                                    ! Q(a)
              x = x * a + q1                                    !
              x = x * a + q0                                    !
              y = - 0.1875 * q3 * q3 + 0.5 * q2                 ! Q''(a) / 2
              z = max (y * y - x, 0.0)                          ! force discriminant to be >= 0
              z = sqrt (z)                                      ! square root of discriminant
              y = y + sign (z,y)                                ! larger magnitude root
              x = x / y                                         ! smaller magnitude root
              c = max (y, 0.0)                                  ! ensure root of biquadratic > 0
              d = max (x, 0.0)                                  ! ensure root of biquadratic > 0
              c = sqrt (c)                                      ! large magnitude imaginary component
              d = sqrt (d)                                      ! small magnitude imaginary component

          end if

          if (a > b) then

              root (1,Re) = a
              root (2,Re) = a
              root (3,Re) = b
              root (4,Re) = b
              root (1,Im) = c
              root (2,Im) = - c
              root (3,Im) = d
              root (4,Im) = - d

          else if (a < b) then

              root (1,Re) = b
              root (2,Re) = b
              root (3,Re) = a
              root (4,Re) = a
              root (1,Im) = d
              root (2,Im) = - d
              root (3,Im) = c
              root (4,Im) = - c

          else

              root (1,Re) = a
              root (2,Re) = a
              root (3,Re) = a
              root (4,Re) = a
              root (1,Im) = c
              root (2,Im) = - c
              root (3,Im) = d
              root (4,Im) = - d

          end if

      end if    ! # of real roots 'if'

  end select    ! quartic type select
!
!
!     ...Ready!
!
!
  return
end subroutine Roots_x4Polynomial
