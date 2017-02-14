!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/2Dellipse/sim_distancePoint2Ellipse2Dxy
!!
!! NAME
!!
!!  sim_distancePoint2Ellipse2Dxy
!!
!! SYNOPSIS
!!
!!  sim_distancePoint2Ellipse2Dxy (real,    intent (in)  :: a,
!!                                 real,    intent (in)  :: b,
!!                                 real,    intent (in)  :: Rot,
!!                                 real,    intent (in)  :: Cx,
!!                                 real,    intent (in)  :: Cy,
!!                                 real,    intent (in)  :: Px,
!!                                 real,    intent (in)  :: Py)
!!
!! DESCRIPTION
!!
!!  Calculates the minimum/maximum distance from a point (Px,Py) to a general ellipse in 2D
!!  cartesian domain. The general ellipse is characterized by the following parameters:
!!
!!            a   ->  length of major semi-axis
!!            b   ->  length of minor semi-axis
!!            Rot ->  angle (degrees) of clockwise rotation of the minor axis wrt the y-axis
!!            Cx  ->  the ellipse center x-coordinate
!!            Cy  ->  the ellipse center y-coordinate
!!
!!  which positions the ellipse in a unique way in the 2D cartesian domain. The parametric
!!  equation of this ellipse is:
!!
!!              | x(t) |   |   cos (Rot)  sin (Rot) |   | a sin (t) |   | Cx |
!!              |      | = |                        | * |           | + |    |
!!              | y(t) |   | - sin (Rot)  cos (Rot) |   | b cos (t) |   | Cy |
!!
!!  where 0 <= t <= 2pi. Derivation of the minimum distance conditions (vector from the
!!  point (Px,Py) to the ellipse curve must the perpendicular to the tangent vector at
!!  the ellipse point) leads to the following equation to be solved:
!!
!!               R * cos (t) + S * sin (t) + U * cos (t) * sin (t) = 0
!!
!!  where:
!!
!!               R = a  * { (Py - Cy) sin (Rot) - (Px - Cx) cos (Rot) }
!!               S = b  * { (Px - Cx) sin (Rot) + (Py - Cy) cos (Rot) }
!!               U = a^2 - b^2
!!
!!  The U term will be > 0 for an ellipse and = 0 for a circle. The elliptical U term
!!  is responsible for the 'cos (t) * sin (t)' term in the above mixed trigonometric
!!  equation and makes finding the solutions in 't' hard. After 't' has been found,
!!  the distance is given by the expression:
!!
!!            distance = sqrt { (Px - x(t))^2 + (Py - y(t))^2 }
!!
!!  There are 2 ways to solve this mixed trigonometric equation:
!!
!!
!!      1) Substitution of either sin or cos. This leads to quartic equations in terms
!!         of cos (t) or sin (t). For example, substituting sin (t) = sqrt (1 - cos^2 (t))
!!         leads to:
!!
!!                  cos^4 (t) + A cos^3 (t) + B cos^2 (t) + C cos(t) + D = 0
!!
!!         with:
!!
!!                    A = 2 (S/U)
!!                    B = (R/U)^2 + (S/U)^2 - 1
!!                    C = - 2 (S/U)
!!                    D = - (S/U)^2
!!
!!         Drawback: Besides the need to solve a quartic, up to 8 solutions have to be
!!                   checked: 4 solutions for cos (t) and each will lead to 2 possible
!!                   sin (t) values via sin (t) = +/- sqrt (1 - cos^2 (t)).
!!
!!         Advantage: The solution of the quartic gives all possible solutions, i.e.
!!                    they contain both the minimum and maximum distance info. Also
!!                    the quartic solver in FLASH is now pretty mature and very stable.
!!
!!
!!      2) Newton-Raphson (NR) iterations, starting from a suitable initial 't'. Here one
!!         has to be careful in choosing a good initial guess.
!!
!!         Drawback: Not easy to choose an initial good guess. The mixed trigonometric
!!                   equation can have up to 4 't' roots and since we don't know which
!!                   of these roots will lead to the minimum or maximum distance, one
!!                   would have to find all roots -> very challenging, as each root
!!                   cannot be factored out as for example for polynomials.
!!
!!         Advantage: If the aforementioned drawbacks can be overcome and the NR can
!!                    be made stable by choosing good initial guesses, then the four sets
!!                    of NR iterations would probably be more economical than solving the
!!                    quartic, although the quartic solver in FLASH is already highly
!!                    optimized.
!!
!!
!!  In the present code we implemented the first approach based on solving the quartic
!!  equation in cos (t). Instead of setting up the equations in the rotated framework,
!!  we use the equations of the unrotated framework, as then the minor/major axis of the
!!  ellipse coincides with the y/x axis of the domain. This procedure eliminates unnecessary
!!  multiplications by the sine and cosine of the rotation angle. The only stage at which
!!  we need to apply the rotation is at the beginning when we need to rotate BOTH the
!!  point (Px,Py) and the ellipse center (Cx,Cy) backwards to the corresponding points
!!  (Px',Py') and (Cx',Cy'). In the unrotated framework the R and S constants become:
!!
!!                                R = - a  * (Px' - Cx')
!!                                S =   b  * (Py' - Cy')
!!
!!  and the unrotated ellipse parametric equation is:
!!
!!                          | x(t) |   | a sin (t) |   | Cx' |
!!                          |      | = |           | + |     |
!!                          | y(t) |   | b cos (t) |   | Cy' |
!!
!!
!! ARGUMENTS
!!
!!  a        : length of major semi-axis
!!  b        : length of minor semi-axis
!!  Rot      : angle (degrees) of clockwise rotation of the minor axis wrt the y-axis
!!  Cx       : the ellipse center x-coordinate
!!  Cy       : the ellipse center y-coordinate
!!  Px       : the x-coordinate of the point (Px,Py)
!!  Py       : the y-coordinate of the point (Px,Py)
!!
!! NOTES
!!
!!  1) This is an array real function with two elements:
!!
!!        sim_distancePoint2Ellipse2Dxy (1) : minimum distance
!!        sim_distancePoint2Ellipse2Dxy (2) : maximum distance
!!
!!  2) The code works also if each of the major/minor semiaxis a/b length is negative or
!!     zero. The following cases can be handled by the code:
!!
!!        i) |a| > 0 and b = 0   (line in space of length 2|a|)
!!       ii)   a = 0 and b = 0   (point in space)
!!
!!***

function sim_distancePoint2Ellipse2Dxy (a, b, Rot, Cx, Cy, Px, Py)

  use Driver_interface,  ONLY : Driver_abortFlash
  use Roots_interface,   ONLY : Roots_x4Polynomial
  use Simulation_data,   ONLY : sim_deg2rad

  implicit none

  real, intent (in)  :: a, b
  real, intent (in)  :: Rot
  real, intent (in)  :: Cx, Cy
  real, intent (in)  :: Px, Py

  logical :: circle

  integer :: i,j,k
  integer :: nReal

  real    :: PCx, PCy
  real    :: U
  real    :: sim_distancePoint2Ellipse2Dxy (1:2)
  real    :: p (1:8,1:2)
!
!
!     ...Establish needed constants. Decide early on, if the major/minor ellipse axes are flipped
!        or if the difference between these axes is so tiny, that the code thinks the ellipse
!        is a circle.
!
!
  U  = a * a - b * b

  if (U < 0.0) then
      call Driver_abortFlash ('[sim_distancePoint2Ellipse2Dxy] ERROR: ellipse major axis < minor axis!')
  end if

  circle = (U == 0.0)

  if (circle) then

      PCx = Px - Cx  
      PCy = Py - Cy

      p (1,1) = sqrt (PCx * PCx + PCy * PCy)

      sim_distancePoint2Ellipse2Dxy (1) = abs (p (1,1) - a)
      sim_distancePoint2Ellipse2Dxy (2) =      p (1,1) + a

  else

      p (1,1) = Rot * sim_deg2rad                      ! rotation angle Rot in radians
      p (2,1) = cos (p (1,1))                          ! cos (Rot)
      p (3,1) = sin (p (1,1))                          ! sin (Rot)

      PCx = p (2,1) * (Px - Cx) - p (3,1) * (Py - Cy)  ! rotate Px - Cx backwards
      PCy = p (3,1) * (Px - Cx) + p (2,1) * (Py - Cy)  ! rotate Py - Cy backwards

      p (5,1) = - a * PCx / U                          ! this is R/U
      p (6,1) = + b * PCy / U                          ! this is S/U
      p (5,2) = p (6,1) + p (6,1)                      ! this is 2(S/U)
      p (8,2) = - p (6,1) * p (6,1)                    ! this is - (S/U)^2
      p (6,2) = p (5,1) * p (5,1) - p (8,2) - 1.0      ! this is (R/U)^2 + (S/U)^2 - 1
      p (7,2) = - p (5,2)                              ! this is - 2 (S/U)

      call  Roots_x4Polynomial (p (5,2), p (6,2), p (7,2), p (8,2),     nReal, p (1:4,1:2))
!
!
!     ...The nReal (2 or 4) quartic real roots correspond to the nReal cos (t) values.
!        Each of these give rise to 2 x nReal corresponding sin (t) values via:
!
!                           sin (t) = +/- (1 - cos^2 (t))
!
!        These represent 2 x nReal points on the ellipse curve, of which only 2 lead to
!        the minimum and maximum distance. Also, due to the numerically intensive quartic
!        root determination, the values obtained for cos (t) from the quartic solver can
!        over/undershoot the range [-1,+1] permissible for the cosine values. The code
!        below adjusts for this, preventing a potential sin (t) evaluation crash.
!
!
!      write (*,*) ' nReal = ',nReal                                                  !
!      write (*,*) ' '                                                                !
!      write (*,*) ' Quartic Solver Root 1 (x + iy) = ',p (1,1),' + ',p (1,2),' i'    ! for debugging
!      write (*,*) ' Quartic Solver Root 2 (x + iy) = ',p (2,1),' + ',p (2,2),' i'    !
!      write (*,*) ' Quartic Solver Root 3 (x + iy) = ',p (3,1),' + ',p (3,2),' i'    ! purpose only
!      write (*,*) ' Quartic Solver Root 4 (x + iy) = ',p (4,1),' + ',p (4,2),' i'    !
!      write (*,*) ' '                                                                !
!
!
      if (nReal == 2 .or. nReal == 4) then

          i = nReal
          j = nReal + 1
          k = nReal + nReal

          p (1:i,1) = min ( 1.0, p (1:i,1))                    ! ensure cos (t) <= +1
          p (1:i,1) = max (-1.0, p (1:i,1))                    ! ensure cos (t) >= -1
          p (1:i,2) = a * sqrt (1.0 - p (1:i,1) * p (1:i,1))   ! generate the a * sin (t)
          p (1:i,1) = PCy - b * p (1:i,1)                      ! generate the PCy - b * cos (t)
          p (1:i,1) = p (1:i,1) * p (1:i,1)                    ! generate the (PCy - b * cos (t))^2
          p (j:k,1) = p (1:i,1)                                ! store extra (PCy - b * cos (t))^2
          p (j:k,2) = PCx - p (1:i,2)                          ! generate the PCx - a * sin (t)
          p (1:i,2) = PCx + p (1:i,2)                          ! generate the PCx + a * sin (t)
          p (1:k,2) = p (1:k,2) * p (1:k,2)                    ! generate all (PCx +/- a * sin (t))^2
          p (1:k,1) = p (1:k,1) + p (1:k,2)                    ! generate all square distances

          sim_distancePoint2Ellipse2Dxy (1) = sqrt (minval (p (1:k,1)))   ! minimum distance
          sim_distancePoint2Ellipse2Dxy (2) = sqrt (maxval (p (1:k,1)))   ! maximum distance

      else
          call Driver_abortFlash ('[sim_distancePoint2Ellipse2Dxy] ERROR: no/bad # of real quartic roots!')
      end if

  end if
!
!
!     ...Ready!
!
!
  return
end function sim_distancePoint2Ellipse2Dxy
