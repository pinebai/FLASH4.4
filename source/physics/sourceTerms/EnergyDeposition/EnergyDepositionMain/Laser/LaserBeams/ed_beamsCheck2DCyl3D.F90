!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beamsCheck2DCyl3D
!!
!! NAME
!!
!!  ed_beamsCheck2DCyl3D
!!
!! SYNOPSIS
!!
!!  call ed_beamsCheck2DCyl3D ()
!!
!! DESCRIPTION
!!
!!  Checks the collected 3D cartesian beams data for 2D cylindrical grids (cartesian). This is
!!  a rather specialized routine, in which the dimensions where the beams are defined is different
!!  from the dimension of the underlaying grid. All checks which depend on domain grid details
!!  should go in here. Currently it contains the following:
!!
!!     1) Check, if all beam elliptical 3D target areas are completely within the
!!        3D cylindrical domain implicated by the 2D cylindrical grid.
!!
!!     2) Check, if all beam elliptical 3D lens areas are completely outside the
!!        3D cylindrical domain implicated by the 2D cylindrical grid.
!!
!!  The 3D cylindrical domain is obtained from the (x=R,z) 2D cylindrical grid by a 360 degree
!!  revolution around the 2D cylindrical z-axis.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsCheck2DCyl3D ()

  use Driver_interface,         ONLY : Driver_abortFlash

  use Roots_interface,          ONLY : Roots_x4Polynomial
 
  use EnergyDeposition_data,    ONLY : ed_beams,               &
                                       ed_numberOfBeams,       &
                                       ed_orthogonalTolerance, &
                                       ed_xminDomain,          &
                                       ed_xmaxDomain,          &
                                       ed_yminDomain,          &
                                       ed_ymaxDomain
  
  implicit none

#include "Flash.h"
#include "EnergyDeposition.h"
#include "constants.h"

  logical :: circular
  logical :: CisZero, DisZero
  logical :: inDomain
  logical :: outOfDomain
  logical :: useCosine

  integer :: beam
  integer :: n
  integer :: nReal
  integer :: nTrig

  integer, parameter :: Re = 1

  real    :: A, B, C, D
  real    :: aa, bb, cc, dd, ac, ad, bc, bd
  real    :: cosine, sine
  real    :: Ex1, Ex2, Ey1, Ey2
  real    :: ez
  real    :: LdotL, s1dotL, s2dotL
  real    :: lensRmin, lensZmin, lensZmax
  real    :: Lx, Ly, Lz
  real    :: maxRsq, minRsq
  real    :: omega
  real    :: q0, q1, q2, q3, q3c, q3s, q4, q4inv
  real    :: root1, root2, root3, root4
  real    :: Rsq1, Rsq2
  real    :: s1dots1, s1dots2, s2dots2
  real    :: s1L, s2L, s1T, s2T
  real    :: s1x, s2x, s1y, s2y
  real    :: targetRmax, targetZmin, targetZmax
  real    :: TdotT, s1dotT, s2dotT
  real    :: Tx, Ty, Tz
  real    :: u1x, u1y, u1z, u2x, u2y, u2z

  real    :: Trig (1:4)

  real    :: root (1:4,1:2)
!
!
!     ...Loop over all 3D beams.
!
!
  do beam = 1, ed_numberOfBeams

     u1x = ed_beams (beam) % semiAxisUnitMajorX
     u1y = ed_beams (beam) % semiAxisUnitMajorY
     u1z = ed_beams (beam) % semiAxisUnitMajorZ
     u2x = ed_beams (beam) % semiAxisUnitMinorX
     u2y = ed_beams (beam) % semiAxisUnitMinorY
     u2z = ed_beams (beam) % semiAxisUnitMinorZ
     s1T = ed_beams (beam) % targetSemiAxisMajor
     s2T = ed_beams (beam) % targetSemiAxisMinor
     s1L = ed_beams (beam) % lensSemiAxisMajor
     s2L = ed_beams (beam) % lensSemiAxisMinor
     Tx  = ed_beams (beam) % targetX
     Ty  = ed_beams (beam) % targetY
     Tz  = ed_beams (beam) % targetZ
     Lx  = ed_beams (beam) % lensX
     Ly  = ed_beams (beam) % lensY
     Lz  = ed_beams (beam) % lensZ
!
!
!     ...We first deal with the (non-radial) z-component of the elliptical target area.
!        Calculate the minimum and maximum target area extension in z-direction.
!
!
     ez = sqrt ( (s1T * u1z) ** 2 + (s2T * u2z) ** 2 )

     targetZmin = Tz - ez
     targetZmax = Tz + ez
!
!
!     ...Next deal with the (x,y) part of the elliptical target area. Form the 2D projections
!        onto the (x,y) plane and get the relevant info.
!
!
     s1x = s1T * u1x
     s1y = s1T * u1y
     s2x = s2T * u2x
     s2y = s2T * u2y

     TdotT   =  Tx * Tx  +  Ty * Ty
     s1dotT  = s1x * Tx  + s1y * Ty
     s2dotT  = s2x * Tx  + s2y * Ty
     s1dots1 = s1x * s1x + s1y * s1y
     s1dots2 = s1x * s2x + s1y * s2y
     s2dots2 = s2x * s2x + s2y * s2y

     C = s1dots1 - s2dots2
     D = s1dots2

     CisZero = (abs (C) / max (s1dots1, s2dots2)) < ed_orthogonalTolerance
     DisZero = (abs (D) / max (s1dots1, s2dots2)) < ed_orthogonalTolerance

     circular = DisZero .and. CisZero
!
!
!     ...Set up the quartic coefficients and solve the appropriate quartic. Handle separately
!        the case when the projected 2D target ellipse is a circle. For the target area we are
!        only interested in knowing if it is contained in the domain, so only the maximum
!        radial value is of interest.
!
!
     if (circular) then

         targetRmax = sqrt (TdotT) + sqrt (s1dots1)      ! add lengths of projected T and s1

     else

         A = s1dotT
         B = s2dotT

         aa = A * A
         bb = B * B
         cc = C * C
         dd = D * D
         ac = A * C
         ad = A * D
         bc = B * C
         bd = B * D

         q4    = cc + dd + dd + dd + dd
         q4inv = 1.0 / q4
         q3c   = bd + bd + bd + bd + ac + ac
         q3s   = ad + ad + ad + ad - bc - bc

         useCosine = abs (q3c) > abs (q3s)

         if (useCosine) then
             q3 = q3c * q4inv
             q2 = (aa + bb) * q4inv - 1.0
             q1 = - (ac + ac + bd + bd) * q4inv
             q0 = (dd - aa) * q4inv
         else
             q3 = q3s * q4inv
             q2 = (aa + bb) * q4inv - 1.0
             q1 = (bc + bc - ad - ad) * q4inv
             q0 = (dd - bb) * q4inv
         end if

         call Roots_x4Polynomial (q3,q2,q1,q0,   nReal,root)

         nTrig = 0

         root1 = root (1,Re)
         root2 = root (2,Re)
         root3 = root (3,Re)
         root4 = root (4,Re)

         if (nReal > 0 .and. abs (root1) <= 1.0) then
             nTrig = nTrig + 1
             Trig (nTrig) = root1
         end if

         if (nReal > 1 .and. abs (root2) <= 1.0) then 
             nTrig = nTrig + 1
             Trig (nTrig) = root2
         end if

         if (nReal > 2 .and. abs (root3) <= 1.0) then 
             nTrig = nTrig + 1
             Trig (nTrig) = root3
         end if

         if (nReal > 3 .and. abs (root4) <= 1.0) then 
             nTrig = nTrig + 1
             Trig (nTrig) = root4
         end if

         if (nTrig == 0) then
             call Driver_abortFlash ("ed_beamsCheck2DCyl3D: Quartic root failure for target!")
         end if
!
!
!     ...For each trigonometric solution calculate the corresponding angle and calculate
!        the corresponding radial distances. We have to distinguish two cases:
!
!            1) we use the cosine solution(s) of the quartic. In this case there are
!               2 possible angles in the range 0 - 2pi. Each of these angles corresponds
!               to different sine values, which are related by a minus sign. Both sine
!               values are tested for extrema.
!
!            2) we use the sine solution(s) of the quartic. In this case there are also
!               2 possible angles in the range 0 - 2pi. Each of these angles corresponds
!               to different cosine values, which are related by a minus sign. Both cosine
!               values are tested for extrema.
!
!        Note, that we are not using the cos(x)^2 + sin(x)^2 = 1 formula to get the other
!        complementary trig function. If any one of cos(x)^2 or sin(x)^2 is close to 1,
!        errors for the other trig function become large and the trig function is inaccurate.
!
!
!
         maxRsq = 0.0

         do n = 1,nTrig

            if (useCosine) then

                cosine = Trig (n)
                omega  = acos (cosine)                       ! 0 - pi is range of omega
                sine   =  sin (omega)

                Ex1 = Tx + s1x * cosine + s2x * sine         ! the 1st sine solution: + sine
                Ey1 = Ty + s1y * cosine + s2y * sine
                Ex2 = Tx + s1x * cosine - s2x * sine         ! the 2nd sine solution: - sine
                Ey2 = Ty + s1y * cosine - s2y * sine

            else

                sine   = Trig (n)
                omega  = asin (sine)                         ! - pi/2 to + pi/2 is range of omega
                cosine =  cos (omega)

                Ex1 = Tx + s1x * cosine + s2x * sine         ! the 1st cosine solution: + cosine
                Ey1 = Ty + s1y * cosine + s2y * sine
                Ex2 = Tx - s1x * cosine + s2x * sine         ! the 2nd cosine solution: - cosine
                Ey2 = Ty - s1y * cosine + s2y * sine

            end if

            Rsq1 = Ex1 * Ex1 + Ey1 * Ey1
            Rsq2 = Ex2 * Ex2 + Ey2 * Ey2

            maxRsq = max (maxRsq , Rsq1, Rsq2)

         end do

         targetRmax = sqrt (maxRsq)

     end if
!
!
!     ...Check, if beam target area is contained completely within domain.
!
!
     inDomain =      (targetZmin >= ed_yminDomain) &
               .and. (targetZmax <= ed_ymaxDomain) &
               .and. (targetRmax <= ed_xmaxDomain)

     if (.not.inDomain) then
         call Driver_abortFlash ("ed_beamsCheck2DCyl3D: Beam target (partially) outside of domain!")
     end if
!
!
!     ...Proceed with the lens area. The (non-radial) z-component of the elliptical lens area
!        is treated first. Calculate the minimum and maximum lens area extension in z-direction.
!
!
     ez = sqrt ( (s1L * u1z) ** 2 + (s2L * u2z) ** 2 )

     lensZmin = Lz - ez
     lensZmax = Lz + ez
!
!
!     ...Handle the (x,y) part of the elliptical lens area. Form the 2D projections onto
!        the (x,y) plane and get the relevant info.
!
!
     s1x = s1L * u1x
     s1y = s1L * u1y
     s2x = s2L * u2x
     s2y = s2L * u2y

     LdotL   =  Lx * Lx  +  Ly * Ly
     s1dotL  = s1x * Lx  + s1y * Ly
     s2dotL  = s2x * Lx  + s2y * Ly
     s1dots1 = s1x * s1x + s1y * s1y
     s1dots2 = s1x * s2x + s1y * s2y
     s2dots2 = s2x * s2x + s2y * s2y

     C = s1dots1 - s2dots2
     D = s1dots2

     CisZero = (abs (C) / max (s1dots1, s2dots2)) < ed_orthogonalTolerance
     DisZero = (abs (D) / max (s1dots1, s2dots2)) < ed_orthogonalTolerance

     circular = DisZero .and. CisZero
!
!
!     ...Set up the quartic coefficients and solve the appropriate quartic. Handle separately
!        the case when the projected 2D lens ellipse is a circle. For the lens area we are
!        only interested in knowing if it is outside of the domain, so only the minimum
!        radial value is of interest.
!
!
     if (circular) then

         lensRmin = sqrt (LdotL) - sqrt (s1dots1)   ! subtract lengths of projected s1 from projected L

     else

         A = s1dotL
         B = s2dotL

         aa = A * A
         bb = B * B
         cc = C * C
         dd = D * D
         ac = A * C
         ad = A * D
         bc = B * C
         bd = B * D

         q4    = cc + dd + dd + dd + dd
         q4inv = 1.0 / q4
         q3c   = bd + bd + bd + bd + ac + ac
         q3s   = ad + ad + ad + ad - bc - bc

         useCosine = abs (q3c) > abs (q3s)

         if (useCosine) then
             q3 = q3c * q4inv
             q2 = (aa + bb) * q4inv - 1.0
             q1 = - (ac + ac + bd + bd) * q4inv
             q0 = (dd - aa) * q4inv
         else
             q3 = q3s * q4inv
             q2 = (aa + bb) * q4inv - 1.0
             q1 = (bc + bc - ad - ad) * q4inv
             q0 = (dd - bb) * q4inv
         end if

         call Roots_x4Polynomial (q3,q2,q1,q0,   nReal,root)

         nTrig = 0

         root1 = root (1,Re)
         root2 = root (2,Re)
         root3 = root (3,Re)
         root4 = root (4,Re)

         if (nReal > 0 .and. abs (root1) <= 1.0) then
             nTrig = nTrig + 1
             Trig (nTrig) = root1
         end if

         if (nReal > 1 .and. abs (root2) <= 1.0) then 
             nTrig = nTrig + 1
             Trig (nTrig) = root2
         end if

         if (nReal > 2 .and. abs (root3) <= 1.0) then 
             nTrig = nTrig + 1
             Trig (nTrig) = root3
         end if

         if (nReal > 3 .and. abs (root4) <= 1.0) then 
             nTrig = nTrig + 1
             Trig (nTrig) = root4
         end if

         if (nTrig == 0) then
             call Driver_abortFlash ("ed_beamsCheck2DCyl3D: Quartic root failure for lens!")
         end if

         do n = 1,nTrig

            if (useCosine) then

                cosine = Trig (n)
                omega  = acos (cosine)
                sine   =  sin (omega)

                Ex1 = Lx + s1x * cosine + s2x * sine
                Ey1 = Ly + s1y * cosine + s2y * sine
                Ex2 = Lx + s1x * cosine - s2x * sine
                Ey2 = Ly + s1y * cosine - s2y * sine

            else

                sine   = Trig (n)
                omega  = asin (sine)
                cosine =  cos (omega)

                Ex1 = Lx + s1x * cosine + s2x * sine
                Ey1 = Ly + s1y * cosine + s2y * sine
                Ex2 = Lx - s1x * cosine + s2x * sine
                Ey2 = Ly - s1y * cosine + s2y * sine

            end if

            Rsq1 = Ex1 * Ex1 + Ey1 * Ey1
            Rsq2 = Ex2 * Ex2 + Ey2 * Ey2

            if (n == 1) then
                minRsq = min (Rsq1, Rsq2)
            else
                minRsq = min (minRsq , Rsq1, Rsq2)
            end if

         end do

         lensRmin = sqrt (minRsq)

     end if
!
!
!     ...Check, if beam lens area is completely outside of the domain.
!
!
     outOfDomain =     (lensZmin > ed_ymaxDomain) &
                  .or. (lensZmax < ed_yminDomain) &
                  .or. (lensRmin > ed_xmaxDomain)

     if (.not.outOfDomain) then
         call Driver_abortFlash ("ed_beamsCheck2DCyl3D: Beam lens (partially) inside the domain!")
     end if

  end do       ! loop over beams
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beamsCheck2DCyl3D
