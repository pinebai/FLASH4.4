!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleSetRadialBinData
!!
!! NAME
!!
!!  gr_mpoleSetRadialBinData
!!
!! SYNOPSIS
!!
!!  gr_mpoleSetRadialBinData ()
!!
!! DESCRIPTION
!!
!!  Sets all radial bin data that is necessary for moment and potential evaluation.
!!  It currently computes the maximum radial value for each bin and the bin
!!  damping factors. The damping factors are used to avoid possible overflow
!!  or underflow when evaluating the regular and irregular solid harmonic
!!  functions via recursion. When using these factors the recursion becomes
!!  'normalized' in the sense that the highest L-valued solid harmonic functions
!!  are close to unity.
!!
!!  For a given radial bin index Q outside of the inner zone, the maximum and
!!  minimum radius associated with that bin is calculated using the radial zone
!!  data provided. From these radii the damping factors are evaluated as follows,
!!  using the maximum angular momentum value L of the problem:
!!
!!                  regular solid harmonic = (  L-th root of L!) / rmax
!!                irregular solid harmonic = (L+1-th root of L!) / rmin
!!
!!  The x,y,z coordinates will be multiplied by the corresponding damping
!!  factors during application of the recursions. The approximations used
!!  for the roots of the factorials follow from Stirlings factorial approximation
!!  and are:
!!
!!      (  L-th root of L!) ~ (L/e) * (  2L-th root of 2*pi*L)
!!      (L+1-th root of L!) ~ (L/e) * (2L+2-th root of 2*pi*e^2/L)
!!
!!  Here e = 2.7182818...and pi = 3.14159....
!!
!!  For a given radial bin index Q inside the inner zone, there is only one
!!  radius (call it rmax) associated with it, which is sitting in the array
!!  'gr_mpoleInnerZoneDrRadii' in units of the inner zone atomic distance. The
!!  damping factors are evaluated as follows, using the maximum angular momentum
!!  value L of the problem:
!!
!!                  regular solid harmonic = (  L-th root of L!) / rmax
!!                irregular solid harmonic = (L+1-th root of L!) / rmax
!!
!!***

subroutine gr_mpoleSetRadialBinData ()

  use gr_mpoleData,   ONLY : gr_mpoleQDampingR,        &
                             gr_mpoleQDampingI,        &
                             gr_mpoleQRadii,           &
                             gr_mpoleTwoPi,            &
                             gr_mpoleEbase,            &
                             gr_mpoleEbaseInv,         &
                             gr_mpoleDr,               &
                             gr_mpoleDrInnerZone,      &
                             gr_mpoleMaxL,             &
                             gr_mpoleMaxQ,             &
                             gr_mpoleMaxRadialZones,   &
                             gr_mpoleMinRadialZone,    &
                             gr_mpoleZoneRmax,         &
                             gr_mpoleZoneQmax,         &
                             gr_mpoleZoneType,         &
                             gr_mpoleZoneScalar,       &
                             gr_mpoleZoneLogNorm,      &
                             gr_mpoleZoneExponent,     &
                             gr_mpoleInnerZoneMaxR,    &
                             gr_mpoleInnerZoneQmax,    &
                             gr_mpoleInnerZoneDrRadii, &
                             gr_mpoleOuterZoneQshift

  implicit none
  
#include "gr_mpole.h"

  real    :: exponent
  real    :: lognorm
  real    :: LthRoot,Lp1thRoot
  real    :: rglobal,rlocal
  real    :: rmin,rmax
  real    :: scalar
  real    :: xL,xLp1,xLoverEbase

  integer :: Q,Qlocal
  integer :: type
  integer :: zone
!
!
!     ...Calculate the needed root values.
!
!
  if (gr_mpoleMaxL == 0) then
      LthRoot     = ONE
      Lp1thRoot   = ONE
  else
      xL          = real (gr_mpoleMaxL)
      xLp1        = real (gr_mpoleMaxL+1)
      xLoverEbase = xL * gr_mpoleEbaseInv
      LthRoot     = xLoverEbase * ((gr_mpoleTwoPi * xL) ** (ONE/(xL+xL)))
      Lp1thRoot   = xLoverEbase * ((gr_mpoleTwoPi * gr_mpoleEbase * gr_mpoleEbase / xL) ** (ONE/(xLp1+xLp1)))
  end if
!
!
!     ...Start with the inner zone bins. If there is no inner zone, give some
!        meaningful value to the first bin (Q = 1). This is necessary, because
!        otherwise rmin will be set to ZERO and the irregular solid harmonic
!        damping factor will become infinite.
!
!
  if (gr_mpoleInnerZoneQmax > 0) then
      do Q = 1,gr_mpoleInnerZoneQmax
         rmax = gr_mpoleInnerZoneDrRadii (Q) * gr_mpoleDrInnerZone
         gr_mpoleQRadii    (Q) = rmax
         gr_mpoleQDampingR (Q) =   LthRoot / rmax
         gr_mpoleQDampingI (Q) = Lp1thRoot / rmax
      end do
  else
      rmin = gr_mpoleDr
      rmax = gr_mpoleZoneScalar (1) * gr_mpoleDr           ! we have zone = 1 and Q = 1 here
      gr_mpoleQDampingR (1) =   LthRoot / rmax
      gr_mpoleQDampingI (1) = Lp1thRoot / rmin
  end if
!
!
!     ...The rest of the bins outside the inner zone. For rmin take the
!        rmax of the previous bin (in inner zone or not).
!
!
  do Q = gr_mpoleInnerZoneQmax+1,gr_mpoleMaxQ

     do zone = gr_mpoleMinRadialZone, gr_mpoleMaxRadialZones
        if (Q - gr_mpoleZoneQmax (zone) - gr_mpoleOuterZoneQshift <= 0) exit
     end do

     Qlocal   = Q - gr_mpoleZoneQmax (zone - 1) - gr_mpoleOuterZoneQshift
     type     = gr_mpoleZoneType     (zone)
     scalar   = gr_mpoleZoneScalar   (zone)
     exponent = gr_mpoleZoneExponent (zone)

     if (type == ZONE_EXPONENTIAL) then
         rlocal  = scalar * (real (Qlocal) ** exponent) * gr_mpoleDr
     else if (type == ZONE_LOGARITHMIC) then
         lognorm = gr_mpoleZoneLogNorm (zone)
         rlocal  = scalar * scalar * gr_mpoleDr * lognorm * (exp (exponent * real (Qlocal)) - ONE)
     end if

     rglobal = gr_mpoleZoneRmax (zone - 1) + rlocal
     rmin    = rmax
     rmax    = rglobal

     gr_mpoleQRadii    (Q) = rmax
     gr_mpoleQDampingR (Q) =   LthRoot / rmax
     gr_mpoleQDampingI (Q) = Lp1thRoot / rmin

  end do
!
!
!     ...Set the needed outer edges of the radii and damping vectors.
!
!        A radius for the 0-th bin is necessary, in case there is no inner zone.
!        Evaluation of the potential then needs to determine fractional radial parts
!        of the outer zone bins. For the radial fractional parts of the Q-th bin we
!        need the radial info of the Q-th and (Q-1)-th bin.
!
!        The outer edges of the damping vectors are necessary to avoid conditional
!        if's inside the moment and potential evaluation routines in the innermost
!        loops. Note, that the edge damping vector for the regular solid harmonics is
!        set equal to the corresponding first bin one. This is done for numerical safety,
!        as the potential evaluation routines will have to evaluate
!
!                   gr_mpoleQDampingR (Q) / gr_mpoleQDampingR (Q-1)
!
!        values, which will then be taken to a power equal to the maximum angular momentum.
!        Since Damping (1) is usually a large number, for gr_mpoleMaxL > 50, the floating
!        bound can be reached if gr_mpoleQDampingR (0) is set carelessly equal to 1.
!
!
  gr_mpoleQRadii    (0)              = ZERO
  gr_mpoleQDampingR (0)              = gr_mpoleQDampingR (1)
  gr_mpoleQDampingI (gr_mpoleMaxQ+1) = ZERO
!
!
!     ...Ready!
!
!
  return
end subroutine gr_mpoleSetRadialBinData
