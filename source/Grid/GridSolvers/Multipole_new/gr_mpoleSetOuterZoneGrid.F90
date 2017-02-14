!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleSetOuterZoneGrid
!!
!! NAME
!!
!!  gr_mpoleSetOuterZoneGrid
!!
!! SYNOPSIS
!!
!!  gr_mpoleSetOuterZoneGrid ()
!!
!! DESCRIPTION
!!
!!  This routine sets up the outer (statistical) zone radial grid.
!!
!!***

subroutine gr_mpoleSetOuterZoneGrid ()

  use gr_mpoleData,    ONLY : gr_mpoleDr,                    &
                              gr_mpoleDrInv,                 &
                              gr_mpoleMaxR,                  &
                              gr_mpoleMaxQ,                  &
                              gr_mpoleMaxRadialZones,        &
                              gr_mpoleMinRadialZone,         &
                              gr_mpoleZoneRmax,              &
                              gr_mpoleZoneQmax,              &
                              gr_mpoleZoneType,              &
                              gr_mpoleZoneScalar,            &
                              gr_mpoleZoneLogNorm,           &
                              gr_mpoleZoneExponent,          &
                              gr_mpoleZoneScalarInv,         &
                              gr_mpoleZoneLogNormInv,        &
                              gr_mpoleZoneExponentInv,       &
                              gr_mpoleZoneMaxRadiusFraction, &
                              gr_mpoleInnerZoneQmax,         &
                              gr_mpoleInnerZoneMaxR,         &
                              gr_mpoleOuterZoneQshift

  implicit none

#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"

  integer :: outerZoneQmin
  integer :: Q,Qlocal,Qglobal
  integer :: type
  integer :: zone

  real    :: rglobal,rlocal
  real    :: scalar,lognorm,exponent
  real    :: sclInv,lgnInv,expInv
!
!
!       ...Set up the characteristic values for each specified radial zone
!          according to the linear and scaling factors given for each of these
!          zones. These characteristic arrays consist of:
!
!                a) maximum global radial values for each zone
!                b) maximum global bin values for each zone 
!
!          Global here means including all the corresponding values of lower
!          zones, as opposed to local values, which are specific to each zone.
!          Set also the overall maximum number of radial bins to be expected.
!
!          The outer (statistical) zones can be of two types: exponential or
!          logarithmic. Both cases use different equations to determine the
!          radial bin boundaries:
!
!               Exponential  =>  maximum r for local Q = s * dr * Q^t
!               Logarithmic  =>  maximum r for local Q = s * dr * [e^(Qt) - 1]/[e^t - 1]
!
!          A note is in place for the logarithmic scaling. It can easily be shown
!          that:
!
!                  [e^(Qt) - 1]/[e^t - 1] >= Q  for all t >= 0 and Q >= 1
!
!          The equality holds only for the two cases:
!
!                                 1) limit t --> 0
!                                 2)       Q  =  1
!
!          The proof of this is by induction starting with the Q = 2 expression for
!          the ratio. Why not use a different logarithmic scaling formula like
!          s * dr * Q * e^(t[Q-1]) ? The reason is that this formula is very hard
!          to invert for Q.
!
!          The pair of scalar/exponential values (s,t) and to which type each zone
!          belongs to can been specified by the user through runtime parameters. 
!          To determine the maximum local Q value for each zone, the above equations
!          must be inverted:
!
!            Exponential  =>  maximum local Q for an r = [(r/(s*dr))^(1/t)]
!            Logarithmic  =>  maximum local Q for an r = [(1/t)*log{r*(e^t-1)/(s*dr) + 1}]
!
!          where [] denotes the ceiling function.
!
!
  gr_mpoleZoneRmax (0) = ZERO
  gr_mpoleZoneQmax (0) = ZERO

  do zone = 1,gr_mpoleMaxRadialZones

     type    = gr_mpoleZoneType        (zone)
     sclInv  = gr_mpoleZoneScalarInv   (zone)
     expInv  = gr_mpoleZoneExponentInv (zone)
     rglobal = gr_mpoleZoneMaxRadiusFraction (zone) * gr_mpoleMaxR
     rlocal  = rglobal - gr_mpoleZoneRmax (zone-1)

     if (type == ZONE_EXPONENTIAL) then
         Qlocal = ceiling ( (rlocal * sclInv * gr_mpoleDrInv) ** expInv )
     else if (type == ZONE_LOGARITHMIC) then
         lgnInv = gr_mpoleZoneLogNormInv (zone)
         Qlocal = ceiling ( expInv * log (rlocal * sclInv * gr_mpoleDrInv * lgnInv + ONE) )
     end if

     Qglobal = gr_mpoleZoneQmax (zone-1) + Qlocal

     gr_mpoleZoneRmax (zone) = rglobal
     gr_mpoleZoneQmax (zone) = Qglobal

  end do
!
!
!       ...Find the first significant outer zone. All outer zones below
!          the first significant outer zone have been 'swallowed' by the
!          size of the inner zone. If no first significant outer zone
!          is found, i.e. if the complete domain is within the inner zone,
!          the first significant outer zone will be equal to ZERO for
!          further processing (see below).
!
!
  gr_mpoleMinRadialZone = 0

  do zone = gr_mpoleMaxRadialZones,1,-1
     rglobal = gr_mpoleZoneRmax (zone)
     if (rglobal > gr_mpoleInnerZoneMaxR) then
         gr_mpoleMinRadialZone = zone
     end if
  end do
!
!
!       ...Determine within the first significant outer zone, which global
!          radial bin (outerZoneQmin) first exceeds the maximum inner zone bin.
!          This gap has to be closed when merging the inner and outer zones
!          so that at their boundary the radial bin numbers continue without
!          interruption. Note, that the gap can be positive or negative:
!
!                   +ve : outerZoneQmin > gr_mpoleInnerZoneQmax
!                   -ve : outerZoneQmin < gr_mpoleInnerZoneQmax
!
!          In either case a radial bin shift value (gr_mpoleOuterZoneQshift) is
!          calculated, which brings the two zones together. This shift value
!          will be added to the outer zone radial bin values later on when
!          evaluating the moment bins.
!
!
  if (gr_mpoleMinRadialZone == 0) then

      gr_mpoleMaxQ = gr_mpoleInnerZoneQmax

  else

      Qlocal   = gr_mpoleZoneQmax     (gr_mpoleMinRadialZone) - gr_mpoleZoneQmax (gr_mpoleMinRadialZone - 1)
      type     = gr_mpoleZoneType     (gr_mpoleMinRadialZone)
      scalar   = gr_mpoleZoneScalar   (gr_mpoleMinRadialZone)
      exponent = gr_mpoleZoneExponent (gr_mpoleMinRadialZone)

      do Q = 1,Qlocal

         if (type == ZONE_EXPONENTIAL) then
             rlocal = scalar * gr_mpoleDr * (real (Q) ** exponent)
         else if (type == ZONE_LOGARITHMIC) then
             lognorm = gr_mpoleZoneLogNorm (gr_mpoleMinRadialZone)
             rlocal  = scalar * gr_mpoleDr * lognorm * (exp (exponent * real (Q)) - ONE)
         end if

         rglobal = gr_mpoleZoneRmax (gr_mpoleMinRadialZone - 1) + rlocal

         if (rglobal > gr_mpoleInnerZoneMaxR) then
             outerZoneQmin = gr_mpoleZoneQmax (gr_mpoleMinRadialZone - 1) + Q
             exit
         end if

      end do

      gr_mpoleOuterZoneQshift = gr_mpoleInnerZoneQmax - outerZoneQmin + 1

      gr_mpoleMaxQ = gr_mpoleZoneQmax (gr_mpoleMaxRadialZones) + gr_mpoleOuterZoneQshift

  end if
!
!
!       ...Ready!
!
!
  return
end subroutine gr_mpoleSetOuterZoneGrid
