!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleInit
!!
!! NAME
!!
!!  gr_mpoleInit
!!
!! SYNOPSIS
!!
!!  gr_mpoleInit ()
!!
!! DESCRIPTION
!!
!!  Initialize the multipole Poisson solver.  Read in any of the
!!  runtime parameters for this solver.
!!
!!***

subroutine gr_mpoleInit ()

  use gr_mpoleData

  use Grid_data,                   ONLY : gr_geometry
  use Driver_interface,            ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface,           ONLY : Logfile_stamp
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"

  character*30 :: labelExponent
  character*30 :: labelFraction
  character*30 :: labelScalar
  character*30 :: labelType
  character*11 :: typeOfZone
  character*3  :: zoneNumber

  logical :: angularYGrid
  logical :: angularZGrid
  logical :: grid1D
  logical :: grid2D
  logical :: grid3D
  logical :: gridCartesian
  logical :: gridCylindrical
  logical :: gridSpherical
  logical :: gridPolar
  logical :: radialXGrid

  integer :: n
  integer :: status
  integer :: zone

  real    :: angularYRange
  real    :: angularZRange
  real    :: lastZoneFraction
!
!
!    ...Get the external parameters. Catch any 'nonsense' parameters.
!
!
  call RuntimeParameters_get  ("mpole_MultiThreading",       gr_mpoleMultiThreading     )
  call RuntimeParameters_get  ("mpole_Lmax",                 gr_mpoleMaxL               )
  call RuntimeParameters_get  ("mpole_2DSymmetryPlane",      gr_mpoleSymmetryPlane2D    )
  call RuntimeParameters_get  ("mpole_3DAxisymmetry",        gr_mpoleSymmetryAxis3D     )
  call RuntimeParameters_get  ("mpole_DumpMoments",          gr_mpoleMomentsDump        )
  call RuntimeParameters_get  ("mpole_PrintRadialInfo",      gr_mpoleRadialInfoPrint    )
  call RuntimeParameters_get  ("mpole_IgnoreInnerZone",      gr_mpoleIgnoreInnerZone    )
  call RuntimeParameters_get  ("mpole_MaxRadialZones",       gr_mpoleMaxRadialZones     )
  call RuntimeParameters_get  ("mpole_InnerZoneResolution",  gr_mpoleInnerZoneResolution)
  call RuntimeParameters_get  ("mpole_InnerZoneSize",        gr_mpoleInnerZoneSize      )
  call RuntimeParameters_get  ("xmin",                       gr_mpoleDomainXmin         )
  call RuntimeParameters_get  ("xmax",                       gr_mpoleDomainXmax         )
  call RuntimeParameters_get  ("ymin",                       gr_mpoleDomainYmin         )
  call RuntimeParameters_get  ("ymax",                       gr_mpoleDomainYmax         )
  call RuntimeParameters_get  ("zmin",                       gr_mpoleDomainZmin         )
  call RuntimeParameters_get  ("zmax",                       gr_mpoleDomainZmax         )

  if (gr_mpoleMaxRadialZones <= 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: no radial zones specified')
  end if

  if (gr_mpoleInnerZoneSize <= 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: no radial inner zone specified')
  end if

  if (gr_mpoleInnerZoneResolution <= ZERO) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: radial inner zone resolution value must be > 0')
  end if

  if (gr_mpoleInnerZoneResolution > 0.1) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: radial inner zone resolution value must be =< 0.1')
  end if
!
!
!    ...A bunch of local constants needed by other routines. Note, that we calculate
!       Pi explicitely here rather than taking it from the pre-defined constants list.
!       This is done in order to ensure highest possible accuracy. A pre-defined
!       constant is only as accurate as it is defined.
!
!
  gr_mpolePi        = acos (-ONE)
  gr_mpoleTwoPi     = gr_mpolePi + gr_mpolePi
  gr_mpoleHalfPi    = HALF * gr_mpolePi
  gr_mpoleThirdPi   = gr_mpolePi / THREE
  gr_mpoleSixthPi   = HALF * gr_mpoleThirdPi
  gr_mpoleFourPi    = gr_mpoleTwoPi + gr_mpoleTwoPi
  gr_mpoleFourPiInv = ONE / gr_mpoleFourPi
  gr_mpoleEbase     = exp (ONE)
  gr_mpoleEbaseInv  = ONE / gr_mpoleEbase
  gr_mpoleRad2Deg   = 180./ gr_mpolePi

  gr_mpoleInnerZoneResolutionInv = ONE / gr_mpoleInnerZoneResolution
!
!
!    ...Create a handle to the current geometry.
!
!
  grid3D = (NDIM == 3)
  grid2D = (NDIM == 2)
  grid1D = (NDIM == 1)

  gridCartesian   = (gr_geometry == CARTESIAN)
  gridCylindrical = (gr_geometry == CYLINDRICAL)
  gridSpherical   = (gr_geometry == SPHERICAL)
  gridPolar       = (gr_geometry == POLAR)

  if (grid3D .and. gridCartesian)   gr_mpoleGeometry = GRID_3DCARTESIAN
  if (grid2D .and. gridCartesian)   gr_mpoleGeometry = GRID_2DCARTESIAN
  if (grid1D .and. gridCartesian)   gr_mpoleGeometry = GRID_1DCARTESIAN
  if (grid3D .and. gridCylindrical) gr_mpoleGeometry = GRID_3DCYLINDRICAL
  if (grid2D .and. gridCylindrical) gr_mpoleGeometry = GRID_2DCYLINDRICAL
  if (grid1D .and. gridCylindrical) gr_mpoleGeometry = GRID_1DCYLINDRICAL
  if (grid3D .and. gridSpherical)   gr_mpoleGeometry = GRID_3DSPHERICAL
  if (grid2D .and. gridSpherical)   gr_mpoleGeometry = GRID_2DSPHERICAL
  if (grid1D .and. gridSpherical)   gr_mpoleGeometry = GRID_1DSPHERICAL
  if (grid3D .and. gridPolar)       gr_mpoleGeometry = GRID_3DPOLAR
  if (grid2D .and. gridPolar)       gr_mpoleGeometry = GRID_2DPOLAR
  if (grid1D .and. gridPolar)       gr_mpoleGeometry = GRID_1DPOLAR
!
!
!    ...Before proceeding, catch the unsupported geometries (not needed or not
!       yet implemented). Also inform the user (and abort the program), if symmetry
!       requirements cannot be honored or make no sense.
!
!
  if (gr_mpoleGeometry == GRID_2DCARTESIAN   .or. &
      gr_mpoleGeometry == GRID_1DCARTESIAN   .or. &
      gr_mpoleGeometry == GRID_1DCYLINDRICAL .or. &
      gr_mpoleGeometry == GRID_3DSPHERICAL   .or. &
      gr_mpoleGeometry == GRID_3DPOLAR       .or. &
      gr_mpoleGeometry == GRID_2DPOLAR       .or. &
      gr_mpoleGeometry == GRID_1DPOLAR) then

      call Driver_abortFlash ('[gr_mpoleInit] ERROR: unsupported geometry')
  end if

  if (.not.grid2D .and. gr_mpoleSymmetryPlane2D) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: plane of symmetry only allowed in 2D geometry')
  end if

  if (gr_mpoleGeometry == GRID_3DCYLINDRICAL .and. gr_mpoleSymmetryAxis3D) then
      call Logfile_stamp ('3D axis symmetry in 3D cylindrical geometry equivalent to 2D cylindrical run!', &
                          '[gr_mpoleInit]')
      call Driver_abortFlash ('[gr_mpoleInit] SUGGESTION: 3D axis symmetry! -> Rerun in 2D cylindrical!')
  end if
!
!
!    ...Check the given domain boundaries. In case we have a geometry with an underlying
!       angular grid, the FLASH convention is that the angular boundaries are set internally
!       in terms of radians. Any value different from 0 or pi/2pi cannot be handled. Also, when
!       specifying a symmetry plane in 2D, we want to make sure that the lower left boundary
!       corresponds to the center of the problem. Likewise, all radial coordinates should
!       start from 0.
!
!
  radialXGrid  = (gridCylindrical .or. &
                  gridSpherical   .or. &
                  gridPolar)

  angularYGrid = (gr_mpoleGeometry == GRID_2DSPHERICAL   .or. &
                  gr_mpoleGeometry == GRID_3DSPHERICAL   .or. &
                  gr_mpoleGeometry == GRID_2DPOLAR       .or. &
                  gr_mpoleGeometry == GRID_3DPOLAR)

  angularZGrid = (gr_mpoleGeometry == GRID_3DCYLINDRICAL .or. &
                  gr_mpoleGeometry == GRID_3DSPHERICAL)

  angularYRange = gr_mpolePi
  angularZRange = gr_mpoleTwoPi

  if (gr_mpoleGeometry == GRID_2DSPHERICAL .and. gr_mpoleSymmetryPlane2D) then
      angularYRange = HALF * gr_mpolePi
  end if

  if (gr_mpoleGeometry == GRID_2DPOLAR .or. gr_mpoleGeometry == GRID_3DPOLAR) then
      angularYRange = gr_mpoleTwoPi
  end if

  if (radialXGrid) then
      if (gr_mpoleDomainXmin /= ZERO) then
          call Driver_abortFlash ('[gr_mpoleInit] ERROR: radial X-coord out of range')
      end if
  end if

  if (angularYGrid) then
      if (gr_mpoleDomainYmin /= ZERO .or. gr_mpoleDomainYmax /= angularYRange) then
          call Driver_abortFlash ('[gr_mpoleInit] ERROR: angular Y-coord out of range')
      end if
  end if

  if (angularZGrid) then
      if (gr_mpoleDomainZmin /= ZERO .or. gr_mpoleDomainZmax /= angularZRange) then
          call Driver_abortFlash ('[gr_mpoleInit] ERROR: angular Z-coord out of range')
      end if
  end if

  if (gr_mpoleSymmetryPlane2D) then
      if (gr_mpoleDomainXmin /= ZERO .or. gr_mpoleDomainYmin /= ZERO) then
          call Driver_abortFlash ('[gr_mpoleInit] ERROR: 2D symmetry plane requires xmin,ymin = 0')
      end if
  end if
!
!
!  ...Set the internal domain limits, if radial/angular geometries are present.
!     (--> improves readability of the code, i.e. avoids the confusing FLASH
!      convention naming).
!
!
  select case (gr_mpoleGeometry)

     case (GRID_3DCYLINDRICAL)

           gr_mpoleDomainRmin     = gr_mpoleDomainXmin
           gr_mpoleDomainPhiMin   = gr_mpoleDomainZmin     ! order is important here
           gr_mpoleDomainZmin     = gr_mpoleDomainYmin     ! otherwise we loose Z info
           gr_mpoleDomainRmax     = gr_mpoleDomainXmax
           gr_mpoleDomainPhiMax   = gr_mpoleDomainZmax     ! order is important here
           gr_mpoleDomainZmax     = gr_mpoleDomainYmax     ! otherwise we loose Z info

     case (GRID_2DCYLINDRICAL)

           gr_mpoleDomainRmin     = gr_mpoleDomainXmin
           gr_mpoleDomainZmin     = gr_mpoleDomainYmin
           gr_mpoleDomainRmax     = gr_mpoleDomainXmax
           gr_mpoleDomainZmax     = gr_mpoleDomainYmax

     case (GRID_2DSPHERICAL)

           gr_mpoleDomainRmin     = gr_mpoleDomainXmin
           gr_mpoleDomainThetaMin = gr_mpoleDomainYmin
           gr_mpoleDomainRmax     = gr_mpoleDomainXmax
           gr_mpoleDomainThetaMax = gr_mpoleDomainYmax

     case (GRID_1DSPHERICAL)

           gr_mpoleDomainRmin     = gr_mpoleDomainXmin
           gr_mpoleDomainRmax     = gr_mpoleDomainXmax

  end select
!
!
!  ...Printout some L,M restriction info (if applicable).
!
!
  select case (gr_mpoleGeometry)

     case (GRID_3DCARTESIAN)

       if (gr_mpoleSymmetryAxis3D) then
           call Logfile_stamp ('3D axissymmetry, ignoring M > 0 moments','[gr_mpoleInit]')
       end if

     case (GRID_2DCYLINDRICAL)
           call Logfile_stamp ('2D cylindrical, ignoring M > 0 moments','[gr_mpoleInit]')
     case (GRID_2DSPHERICAL)
           call Logfile_stamp ('2D spherical, ignoring M > 0 moments','[gr_mpoleInit]')
     case (GRID_1DSPHERICAL)
           call Logfile_stamp ('1D spherical, ignoring L > 0 moments','[gr_mpoleInit]')
  end select
!
!
!  ...Set the M and combined LM dimensions from the given maximum L.
!     Override maximum L only for the 1D spherical case.
!     Set values of constants depending on this data.
!
!
  select case (gr_mpoleGeometry)

     case (GRID_3DCARTESIAN)

       if (gr_mpoleSymmetryAxis3D) then

           gr_mpoleMaxM  = 0
           gr_mpoleMaxLM = gr_mpoleMaxL + 1
           gr_mpoleTotalNrCosineMoments = gr_mpoleMaxLM

       else

           gr_mpoleMaxM  = gr_mpoleMaxL
           gr_mpoleMaxLM = (gr_mpoleMaxL + 1) * (gr_mpoleMaxL + 1)
           gr_mpoleTotalNrCosineMoments = (gr_mpoleMaxL * (gr_mpoleMaxL + 1) / 2) + gr_mpoleMaxL + 1

       end if

     case (GRID_3DCYLINDRICAL)

           gr_mpoleMaxM  = gr_mpoleMaxL
           gr_mpoleMaxLM = (gr_mpoleMaxL + 1) * (gr_mpoleMaxL + 1)
           gr_mpoleTotalNrCosineMoments = (gr_mpoleMaxL * (gr_mpoleMaxL + 1) / 2) + gr_mpoleMaxL + 1

     case (GRID_2DCYLINDRICAL , GRID_2DSPHERICAL)

           gr_mpoleMaxM  = 0
           gr_mpoleMaxLM = gr_mpoleMaxL + 1
           gr_mpoleTotalNrCosineMoments = gr_mpoleMaxLM

     case (GRID_1DSPHERICAL)

           gr_mpoleMaxM  = 0
           gr_mpoleMaxL  = 0
           gr_mpoleMaxLM = 1
           gr_mpoleTotalNrCosineMoments = 1

     case default
           call Driver_abortFlash ('[gr_mpoleInit] PROGRAMMER ERROR 2: We should never be here!')
  end select

  gr_mpoleMax2L = gr_mpoleMaxL + gr_mpoleMaxL
!
!
!    ...Allocate the arrays related to the leaf blocks.
!
!
  allocate (gr_mpoleBlockList (1:MAXBLOCKS), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleBlockList allocate failed')
  end if
!
!
!    ...Allocate other needed data structures.
!
!
  if (gr_mpoleMaxL > 0) then

      allocate (gr_mpoleNumberInv (1:gr_mpoleMax2L), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleNumberInv allocate failed')
      end if

  end if
!
!
!    ...Allocate the radial zone arrays.
!
!
  allocate (gr_mpoleZoneRmax (0:gr_mpoleMaxRadialZones), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleZoneRmax allocate failed')
  end if

  allocate (gr_mpoleZoneQmax (0:gr_mpoleMaxRadialZones), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleZoneQmax allocate failed')
  end if

  allocate (gr_mpoleZoneType (1:gr_mpoleMaxRadialZones), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleZoneType allocate failed')
  end if

  allocate (gr_mpoleZoneScalar (1:gr_mpoleMaxRadialZones), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleZoneScalar allocate failed')
  end if

  allocate (gr_mpoleZoneLogNorm (1:gr_mpoleMaxRadialZones), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleZoneLogNorm allocate failed')
  end if

  allocate (gr_mpoleZoneExponent (1:gr_mpoleMaxRadialZones), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleZoneExponent allocate failed')
  end if

  allocate (gr_mpoleZoneScalarInv (1:gr_mpoleMaxRadialZones), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleZoneScalarInv allocate failed')
  end if

  allocate (gr_mpoleZoneLogNormInv (1:gr_mpoleMaxRadialZones), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleZoneLogNormInv allocate failed')
  end if

  allocate (gr_mpoleZoneExponentInv (1:gr_mpoleMaxRadialZones), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleZoneExponentInv allocate failed')
  end if

  allocate (gr_mpoleZoneMaxRadiusFraction (1:gr_mpoleMaxRadialZones), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[gr_mpoleInit] ERROR: gr_mpoleZoneMaxRadiusFraction allocate failed')
  end if
!
!
!    ...Fill the needed arrays.
!
!
  do n = 1,gr_mpoleMax2L
     gr_mpoleNumberInv (n) = ONE / real (n)
  end do

  do zone = 1,gr_mpoleMaxRadialZones

     write (zoneNumber,'(I3)') zone

     labelFraction = "mpole_ZoneRadiusFraction_"//adjustl (zoneNumber)
     labelScalar   = "mpole_ZoneScalar_"        //adjustl (zoneNumber)
     labelExponent = "mpole_ZoneExponent_"      //adjustl (zoneNumber)
     labelType     = "mpole_ZoneType_"          //adjustl (zoneNumber)

     call RuntimeParameters_get  (labelFraction,  gr_mpoleZoneMaxRadiusFraction (zone))
     call RuntimeParameters_get  (labelScalar,    gr_mpoleZoneScalar            (zone))
     call RuntimeParameters_get  (labelExponent,  gr_mpoleZoneExponent          (zone))
     call RuntimeParameters_get  (labelType,      typeOfZone                          )
     
     gr_mpoleZoneLogNorm     (zone) = ONE / (exp (gr_mpoleZoneExponent (zone)) - ONE)
     gr_mpoleZoneScalarInv   (zone) = ONE / gr_mpoleZoneScalar   (zone)
     gr_mpoleZoneLogNormInv  (zone) = ONE / gr_mpoleZoneLogNorm  (zone)
     gr_mpoleZoneExponentInv (zone) = ONE / gr_mpoleZoneExponent (zone)

     if (typeOfZone == 'exponential') then
         gr_mpoleZoneType (zone) = ZONE_EXPONENTIAL
     else if (typeOfZone == 'logarithmic') then
         gr_mpoleZoneType (zone) = ZONE_LOGARITHMIC
     else
         call Driver_abortFlash ('[gr_mpoleInit] ERROR: unvalid radial zone type')
     end if

  end do
!
!
!    ...Check for bad radial zone data.
!
!
  lastZoneFraction = gr_mpoleZoneMaxRadiusFraction (gr_mpoleMaxRadialZones)

  if (lastZoneFraction /= ONE) then
      gr_mpoleZoneMaxRadiusFraction (gr_mpoleMaxRadialZones) = ONE
      call Logfile_stamp('last radial zone fraction reset to 1','[gr_mpoleInit]')
  end if

  do zone = 2,gr_mpoleMaxRadialZones
     if ( gr_mpoleZoneMaxRadiusFraction (zone) < gr_mpoleZoneMaxRadiusFraction (zone-1) ) then
          call Driver_abortFlash ('[gr_mpoleInit] ERROR: radial zone fractions out of order!')
     end if
  end do

  do zone = 1,gr_mpoleMaxRadialZones
     if ( gr_mpoleZoneScalar (zone) <= ZERO ) then
          call Driver_abortFlash ('[gr_mpoleInit] ERROR: radial zone scalar =< 0 !')
     end if
  end do

  do zone = 1,gr_mpoleMaxRadialZones
     if ( gr_mpoleZoneExponent (zone) <= ZERO ) then
          call Driver_abortFlash ('[gr_mpoleInit] ERROR: radial zone exponent =< 0 !')
     end if
  end do
!
!
!       Done.
!
!
  return
end subroutine gr_mpoleInit
