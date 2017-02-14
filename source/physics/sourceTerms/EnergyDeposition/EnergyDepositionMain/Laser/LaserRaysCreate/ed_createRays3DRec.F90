!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRaysCreate/ed_createRays3DRec
!!
!! NAME
!!
!!  ed_createRays3DRec
!!
!! SYNOPSIS
!!
!!  call ed_createRays3DRec (integer, intent (in) :: blockCount,
!!                           integer, intent (in) :: blockList (:),
!!                           real,    intent (in) :: timeStep,
!!                           real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates rays and places them in their initial blocks for those geometries consisting
!!  formally of 3D rectangular grids (cartesian). On exit, all rays hitting the domain boundary
!!  have been generated for the current processor. Their block IDs are not ordered as the
!!  outer loop is over all beams.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeStep       : current timestep value
!!  timeSimulation : current simulation time
!!
!! NOTES
!!
!!***

subroutine ed_createRays3DRec (blockCount, blockList, timeStep, timeSimulation)

  use Driver_interface,             ONLY : Driver_abortFlash
 
  use Grid_interface,               ONLY : Grid_getBlkBC,      &
                                           Grid_getBlkBoundBox

  use ed_interface,                 ONLY : ed_beam3DGridPointsDelta,       &
                                           ed_beam3DGridPointsRadial,      &
                                           ed_beam3DGridPointsRecBeam,     &
                                           ed_beam3DGridPointsSquare,      &
                                           ed_beam3DGridPointsStatistical, &
                                           ed_beam3DGridSetupStatistical,  &
                                           ed_beam3DGridWeightStatistical, &
                                           ed_computeBeamPower
                                           
  use ed_beamCrossSectionFunctions, ONLY : ed_beamCrossSectionWeight
  
  use EnergyDeposition_data,        ONLY : ed_beams,              &
                                           ed_cellWallThickness,  &
                                           ed_domainErrorMarginX, &
                                           ed_domainErrorMarginY, &
                                           ed_domainErrorMarginZ, &
                                           ed_electronCharge,     &
                                           ed_electronMass,       &
                                           ed_energyInTimestep,   &
                                           ed_globalMe,           &
                                           ed_maxRayCount,        &
                                           ed_numberOfBeams,      &
                                           ed_rayCount,           &
                                           ed_rays,               &
                                           ed_rayZeroPower,       &
                                           ed_speedOfLight,       &
                                           ed_xminDomain,         &
                                           ed_xmaxDomain,         &
                                           ed_yminDomain,         &
                                           ed_ymaxDomain,         &
                                           ed_zminDomain,         &
                                           ed_zmaxDomain
  
  implicit none

#include "EnergyDeposition.h"
#include "constants.h"

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation

  character (len = BEAM_STRING_LENGTH) :: functionType
  character (len = BEAM_STRING_LENGTH) :: gridType

  logical :: beamActive
  logical :: blockFaceMinX, blockFaceMaxX
  logical :: blockFaceMinY, blockFaceMaxY
  logical :: blockFaceMinZ, blockFaceMaxZ
  logical :: blockReflectMinX, blockReflectMaxX
  logical :: blockReflectMinY, blockReflectMaxY
  logical :: blockReflectMinZ, blockReflectMaxZ
  logical :: ignoreBoundary
  logical :: inBlock
  logical :: moreGridPoints
  logical :: onDomainFace
  logical :: parallel2XYplane
  logical :: parallel2XZplane
  logical :: parallel2YZplane
  logical :: proceed
  logical :: rayNotOnBoundary, rayReflects
  logical :: seedIncrement, seedInitialize
  logical :: startGrid
  logical :: xmaxLimit, ymaxLimit, zmaxLimit
  logical :: xminLimit, yminLimit, zminLimit

  integer :: beam
  integer :: blockID
  integer :: block
  integer :: maxGridPoints
  integer :: n
  integer :: nGridPoints
  integer :: nRaysBeam
  integer :: nTics1stDim, nTics2ndDim
  integer :: pulseNumber
  integer :: seed, seedMaximum, seedStepping
  integer :: usedGridPoints

  real    :: beamAspectRatio
  real    :: beamCriticalDensity
  real    :: beamPower
  real    :: cellWallThicknessHalf
  real    :: compDomainXmin, compDomainXmax
  real    :: compDomainYmin, compDomainYmax
  real    :: compDomainZmin, compDomainZmax
  real    :: criticalDensityFactor
  real    :: delta1stDim, delta2ndDim
  real    :: distToFaceMinX, distToFaceMinY, distToFaceMinZ
  real    :: distToFaceMaxX, distToFaceMaxY, distToFaceMaxZ
  real    :: eTx, eTy, eLx, eLy
  real    :: firstTic1stDim, firstTic2ndDim
  real    :: frequency
  real    :: gaussianCenterMajor, gaussianCenterMinor
  real    :: gaussianExponent
  real    :: gaussianRadiusMajor, gaussianRadiusMinor
  real    :: gridWeight, gridWeightInv
  real    :: initialRaySpeed
  real    :: lensX, lensY, lensZ
  real    :: Lx, Ly, Lz
  real    :: magnifyT2L
  real    :: powerFraction
  real    :: Px, Py, Pz
  real    :: rayPower, raySpeed, rayWeight
  real    :: rayX, rayY, rayZ
  real    :: RsizeInv
  real    :: Rx, Ry, Rz
  real    :: RxInv, RyInv, RzInv
  real    :: targetSemiAxisMajor, targetSemiAxisMinor
  real    :: targetX, targetY, targetZ
  real    :: Tx, Ty, Tz
  real    :: u1x, u1y, u1z
  real    :: u2x, u2y, u2z
  real    :: uRx, uRy, uRz
  real    :: userDomainXmin, userDomainXmax
  real    :: userDomainYmin, userDomainYmax
  real    :: userDomainZmin, userDomainZmax
  real    :: velX, velY, velZ
  real    :: w, wtest, wNotValid
  real    :: x, y
  real    :: xminBlock, yminBlock, zminBlock
  real    :: xmaxBlock, ymaxBlock, zmaxBlock

  integer :: faces  (LOW:HIGH,1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)

  real,    allocatable :: xGrid         (:)
  real,    allocatable :: yGrid         (:)
  integer, allocatable :: blockBoundary (:,:)
  real,    allocatable :: blockBndBox   (:,:)
!
!
!     ...Allocate arrays for holding the beam grid points.
!
!
  maxGridPoints = BEAM_GRID_ARRAYSIZE

  allocate (xGrid (1:maxGridPoints))
  allocate (yGrid (1:maxGridPoints))
!
!
!     ...In order to avoid repeated calls to the grid get bounding box and boundary
!        condition function inside the innermost loop for each ray, we determine the
!        bounding box and boundary conditions for each block beforehand and store
!        this info away for future reference.
!
!
  allocate (blockBndBox   (1:6,1:blockCount))
  allocate (blockBoundary (1:6,1:blockCount))

  do block = 1, blockCount

     blockID = blockList (block)

     call Grid_getBlkBC       (blockID, faces)
     call Grid_getBlkBoundBox (blockID, bndBox)

     blockBndBox   (1,block) = bndBox (LOW ,IAXIS)
     blockBndBox   (2,block) = bndBox (HIGH,IAXIS)
     blockBndBox   (3,block) = bndBox (LOW ,JAXIS)
     blockBndBox   (4,block) = bndBox (HIGH,JAXIS)
     blockBndBox   (5,block) = bndBox (LOW ,KAXIS)
     blockBndBox   (6,block) = bndBox (HIGH,KAXIS)

     blockBoundary (1,block) = faces  (LOW ,IAXIS)
     blockBoundary (2,block) = faces  (HIGH,IAXIS)
     blockBoundary (3,block) = faces  (LOW ,JAXIS)
     blockBoundary (4,block) = faces  (HIGH,JAXIS)
     blockBoundary (5,block) = faces  (LOW ,KAXIS)
     blockBoundary (6,block) = faces  (HIGH,KAXIS)

  end do
!
!
!     ...Set the computational domain boundaries. These will be used to decide, whether
!        rays hit the domain boundaries. The computational domain is slightly larger
!        in each dimensional direction to allow for computational rounding errors. Usage
!        of the user defined domain only results in potential loss of rays due to boundary
!        misses.
!
!
  userDomainXmin = ed_xminDomain
  userDomainXmax = ed_xmaxDomain
  userDomainYmin = ed_yminDomain
  userDomainYmax = ed_ymaxDomain
  userDomainZmin = ed_zminDomain
  userDomainZmax = ed_zmaxDomain

  compDomainXmin = ed_xminDomain - ed_domainErrorMarginX
  compDomainXmax = ed_xmaxDomain + ed_domainErrorMarginX
  compDomainYmin = ed_yminDomain - ed_domainErrorMarginY
  compDomainYmax = ed_ymaxDomain + ed_domainErrorMarginY
  compDomainZmin = ed_zminDomain - ed_domainErrorMarginZ
  compDomainZmax = ed_zmaxDomain + ed_domainErrorMarginZ
!
!
!     ...Set some extra needed data.
!
!
  criticalDensityFactor = ed_electronMass * PI / (ed_electronCharge * ed_electronCharge)
  cellWallThicknessHalf = 0.5 * ed_cellWallThickness
!
!
!     ...Outer loop over all (still active) beams. A beam is considered active, if its
!        power is such that it will produce on average rays that are considered to have
!        power above the ray zero power threshold.
!
!
  do beam = 1, ed_numberOfBeams

     nRaysBeam   = ed_beams (beam) % numberOfRays
     pulseNumber = ed_beams (beam) % pulseNumber

     call ed_computeBeamPower (timeStep, timeSimulation, pulseNumber,   beamPower)

     beamActive = (beamPower > real (nRaysBeam) * ed_rayZeroPower)

     if (beamActive) then

         delta1stDim         = ed_beams (beam) % gridDelta1stDim
         delta2ndDim         = ed_beams (beam) % gridDelta2ndDim
         firstTic1stDim      = ed_beams (beam) % gridFirstTic1stDim
         firstTic2ndDim      = ed_beams (beam) % gridFirstTic2ndDim
         frequency           = ed_beams (beam) % frequency
         functionType        = ed_beams (beam) % crossSectionFunctionType
         gaussianCenterMajor = ed_beams (beam) % gaussianCenterMajor
         gaussianCenterMinor = ed_beams (beam) % gaussianCenterMinor
         gaussianExponent    = ed_beams (beam) % gaussianExponent
         gaussianRadiusMajor = ed_beams (beam) % gaussianRadiusMajor
         gaussianRadiusMinor = ed_beams (beam) % gaussianRadiusMinor
         gridType            = ed_beams (beam) % gridType
         ignoreBoundary      = ed_beams (beam) % ignoreBoundaryCondition
         initialRaySpeed     = ed_beams (beam) % initialRaySpeed
         lensX               = ed_beams (beam) % lensX
         lensY               = ed_beams (beam) % lensY
         lensZ               = ed_beams (beam) % lensZ
         magnifyT2L          = ed_beams (beam) % target2LensMagnification
         nTics1stDim         = ed_beams (beam) % gridnTics1stDim
         nTics2ndDim         = ed_beams (beam) % gridnTics2ndDim
         seed                = ed_beams (beam) % gridSeed
         seedMaximum         = ed_beams (beam) % gridSeedMaximum
         seedStepping        = ed_beams (beam) % gridSeedStepping
         targetSemiAxisMajor = ed_beams (beam) % targetSemiAxisMajor
         targetSemiAxisMinor = ed_beams (beam) % targetSemiAxisMinor
         targetX             = ed_beams (beam) % targetX
         targetY             = ed_beams (beam) % targetY
         targetZ             = ed_beams (beam) % targetZ
         u1x                 = ed_beams (beam) % semiAxisUnitMajorX
         u1y                 = ed_beams (beam) % semiAxisUnitMajorY
         u1z                 = ed_beams (beam) % semiAxisUnitMajorZ
         u2x                 = ed_beams (beam) % semiAxisUnitMinorX
         u2y                 = ed_beams (beam) % semiAxisUnitMinorY
         u2z                 = ed_beams (beam) % semiAxisUnitMinorZ

         beamAspectRatio     = targetSemiAxisMajor / targetSemiAxisMinor
         beamCriticalDensity = criticalDensityFactor * frequency * frequency
!
!
!     ...Get the proper grid weight. In case of statistical grids, the grid weight has
!        to be recalculated every time step. For regular grids, we use the value evaluated
!        during setup of the beams. The statistical grids need a setup call every time step,
!        in order to set a different seed value for the random number generator. This
!        new seed value needs to be stored into the beams array.
!
!
         if (gridType == 'statistical2D') then

             seedInitialize = .false.
             seedIncrement  = .true.

             call ed_beam3DGridSetupStatistical  (seedInitialize,                 &
                                                  seedIncrement,                  &
                                                  seedMaximum,                    &
                                                  seedStepping,                   &
                                                                             seed )

             call ed_beam3DGridWeightStatistical (targetSemiAxisMajor,            &
                                                  targetSemiAxisMinor,            &
                                                  functionType,                   &
                                                  gaussianExponent,               &
                                                  gaussianRadiusMajor,            &
                                                  gaussianRadiusMinor,            &
                                                  gaussianCenterMajor,            &
                                                  gaussianCenterMinor,            &
                                                  seed,                           &
                                                  nRaysBeam,                      &
                                                                       gridWeight )

             ed_beams (beam) % gridSeed   = seed
             ed_beams (beam) % gridWeight = gridWeight
         else
             gridWeight = ed_beams (beam) % gridWeight
         end if

         gridWeightInv = 1.0 / gridWeight
!
!
!     ...Create all rays for the current beam. This is done by looping over all possible
!        grid points in both the target and lens areas.
!
!        For each ray we obtain a global x,y,z position for both the lens and target area.
!        The ray can be represented by a vector 'R' with tail point 'L' on the lens and head
!        point 'T' at the target. It therefore lays on a line in 3D, given by the following
!        parametric line vector equation:
!
!                                      P = L + w * R
!                                        = L + w * (T - L)
!
!        where 'w' is a real parameter that can attain all possible values. In particular we
!        must have:
!
!                              w >= 0    -->   ray moves from L to T
!
!        Having the ray line equation, we now must see if and where the domain boundaries
!        are crossed. For each domain boundary we have very simple plane equations:
!
!                     domain yz-faces  -->  plane equations: x = xminDomain and x = xmaxDomain
!                     domain xz-faces  -->  plane equations: y = yminDomain and y = ymaxDomain
!                     domain xy-faces  -->  plane equations: z = zminDomain and z = zmaxDomain
!
!        Lets assume we have a plane equation x = a. Then we find the 'w' value, such that Px = a
!        is on the plane. We obtain from the ray line equation:
!
!                                    w = (a - Lx) / (Tx - Lx)
!
!        From this we get the corresponding Py and Pz coordinates:
!
!                                   Py = Ly + w * (Ty - Ly)
!                                   Pz = Lz + w * (Tz - Lz)
!
!        Significance of the 'w' values:
!
!                         i)     w < 0  -->  ray moves away from the domain         (exclude)
!                        ii)     w > 1  -->  target point is outside of the domain  (exclude)
!                       iii)     w = 0  -->  lens area touches domain boundary      (exclude)
!                        iv)     w = 1  -->  target area touches domain boundary    (ok,can happen)
!                         v) 0 < w < 1  -->  proper lens and target positions       (ok)
!
!        If the 'w' value is ok, we have to check, if the point (Py,Pz) is actually contained on the
!        domain yz-face. If yes, add this 'w' to the allowed w-list. To find the relevant 'w' we
!        need to pick the smallest 'w' from the allowed w-list. If the list turns out to be
!        empty, then we know that the ray is moving away from the domain. This should never
!        happen as the target area is supposed to be entirely contained in the domain.
!
!
!
         startGrid      = .true.
         moreGridPoints = .true.
         usedGridPoints = 0

         do while (moreGridPoints)

            if (gridType == 'delta2D') then

                call ed_beam3DGridPointsDelta       (targetSemiAxisMajor,            &
                                                     targetSemiAxisMinor,            &
                                                     nTics1stDim,                    &
                                                     nTics2ndDim,                    &
                                                     delta1stDim,                    &
                                                     delta2ndDim,                    &
                                                     firstTic1stDim,                 &
                                                     firstTic2ndDim,                 &
                                                     startGrid,                      &
                                                     maxGridPoints,                  &
                                                                     moreGridPoints, &
                                                                     nGridPoints,    &
                                                                     xGrid,          &
                                                                     yGrid           )
            else if (gridType == 'square2D') then

                call ed_beam3DGridPointsSquare      (targetSemiAxisMajor,            &
                                                     targetSemiAxisMinor,            &
                                                     nTics1stDim,                    &
                                                     nTics2ndDim,                    &
                                                     delta1stDim,                    &
                                                     startGrid,                      &
                                                     maxGridPoints,                  &
                                                                     moreGridPoints, &
                                                                     nGridPoints,    &
                                                                     xGrid,          &
                                                                     yGrid           )
            else if (gridType == 'radial2D') then

                call ed_beam3DGridPointsRadial      (targetSemiAxisMajor,            &
                                                     targetSemiAxisMinor,            &
                                                     nTics1stDim,                    &
                                                     nTics2ndDim,                    &
                                                     delta1stDim,                    &
                                                     delta2ndDim,                    &
                                                     firstTic1stDim,                 &
                                                     firstTic2ndDim,                 &
                                                     startGrid,                      &
                                                     maxGridPoints,                  &
                                                                     moreGridPoints, &
                                                                     nGridPoints,    &
                                                                     xGrid,          &
                                                                     yGrid           )
            else if (gridType == 'rectangular2D') then

                call ed_beam3DGridPointsRecBeam     (nTics1stDim,                    &
                                                     nTics2ndDim,                    &
                                                     delta1stDim,                    &
                                                     delta2ndDim,                    &
                                                     startGrid,                      &
                                                     maxGridPoints,                  &
                                                                     moreGridPoints, &
                                                                     nGridPoints,    &
                                                                     xGrid,          &
                                                                     yGrid           )
            else if (gridType == 'statistical2D') then

                call ed_beam3DGridPointsStatistical (targetSemiAxisMajor,            &
                                                     targetSemiAxisMinor,            &
                                                     seed,                           &
                                                     nRaysBeam,                      &
                                                     startGrid,                      &
                                                     maxGridPoints,                  &
                                                                     moreGridPoints, &
                                                                     nGridPoints,    &
                                                                     xGrid,          &
                                                                     yGrid           )
            else
                call Driver_abortFlash ('[ed_createRays3DRec] ERROR: unknown 3D beam grid type')
            end if

            do n = 1, nGridPoints

               eTx = xGrid (n)                                     ! local x-coordinate in target ellipse
               eTy = yGrid (n)                                     ! local y-coordinate in target ellipse
               eLx = eTx * magnifyT2L                              ! local x-coordinate in lens ellipse
               eLy = eTy * magnifyT2L                              ! local y-coordinate in lens ellipse

               Lx = eLx * u1x  +  eLy * u2x  +  lensX              ! global lens ray x-coordinate
               Ly = eLx * u1y  +  eLy * u2y  +  lensY              ! global lens ray y-coordinate
               Lz = eLx * u1z  +  eLy * u2z  +  lensZ              ! global lens ray z-coordinate
               Tx = eTx * u1x  +  eTy * u2x  +  targetX            ! global target ray x-coordinate
               Ty = eTx * u1y  +  eTy * u2y  +  targetY            ! global target ray y-coordinate
               Tz = eTx * u1z  +  eTy * u2z  +  targetZ            ! global target ray z-coordinate

               Rx = Tx - Lx                                        ! ray vector x-coordinate
               Ry = Ty - Ly                                        ! ray vector y-coordinate
               Rz = Tz - Lz                                        ! ray vector z-coordinate

               parallel2YZplane = (Rx == 0.0)
               parallel2XZplane = (Ry == 0.0)
               parallel2XYplane = (Rz == 0.0)
!
!
!     ...Start the search for the smallest 'w' such that: 0 < w <= 1. The initial condition on
!        the 'w' is such that it corresponds to a nonvalid positive number, which can be checked
!        against after the search is completed.
!
!        Check first, if ray hits any of the two boundary yz-faces of the domain. Note, that
!        while the check is being performed using the computational domain, the rays will be
!        placed exactly on the user defined domain.
!
!
               wNotValid = 2.0
               w = wNotValid

               if (.not. parallel2YZplane) then

                   RxInv   = 1.0 / Rx
                   wtest   = (userDomainXmin - Lx) * RxInv
                   proceed = (wtest > 0.0) .and. (wtest <= 1.0) .and. (wtest <= w)

                   if (proceed) then
                       Py = Ly + wtest * Ry
                       Pz = Lz + wtest * Rz
                       onDomainFace =      (Py >= compDomainYmin) .and. (Py <= compDomainYmax) &
                                     .and. (Pz >= compDomainZmin) .and. (Pz <= compDomainZmax)
                       if (onDomainFace) then
                           w = wtest
                           rayX = userDomainXmin
                           rayY = min (  Py, userDomainYmax)          !
                           rayY = max (rayY, userDomainYmin)          ! force rays on user defined domain
                           rayZ = min (  Pz, userDomainZmax)          !
                           rayZ = max (rayZ, userDomainZmin)          !
                       end if
                   end if

                   wtest   = (userDomainXmax - Lx) * RxInv
                   proceed = (wtest > 0.0) .and. (wtest <= 1.0) .and. (wtest <= w)

                   if (proceed) then
                       Py = Ly + wtest * Ry
                       Pz = Lz + wtest * Rz
                       onDomainFace =      (Py >= compDomainYmin) .and. (Py <= compDomainYmax) &
                                     .and. (Pz >= compDomainZmin) .and. (Pz <= compDomainZmax)
                       if (onDomainFace) then
                           w = wtest
                           rayX = userDomainXmax
                           rayY = min (  Py, userDomainYmax)
                           rayY = max (rayY, userDomainYmin)
                           rayZ = min (  Pz, userDomainZmax)
                           rayZ = max (rayZ, userDomainZmin)
                       end if
                   end if

               end if
!
!
!     ...Next check, if ray hits any of the two boundary xz-faces of the domain.
!
!
               if (.not. parallel2XZplane) then

                   RyInv   = 1.0 / Ry
                   wtest   = (userDomainYmin - Ly) * RyInv
                   proceed = (wtest > 0.0) .and. (wtest <= 1.0) .and. (wtest <= w)

                   if (proceed) then
                       Px = Lx + wtest * Rx
                       Pz = Lz + wtest * Rz
                       onDomainFace =      (Px >= compDomainXmin) .and. (Px <= compDomainXmax) &
                                     .and. (Pz >= compDomainZmin) .and. (Pz <= compDomainZmax)
                       if (onDomainFace) then
                           w = wtest
                           rayY = userDomainYmin
                           rayX = min (  Px, userDomainXmax)
                           rayX = max (rayX, userDomainXmin)
                           rayZ = min (  Pz, userDomainZmax)
                           rayZ = max (rayZ, userDomainZmin)
                       end if
                   end if

                   wtest   = (userDomainYmax - Ly) * RyInv
                   proceed = (wtest > 0.0) .and. (wtest <= 1.0) .and. (wtest <= w)

                   if (proceed) then
                       Px = Lx + wtest * Rx
                       Pz = Lz + wtest * Rz
                       onDomainFace =      (Px >= compDomainXmin) .and. (Px <= compDomainXmax) &
                                     .and. (Pz >= compDomainZmin) .and. (Pz <= compDomainZmax)
                       if (onDomainFace) then
                           w = wtest
                           rayY = userDomainYmax
                           rayX = min (  Px, userDomainXmax)
                           rayX = max (rayX, userDomainXmin)
                           rayZ = min (  Pz, userDomainZmax)
                           rayZ = max (rayZ, userDomainZmin)
                       end if
                   end if

               end if
!
!
!     ...Finally, check, if ray hits any of the two boundary xy-faces of the domain.
!
!
               if (.not. parallel2XYplane) then

                   RzInv   = 1.0 / Rz
                   wtest   = (userDomainZmin - Lz) * RzInv
                   proceed = (wtest > 0.0) .and. (wtest <= 1.0) .and. (wtest <= w)

                   if (proceed) then
                       Px = Lx + wtest * Rx
                       Py = Ly + wtest * Ry
                       onDomainFace =      (Px >= compDomainXmin) .and. (Px <= compDomainXmax) &
                                     .and. (Py >= compDomainYmin) .and. (Py <= compDomainYmax)
                       if (onDomainFace) then
                           w = wtest
                           rayZ = userDomainZmin
                           rayX = min (  Px, userDomainXmax)
                           rayX = max (rayX, userDomainXmin)
                           rayY = min (  Py, userDomainYmax)
                           rayY = max (rayY, userDomainYmin)
                       end if
                   end if

                   wtest   = (userDomainZmax - Lz) * RzInv
                   proceed = (wtest > 0.0) .and. (wtest <= 1.0) .and. (wtest <= w)

                   if (proceed) then
                       Px = Lx + wtest * Rx
                       Py = Ly + wtest * Ry
                       onDomainFace =      (Px >= compDomainXmin) .and. (Px <= compDomainXmax) &
                                     .and. (Py >= compDomainYmin) .and. (Py <= compDomainYmax)
                       if (onDomainFace) then
                           w = wtest
                           rayZ = userDomainZmax
                           rayX = min (  Px, userDomainXmax)
                           rayX = max (rayX, userDomainXmin)
                           rayY = min (  Py, userDomainYmax)
                           rayY = max (rayY, userDomainYmin)
                       end if
                   end if

               end if
!
!
!     ...Catch the situation in which the ray did not make it. This should never happen,
!        because the entire target area is inside the domain. Note, that if it passes this
!        test, then we know for sure (from the code structure above) that the ray's position
!        is on one of the boundaries and hence on a block.
!
!
               rayNotOnBoundary = (w == wNotValid)

               if (rayNotOnBoundary) then
                   call Driver_abortFlash("ed_createRays3DRec: ray does not hit domain boundary!")
               end if
!
!
!     ...loop over all blocks and see, if ray is contained in one of it. As soon as that block
!        is found, exit the block loop, since each ray can be assigned to only one block.
!        A ray belongs to a block, if its x,y,z coordinate is such that:
!
!                        x,y,z block lower limit  <=  ray x,y,z  <  x,y,z block upper limit
!
!        except when any of the upper limits of the block coincide with the domain boundaries.
!        In that case the less '<' sign on the right must be replaced by a <= sign, otherwise
!        rays will dissapear and not accounted for, resulting in beam power loss.
!
!
               do block = 1, blockCount

                  xminBlock = blockBndBox (1,block)
                  xmaxBlock = blockBndBox (2,block)
                  yminBlock = blockBndBox (3,block)
                  ymaxBlock = blockBndBox (4,block)
                  zminBlock = blockBndBox (5,block)
                  zmaxBlock = blockBndBox (6,block)

                  xminLimit = (rayX >= xminBlock)
                  yminLimit = (rayY >= yminBlock)
                  zminLimit = (rayZ >= zminBlock)
                  xmaxLimit = (rayX <  xmaxBlock) .or. ((rayX == xmaxBlock) .and. (xmaxBlock == userDomainXmax))
                  ymaxLimit = (rayY <  ymaxBlock) .or. ((rayY == ymaxBlock) .and. (ymaxBlock == userDomainYmax))
                  zmaxLimit = (rayZ <  zmaxBlock) .or. ((rayZ == zmaxBlock) .and. (zmaxBlock == userDomainZmax))

                  inBlock  =      xminLimit &
                            .and. xmaxLimit &
                            .and. yminLimit &
                            .and. ymaxLimit &
                            .and. zminLimit &
                            .and. zmaxLimit

                  if (inBlock) then

                      distToFaceMinX = abs (xminBlock - rayX)
                      distToFaceMaxX = abs (xmaxBlock - rayX)
                      distToFaceMinY = abs (yminBlock - rayY)
                      distToFaceMaxY = abs (ymaxBlock - rayY)
                      distToFaceMinZ = abs (zminBlock - rayZ)
                      distToFaceMaxZ = abs (zmaxBlock - rayZ)

                      blockFaceMinX = (distToFaceMinX < cellWallThicknessHalf)
                      blockFaceMaxX = (distToFaceMaxX < cellWallThicknessHalf)
                      blockFaceMinY = (distToFaceMinY < cellWallThicknessHalf)
                      blockFaceMaxY = (distToFaceMaxY < cellWallThicknessHalf)
                      blockFaceMinZ = (distToFaceMinZ < cellWallThicknessHalf)
                      blockFaceMaxZ = (distToFaceMaxZ < cellWallThicknessHalf)

                      blockReflectMinX = (blockBoundary (1,block) == REFLECTING)
                      blockReflectMaxX = (blockBoundary (2,block) == REFLECTING)
                      blockReflectMinY = (blockBoundary (3,block) == REFLECTING)
                      blockReflectMaxY = (blockBoundary (4,block) == REFLECTING)
                      blockReflectMinZ = (blockBoundary (5,block) == REFLECTING)
                      blockReflectMaxZ = (blockBoundary (6,block) == REFLECTING)

                      rayReflects =     (blockFaceMinX .and. blockReflectMinX) &
                                   .or. (blockFaceMaxX .and. blockReflectMaxX) &
                                   .or. (blockFaceMinY .and. blockReflectMinY) &
                                   .or. (blockFaceMaxY .and. blockReflectMaxY) &
                                   .or. (blockFaceMinZ .and. blockReflectMinZ) &
                                   .or. (blockFaceMaxZ .and. blockReflectMaxZ)

                      if (ignoreBoundary .or. .not.rayReflects) then

                          ed_rayCount = ed_rayCount + 1              ! this count is per processor, not block

                          if (ed_rayCount > ed_maxRayCount) then
                              call Driver_abortFlash ("ed_createRays3DRec: Not enough storage for rays array")
                          end if

                          blockID = blockList (block)

                          RsizeInv = 1.0 / sqrt (Rx * Rx + Ry * Ry + Rz * Rz)

                          uRx = Rx * RsizeInv                        ! unit vector along ray vector
                          uRy = Ry * RsizeInv
                          uRz = Rz * RsizeInv

                          rayWeight = ed_beamCrossSectionWeight (functionType,                   &
                                                                 x        = eTx,                 &
                                                                 y        = eTy,                 &
                                                                 Cx       = gaussianCenterMajor, &
                                                                 Cy       = gaussianCenterMinor, &
                                                                 Rx       = gaussianRadiusMajor, &
                                                                 Ry       = gaussianRadiusMinor, &
                                                                 Exponent = gaussianExponent     )

                          if (gridType == 'radial2D') then
                              x = eTx
                              y = eTy * beamAspectRatio
                              rayWeight = rayWeight * sqrt (x * x + y * y)
                          end if

                          powerFraction = rayWeight * gridWeightInv

                          rayPower = beamPower * powerFraction
                          raySpeed = initialRaySpeed * ed_speedOfLight

                          velX = uRx * raySpeed
                          velY = uRy * raySpeed
                          velZ = uRz * raySpeed
!
!
!     ...We now need to make sure, that the ray is placed inside the block and not within the block
!        wall. If this is not done, then initial cell index determination can suffer from roundoff
!        errors with placement of the ray in a wrong cell ouside the current block. Note, that even
!        after pushing the ray inside the block is done, there is still the possibility, that the ray
!        sits exactly on a cell face boundary, but this will happen inside the block, so there is no
!        danger of placing the ray initially in the wrong block. The appropriate (initial) nudging for
!        the cell will be done in the corresponding ray tracing routine, as here we do not have the
!        individual cell information.
!
!
                          if (blockFaceMinX) then
                              rayX = xminBlock + cellWallThicknessHalf
                          else if (blockFaceMaxX) then
                              rayX = xmaxBlock - cellWallThicknessHalf
                          end if

                          if (blockFaceMinY) then
                              rayY = yminBlock + cellWallThicknessHalf
                          else if (blockFaceMaxY) then
                              rayY = ymaxBlock - cellWallThicknessHalf
                          end if

                          if (blockFaceMinZ) then
                              rayZ = zminBlock + cellWallThicknessHalf
                          else if (blockFaceMaxZ) then
                              rayZ = zmaxBlock - cellWallThicknessHalf
                          end if

                          ed_rays (RAY_BLCK,ed_rayCount) = real (blockID)
                          ed_rays (RAY_PROC,ed_rayCount) = real (ed_globalMe)
                          ed_rays (RAY_POSX,ed_rayCount) = rayX
                          ed_rays (RAY_POSY,ed_rayCount) = rayY
                          ed_rays (RAY_POSZ,ed_rayCount) = rayZ
                          ed_rays (RAY_VELX,ed_rayCount) = velX
                          ed_rays (RAY_VELY,ed_rayCount) = velY
                          ed_rays (RAY_VELZ,ed_rayCount) = velZ
                          ed_rays (RAY_POWR,ed_rayCount) = rayPower
                          ed_rays (RAY_DENC,ed_rayCount) = beamCriticalDensity

                          ed_energyInTimestep = ed_energyInTimestep + rayPower * timeStep

                          exit

                      end if                                      ! ray reflection condition
                  end if                                          ! ray in block
               end do                                             ! block loop
            end do                                                ! individual grid points loop

            usedGridPoints = usedGridPoints + nGridPoints

         end do                                                   ! array of grid points loop

         if (usedGridPoints /= nRaysBeam) then
             call Driver_abortFlash ("ed_createRays3DRec: # of rays / beam grid points mismatch !")
         end if

     end if                                                       ! active beam
  end do                                                          ! beam loop
!
!
!     ...Deallocate all allocated arrays.
!
!
  deallocate (xGrid)
  deallocate (yGrid)
  deallocate (blockBndBox)
  deallocate (blockBoundary)
!
!
!     ...Ready!
!
!
  return
end subroutine ed_createRays3DRec
