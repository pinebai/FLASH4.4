!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRaysCreate/ed_createRays2DRec
!!
!! NAME
!!
!!  ed_createRays2DRec
!!
!! SYNOPSIS
!!
!!  call ed_createRays2DRec (integer, intent (in) :: blockCount, 
!!                           integer, intent (in) :: blockList (:), 
!!                           real,    intent (in) :: timeStep,   
!!                           real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates rays and places them in their initial blocks for those geometries consisting
!!  formally of 2D rectangular grids (cartesian + cylindrical). On exit, all rays hitting
!!  the domain boundary have been generated for the current processor. Their block ID's are
!!  not ordered as the outer loop is over all beams.
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

subroutine ed_createRays2DRec (blockCount, blockList, timeStep, timeSimulation)

  use Driver_interface,             ONLY : Driver_abortFlash
 
  use Grid_interface,               ONLY : Grid_getBlkBC,      &
                                           Grid_getBlkBoundBox

  use ed_interface,                 ONLY : ed_beam2DGridPointsRegular,     &
                                           ed_beam2DGridPointsStatistical, &
                                           ed_beam2DGridSetupStatistical,  &
                                           ed_beam2DGridWeightStatistical, &
                                           ed_computeBeamPower

  use ed_beamCrossSectionFunctions, ONLY : ed_beamCrossSectionWeight
  
  use EnergyDeposition_data,        ONLY : ed_beams,              &
                                           ed_cellWallThickness,  &
                                           ed_domainErrorMarginX, &
                                           ed_domainErrorMarginY, &
                                           ed_electronCharge,     &
                                           ed_electronMass,       &
                                           ed_energyInTimestep,   &
                                           ed_globalMe,           &
                                           ed_gridGeometry,       &
                                           ed_maxRayCount,        &
                                           ed_numberOfBeams,      &
                                           ed_rayCount,           &
                                           ed_rays,               &
                                           ed_rayZeroPower,       &
                                           ed_speedOfLight,       &
                                           ed_xminDomain,         &
                                           ed_xmaxDomain,         &
                                           ed_yminDomain,         &
                                           ed_ymaxDomain
  
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
  logical :: blockReflectMinX, blockReflectMaxX
  logical :: blockReflectMinY, blockReflectMaxY
  logical :: grid2Dcylindrical
  logical :: ignoreBoundary
  logical :: inBlock
  logical :: moreGridPoints
  logical :: onDomainLine
  logical :: parallel2Xaxis
  logical :: parallel2Yaxis
  logical :: proceed
  logical :: rayNotOnBoundary, rayReflects
  logical :: seedIncrement, seedInitialize
  logical :: startGrid
  logical :: xmaxLimit, ymaxLimit
  logical :: xminLimit, yminLimit

  integer :: beam
  integer :: blockID
  integer :: block
  integer :: maxGridPoints
  integer :: n
  integer :: nGridPoints
  integer :: nRaysBeam
  integer :: nTics
  integer :: pulseNumber
  integer :: seed, seedMaximum, seedStepping
  integer :: usedGridPoints

  real    :: beamCriticalDensity
  real    :: beamPower
  real    :: cellWallThicknessHalf
  real    :: compDomainXmin, compDomainXmax
  real    :: compDomainYmin, compDomainYmax
  real    :: criticalDensityFactor
  real    :: delta
  real    :: distToFaceMinX, distToFaceMinY
  real    :: distToFaceMaxX, distToFaceMaxY
  real    :: firstTic
  real    :: frequency
  real    :: gaussianCenter
  real    :: gaussianExponent
  real    :: gaussianRadius
  real    :: gridWeight, gridWeightInv
  real    :: initialRaySpeed
  real    :: lensX, lensY
  real    :: lineTx, lineLx
  real    :: Lx, Ly
  real    :: magnifyT2L
  real    :: powerFraction
  real    :: Px, Py
  real    :: rayPower, raySpeed
  real    :: rayWeight
  real    :: rayX, rayY
  real    :: RsizeInv
  real    :: Rx, Ry
  real    :: RxInv, RyInv
  real    :: targetSemiAxis
  real    :: targetX, targetY
  real    :: Tx, Ty
  real    :: uRx, uRy
  real    :: userDomainXmin, userDomainXmax
  real    :: userDomainYmin, userDomainYmax
  real    :: ux, uy
  real    :: velX, velY
  real    :: w, wtest, wNotValid
  real    :: xminBlock, yminBlock
  real    :: xmaxBlock, ymaxBlock

  integer :: faces  (LOW:HIGH,1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)

  real,    allocatable :: xGrid         (:)
  integer, allocatable :: blockBoundary (:,:)
  real,    allocatable :: blockBndBox   (:,:)
!
!
!     ...Allocate arrays for holding the beam grid points.
!
!
  maxGridPoints = BEAM_GRID_ARRAYSIZE

  allocate (xGrid (1:maxGridPoints))
!
!
!     ...In order to avoid repeated calls to the grid get bounding box and boundary
!        condition function inside the innermost loop for each ray, we determine the
!        bounding box and boundary conditions for each block beforehand and store
!        this info away for future reference.
!
!
  allocate (blockBndBox   (1:4,1:blockCount))
  allocate (blockBoundary (1:4,1:blockCount))

  do block = 1, blockCount

     blockID = blockList (block)

     call Grid_getBlkBC       (blockID, faces)
     call Grid_getBlkBoundBox (blockID, bndBox)

     blockBndBox   (1,block) = bndBox (LOW ,IAXIS)
     blockBndBox   (2,block) = bndBox (HIGH,IAXIS)
     blockBndBox   (3,block) = bndBox (LOW ,JAXIS)
     blockBndBox   (4,block) = bndBox (HIGH,JAXIS)

     blockBoundary (1,block) = faces  (LOW ,IAXIS)
     blockBoundary (2,block) = faces  (HIGH,IAXIS)
     blockBoundary (3,block) = faces  (LOW ,JAXIS)
     blockBoundary (4,block) = faces  (HIGH,JAXIS)

  end do
!
!
!     ...Set some extra needed data.
!
!
  grid2Dcylindrical     = (ed_gridGeometry == GRID_2DCYLINDRICAL)
  criticalDensityFactor =  ed_electronMass * PI / (ed_electronCharge * ed_electronCharge)
  cellWallThicknessHalf =  0.5 * ed_cellWallThickness
!
!
!     ...Set the computational domain boundaries. These will be used to decide, whether
!        rays hit the domain boundaries. The computational domain is slightly larger
!        in each dimensional direction to allow for computational rounding errors. Usage
!        of the user defined domain only results in potential loss of rays due to boundary
!        misses.
!
!        In case of 2D cylindrical geometries we expand the computational R-domain
!        (x-coordinate) temporarily to allow for the -R part of the cylindrical domain.
!        This is done to accommodate possible beams defined within the entire cylindrical
!        range. Any ray that hits on the -R portion must then be mirrored to the [0,R]
!        portion of the computational 2D cylindrical domain. This is done by making the
!        following transformations, if rayX is < 0:
!
!                rayX -> - rayX
!                velX -> - velX   (done by inverting the ray vector component Rx -> - Rx)
!
!
  userDomainXmin = ed_xminDomain
  userDomainXmax = ed_xmaxDomain
  userDomainYmin = ed_yminDomain
  userDomainYmax = ed_ymaxDomain

  compDomainXmin = ed_xminDomain - ed_domainErrorMarginX
  compDomainXmax = ed_xmaxDomain + ed_domainErrorMarginX
  compDomainYmin = ed_yminDomain - ed_domainErrorMarginY
  compDomainYmax = ed_ymaxDomain + ed_domainErrorMarginY

  if (grid2Dcylindrical) then
      userDomainXmin = - userDomainXmax         ! expand to [-R,0] range
      compDomainXmin = - compDomainXmax         ! expand to [-R,0] range
  end if
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

         delta            = ed_beams (beam) % gridDelta1stDim
         firstTic         = ed_beams (beam) % gridFirstTic1stDim
         frequency        = ed_beams (beam) % frequency
         functionType     = ed_beams (beam) % crossSectionFunctionType
         gaussianCenter   = ed_beams (beam) % gaussianCenterMajor
         gaussianExponent = ed_beams (beam) % gaussianExponent
         gaussianRadius   = ed_beams (beam) % gaussianRadiusMajor
         gridType         = ed_beams (beam) % gridType
         ignoreBoundary   = ed_beams (beam) % ignoreBoundaryCondition
         initialRaySpeed  = ed_beams (beam) % initialRaySpeed
         lensX            = ed_beams (beam) % lensX
         lensY            = ed_beams (beam) % lensY
         magnifyT2L       = ed_beams (beam) % target2LensMagnification
         nTics            = ed_beams (beam) % gridnTics1stDim
         seed             = ed_beams (beam) % gridSeed
         seedMaximum      = ed_beams (beam) % gridSeedMaximum
         seedStepping     = ed_beams (beam) % gridSeedStepping
         targetSemiAxis   = ed_beams (beam) % targetSemiAxisMajor
         targetX          = ed_beams (beam) % targetX
         targetY          = ed_beams (beam) % targetY
         ux               = ed_beams (beam) % semiAxisUnitMajorX
         uy               = ed_beams (beam) % semiAxisUnitMajorY

         beamCriticalDensity  = criticalDensityFactor * frequency * frequency
!
!
!     ...Get the proper grid weight. In case of statistical grids, the grid weight has
!        to be recalculated every time step. For regular grids, we use the value evaluated
!        during setup of the beams. The statistical grids need a setup call every time step,
!        in order to set a different seed value for the random number generator. This
!        new seed value needs to be stored into the beams array.
!
!
         if (gridType == 'statistical1D') then

             seedInitialize = .false.
             seedIncrement  = .true.

             call ed_beam2DGridSetupStatistical  (seedInitialize,                 &
                                                  seedIncrement,                  &
                                                  seedMaximum,                    &
                                                  seedStepping,                   &
                                                                             seed )

             call ed_beam2DGridWeightStatistical (targetX,                        &
                                                  targetSemiAxis,                 &
                                                  functionType,                   &
                                                  gaussianExponent,               &
                                                  gaussianRadius,                 &
                                                  gaussianCenter,                 &
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
!        grid points in both the target and lens linear areas.
!
!        For each ray we obtain a global x,y position for both the lens and target line.
!        The ray can be represented by a vector 'R' with tail point 'L' on the lens and head
!        point 'T' at the target. It therefore lays on a line in 2D, given by the following
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
!        are crossed. For each domain boundary we have very simple line equations:
!
!                     domain y-lines  -->  line equations: x = xminDomain and x = xmaxDomain
!                     domain x-lines  -->  line equations: y = yminDomain and y = ymaxDomain
!
!        Lets assume we have a line equation x = a. Then we find the 'w' value, such that Px = a
!        is on the plane. We obtain from the ray line equation:
!
!                                    w = (a - Lx) / (Tx - Lx)
!
!        From this we get the corresponding Py coordinate:
!
!                                   Py = Ly + w * (Ty - Ly)
!
!        Significance of the 'w' values:
!
!                         i)     w < 0  -->  ray moves away from the domain         (exclude)
!                        ii)     w > 1  -->  target point is outside of the domain  (exclude)
!                       iii)     w = 0  -->  lens area touches domain boundary      (exclude)
!                        iv)     w = 1  -->  target area touches domain boundary    (ok,can happen)
!                         v) 0 < w < 1  -->  proper lens and target positions       (ok)
!
!        If the 'w' value is ok, we have to check, if the point (Py) is actually contained on the
!        domain y-face. If yes, add this 'w' to the allowed w-list. To find the relevant 'w' we
!        need to pick the smallest 'w' from the allowed w-list. If the list turns out to be
!        empty, then we know that the ray is moving away from the domain. This should never
!        happen as the target area is supposed to be entirely contained in the domain.
!
!
         startGrid      = .true.
         moreGridPoints = .true.
         usedGridPoints = 0

         do while (moreGridPoints)

            if (gridType == 'regular1D') then

                call ed_beam2DGridPointsRegular     (targetSemiAxis,                 &
                                                     nTics,                          &
                                                     delta,                          &
                                                     firstTic,                       &
                                                     startGrid,                      &
                                                     maxGridPoints,                  &
                                                                     moreGridPoints, &
                                                                     nGridPoints,    &
                                                                     xGrid           )
            else if (gridType == 'statistical1D') then

                call ed_beam2DGridPointsStatistical (targetSemiAxis,                 &
                                                     seed,                           &
                                                     nRaysBeam,                      &
                                                     startGrid,                      &
                                                     maxGridPoints,                  &
                                                                     moreGridPoints, &
                                                                     nGridPoints,    &
                                                                     xGrid           )
            else
                call Driver_abortFlash ('[ed_createRays2DRec] ERROR: unknown 2D beam grid type')
            end if

            do n = 1, nGridPoints

               lineTx = xGrid (n)                                 ! local x-coordinate in target line
               lineLx = lineTx * magnifyT2L                       ! local x-coordinate in lens line

               Lx = lineLx * ux  +  lensX                         ! global lens ray x-coordinate
               Ly = lineLx * uy  +  lensY                         ! global lens ray y-coordinate
               Tx = lineTx * ux  +  targetX                       ! global target ray x-coordinate
               Ty = lineTx * uy  +  targetY                       ! global target ray y-coordinate

               Rx = Tx - Lx                                       ! ray vector x-component
               Ry = Ty - Ly                                       ! ray vector y-component

               parallel2Xaxis = (Ry == 0.0)
               parallel2Yaxis = (Rx == 0.0)
!
!
!     ...Start the search for the smallest 'w' such that: 0 < w <= 1. The initial condition on
!        the 'w' is such that it corresponds to a nonvalid positive number, which can be checked
!        against after the search is completed.
!
!        Check first, if ray hits any of the two boundary y-lines of the domain. Note, that
!        while the check is being performed using the computational domain, the rays will be
!        placed exactly on the user defined domain.
!
!
               wNotValid = 2.0
               w = wNotValid

               if (.not. parallel2Yaxis) then

                   RxInv   = 1.0 / Rx
                   wtest   = (userDomainXmin - Lx) * RxInv
                   proceed = (wtest > 0.0) .and. (wtest <= 1.0) .and. (wtest <= w)

                   if (proceed) then
                       Py = Ly + wtest * Ry
                       onDomainLine = (Py >= compDomainYmin) .and. (Py <= compDomainYmax)

                       if (onDomainLine) then
                           w = wtest
                           rayX = userDomainXmin
                           rayY = min (  Py, userDomainYmax)          ! force rays on user defined domain
                           rayY = max (rayY, userDomainYmin)          !
                       end if
                   end if

                   wtest   = (userDomainXmax - Lx) * RxInv
                   proceed = (wtest > 0.0) .and. (wtest <= 1.0) .and. (wtest <= w)

                   if (proceed) then
                       Py = Ly + wtest * Ry
                       onDomainLine = (Py >= compDomainYmin) .and. (Py <= compDomainYmax)

                       if (onDomainLine) then
                           w = wtest
                           rayX = userDomainXmax
                           rayY = min (  Py, userDomainYmax)
                           rayY = max (rayY, userDomainYmin)
                       end if
                   end if

               end if
!
!
!     ...Next check, if ray hits any of the two boundary x-lines of the domain.
!
!
               if (.not. parallel2Xaxis) then

                   RyInv   = 1.0 / Ry
                   wtest   = (userDomainYmin - Ly) * RyInv
                   proceed = (wtest > 0.0) .and. (wtest <= 1.0) .and. (wtest <= w)

                   if (proceed) then
                       Px = Lx + wtest * Rx
                       onDomainLine = (Px >= compDomainXmin) .and. (Px <= compDomainXmax)

                       if (onDomainLine) then
                           w = wtest
                           rayY = userDomainYmin
                           rayX = min (  Px, userDomainXmax)
                           rayX = max (rayX, userDomainXmin)
                       end if
                   end if

                   wtest = (userDomainYmax - Ly) * RyInv
                   proceed = (wtest > 0.0) .and. (wtest <= 1.0) .and. (wtest <= w)

                   if (proceed) then
                       Px = Lx + wtest * Rx
                       onDomainLine = (Px >= compDomainXmin) .and. (Px <= compDomainXmax)

                       if (onDomainLine) then
                           w = wtest
                           rayY = userDomainYmax
                           rayX = min (  Px, userDomainXmax)
                           rayX = max (rayX, userDomainXmin)
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
                   call Driver_abortFlash("ed_createRays2DRec: ray does not hit domain boundary!")
               end if
!
!
!     ...In case of 2D cylindrical geometries, handle the [-R,0] rays. They must be placed
!        on the [0,R] part (i.e. the computational domain) with the ray vector mirrored against
!        the z-axis.
!
!
               if (grid2Dcylindrical .and. rayX < 0.0) then
                   Rx   = - Rx
                   rayX = - rayX
               end if
!
!
!     ...loop over all blocks and see, if ray is contained in one of it. As soon as that block
!        is found, exit the block loop, since each ray can be assigned to only one block.
!        A ray belongs to a block, if its x,y coordinate is such that:
!
!                        x,y block lower limit  <=  ray x,y  <  x,y block upper limit
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

                  xminLimit = (rayX >= xminBlock)
                  yminLimit = (rayY >= yminBlock)
                  xmaxLimit = (rayX <  xmaxBlock) .or. ((rayX == xmaxBlock) .and. (xmaxBlock == userDomainXmax))
                  ymaxLimit = (rayY <  ymaxBlock) .or. ((rayY == ymaxBlock) .and. (ymaxBlock == userDomainYmax))

                  inBlock  =      xminLimit &
                            .and. xmaxLimit &
                            .and. yminLimit &
                            .and. ymaxLimit

                  if (inBlock) then

                      distToFaceMinX = abs (xminBlock - rayX)
                      distToFaceMaxX = abs (xmaxBlock - rayX)
                      distToFaceMinY = abs (yminBlock - rayY)
                      distToFaceMaxY = abs (ymaxBlock - rayY)

                      blockFaceMinX = (distToFaceMinX < cellWallThicknessHalf)
                      blockFaceMaxX = (distToFaceMaxX < cellWallThicknessHalf)
                      blockFaceMinY = (distToFaceMinY < cellWallThicknessHalf)
                      blockFaceMaxY = (distToFaceMaxY < cellWallThicknessHalf)

                      blockReflectMinX = (blockBoundary (1,block) == REFLECTING)
                      blockReflectMaxX = (blockBoundary (2,block) == REFLECTING)
                      blockReflectMinY = (blockBoundary (3,block) == REFLECTING)
                      blockReflectMaxY = (blockBoundary (4,block) == REFLECTING)

                      rayReflects =     (blockFaceMinX .and. blockReflectMinX) &
                                   .or. (blockFaceMaxX .and. blockReflectMaxX) &
                                   .or. (blockFaceMinY .and. blockReflectMinY) &
                                   .or. (blockFaceMaxY .and. blockReflectMaxY)

                      if (ignoreBoundary .or. .not.rayReflects) then

                          ed_rayCount = ed_rayCount + 1                ! this count is per processor, not block

                          if (ed_rayCount > ed_maxRayCount) then
                              call Driver_abortFlash ("ed_createRays2DRec: Not enough storage for rays array")
                          end if

                          blockID = blockList (block)

                          RsizeInv = 1.0 / sqrt (Rx * Rx + Ry * Ry)

                          uRx = Rx * RsizeInv                          ! unit vector along ray vector
                          uRy = Ry * RsizeInv

                          rayWeight = ed_beamCrossSectionWeight (functionType,               &
                                                                 x        = lineTx,          &
                                                                 Cx       = gaussianCenter,  &
                                                                 Rx       = gaussianRadius,  &
                                                                 Exponent = gaussianExponent )

                          if (grid2Dcylindrical .AND. magnifyT2L*ux.NE.0.0) then
                              rayWeight = rayWeight * abs (Lx)
                          end if

                          powerFraction = rayWeight * gridWeightInv

                          rayPower = beamPower * powerFraction
                          raySpeed = initialRaySpeed * ed_speedOfLight

                          velX = uRx * raySpeed
                          velY = uRy * raySpeed
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

                          ed_rays (RAY_BLCK,ed_rayCount) = real (blockID)
                          ed_rays (RAY_PROC,ed_rayCount) = real (ed_globalMe)
                          ed_rays (RAY_POSX,ed_rayCount) = rayX
                          ed_rays (RAY_POSY,ed_rayCount) = rayY
                          ed_rays (RAY_POSZ,ed_rayCount) = 0.0
                          ed_rays (RAY_VELX,ed_rayCount) = velX
                          ed_rays (RAY_VELY,ed_rayCount) = velY
                          ed_rays (RAY_VELZ,ed_rayCount) = 0.0
                          ed_rays (RAY_POWR,ed_rayCount) = rayPower
                          ed_rays (RAY_DENC,ed_rayCount) = beamCriticalDensity

                          ed_energyInTimestep = ed_energyInTimestep + rayPower * timeStep

                          exit

                      end if      ! ray reflection condition
                   end if         ! ray in block
               end do             ! block loop
            end do                ! individual grid points loop

            usedGridPoints = usedGridPoints + nGridPoints

         end do                   ! array of grid points loop

         if (usedGridPoints /= nRaysBeam) then
             call Driver_abortFlash ("ed_createRays2DRec: # of rays / beam grid points mismatch !")
         end if

     end if                       ! active beam
  end do                          ! beam loop
!
!
!     ...Deallocate block bounding box and block boundary storage array.
!
!
  deallocate (xGrid)
  deallocate (blockBndBox)
  deallocate (blockBoundary)
!
!
!     ...Ready!
!
!
  return
end subroutine ed_createRays2DRec
