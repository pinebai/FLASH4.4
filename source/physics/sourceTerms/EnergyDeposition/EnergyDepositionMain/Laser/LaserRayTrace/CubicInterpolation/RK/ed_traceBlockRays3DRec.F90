!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/RK/ed_traceBlockRays3DRec
!!
!! NAME
!!
!!  ed_traceBlockRays3DRec
!!
!! SYNOPSIS
!!
!!  call ed_traceBlockRays3DRec (real    (in)    :: timeStep,
!!                               integer (in)    :: rayFirst
!!                               integer (in)    :: rayLast,
!!                               integer (in)    :: iminBlock,
!!                               integer (in)    :: imaxBlock,
!!                               integer (in)    :: jminBlock,
!!                               integer (in)    :: jmaxBlock,
!!                               integer (in)    :: kminBlock,
!!                               integer (in)    :: kmaxBlock,
!!                               real    (in)    :: xminBlock,
!!                               real    (in)    :: xmaxBlock,
!!                               real    (in)    :: yminBlock,
!!                               real    (in)    :: ymaxBlock,
!!                               real    (in)    :: zminBlock,
!!                               real    (in)    :: zmaxBlock,
!!                               real    (in)    :: deltaX,
!!                               real    (in)    :: deltaY,
!!                               real    (in)    :: deltaZ,
!!                               real    (in)    :: deltaInvX,
!!                               real    (in)    :: deltaInvY,
!!                               real    (in)    :: deltaInvZ,
!!                               logical (in)    :: blockReflectMinX,
!!                               logical (in)    :: blockReflectMaxX,
!!                               logical (in)    :: blockReflectMinY,
!!                               logical (in)    :: blockReflectMaxY,
!!                               logical (in)    :: blockReflectMinZ,
!!                               logical (in)    :: blockReflectMaxZ,
!!                               real    (inout) :: cellEnergyDepot (:,:,:),
!!                               real    (inout) :: cellIntensityDepot (:,:,:))
!!
!! DESCRIPTION
!!
!!  Traces the movement of the current collection of active rays through one block
!!  for those geometries consisting formally of 3D rectangular grids (cartesian).
!!  On exit, each ray has either:
!!
!!            i)  reached a different (yet unknown) block
!!           ii)  has been absorbed -> is nonexistent
!!          iii)  has reached the domain boundary and exited -> is nonexistent.
!!
!!  The ray tracing is performed using Runge Kutta integration within each cell.
!!
!! ARGUMENTS
!!
!!  timeStep         : current timestep value
!!  rayFirst         : first ray index to be considered
!!  rayLast          : last ray index to be considered
!!  iminBlock        : minimum cell i-index limit defining the interior block
!!  imaxBlock        : maximum cell i-index limit defining the interior block
!!  jminBlock        : minimum cell j-index limit defining the interior block
!!  jmaxBlock        : maximum cell j-index limit defining the interior block
!!  kminBlock        : minimum cell k-index limit defining the interior block
!!  kmaxBlock        : maximum cell k-index limit defining the interior block
!!  xminBlock        : minimum x-coordinate limit of the block
!!  xmaxBlock        : maximum x-coordinate limit of the block
!!  yminBlock        : minimum y-coordinate limit of the block
!!  ymaxBlock        : maximum y-coordinate limit of the block
!!  zminBlock        : minimum z-coordinate limit of the block
!!  zmaxBlock        : maximum z-coordinate limit of the block
!!  deltaX           : the cell's x-dimension
!!  deltaY           : the cell's y-dimension
!!  deltaZ           : the cell's z-dimension
!!  deltaInvX        : inverse of the cell's x-dimension
!!  deltaInvY        : inverse of the cell's y-dimension
!!  deltaInvZ        : inverse of the cell's z-dimension
!!  blockReflectMinX : is the block boundary on the minimum x-side reflective ?
!!  blockReflectMaxX : is the block boundary on the maximum x-side reflective ?
!!  blockReflectMinY : is the block boundary on the minimum y-side reflective ?
!!  blockReflectMaxY : is the block boundary on the maximum y-side reflective ?
!!  blockReflectMinZ : is the block boundary on the minimum z-side reflective ?
!!  blockReflectMaxZ : is the block boundary on the maximum z-side reflective ?
!!  cellEnergyDepot  : array collecting the ray energy deposition for each cell
!!
!! NOTES
!!        
!!  The code allows for threading to be used on the outer ray trace loop.
!!  The paths of the rays are computed using the geometric optics approximation.
!!  Inside each cell, the paths of the rays are evaluated based on Runge Kutta
!!  integration, using the tricubic expansions of the number electron density and
!!  electron temperature grid.
!!
!!***

subroutine ed_traceBlockRays3DRec (timeStep,                          &
                                   rayFirst,  rayLast,                &
                                   iminBlock, imaxBlock,              &
                                   jminBlock, jmaxBlock,              &
                                   kminBlock, kmaxBlock,              &
                                   xminBlock, xmaxBlock,              &
                                   yminBlock, ymaxBlock,              &
                                   zminBlock, zmaxBlock,              &
                                   deltaX,    deltaY,    deltaZ,      &
                                   deltaInvX, deltaInvY, deltaInvZ,   &
                                   blockReflectMinX,                  &
                                   blockReflectMaxX,                  &
                                   blockReflectMinY,                  &
                                   blockReflectMaxY,                  &
                                   blockReflectMinZ,                  &
                                   blockReflectMaxZ,                  &
                                                     cellEnergyDepot, &
                                                   cellIntensityDepot ) 

  use EnergyDeposition_data,      ONLY : ed_Boltzmann,                   &
                                         ed_cellCubicNele,               &
                                         ed_cellCubicTele,               &
                                         ed_cellDensity,                 &
                                         ed_cellEdges,                   &
                                         ed_cellStepTolerance,           &
                                         ed_cellWallThickness,           &
                                         ed_cellZbar,                    &
                                         ed_depoVarIsPerMass,            &
                                         ed_electronMass,                &
                                         ed_electronCharge,              &
                                         ed_energyOutTimestep,           &
                                         ed_infiniteSpeed,               &
                                         ed_infiniteTime,                &
                                     ed_irradVar,                    &
                                         ed_laserIOMaxNumberOfPositions, &
                                         ed_laserIOMaxNumberOfRays,      &
                                         ed_laserIONumberOfPositions,    &
                                         ed_laserIONumberOfRaysWritten,  &
                                         ed_laserIORayFrequency,         &
                                         ed_laserIORayPositions,         &
                                         ed_laserIORayPower,             &
                                         ed_laserIORayTags,              &
                                         ed_laserIOWrite,                &
                                         ed_powerStepTolerance,          &
                                         ed_rays,                        &
                                         ed_raysMovedIntoDomain,         &
                                         ed_rayZeroPower,                &
                                         ed_RungeKuttaMethod,            &
                                         ed_speedOfLightSquared,         &
                                         ed_unitRoundoff,                &
                                         ed_xminDomain,                  &
                                         ed_xmaxDomain,                  &
                                         ed_yminDomain,                  &
                                         ed_ymaxDomain,                  &
                                         ed_zminDomain,                  &
                                         ed_zmaxDomain

  use ed_raytraceODEfunctionData, ONLY : accX, accY, accZ,                         &
                                         accFactorX, accFactorY, accFactorZ,       &
                                         cellEdgeInvX, cellEdgeInvY, cellEdgeInvZ, &
                                         cellZbar,                                 &
                                         i, j, k,                                  &
                                         Nele,                                     &
                                         rayCritDens,                              &
                                         saveComputations,                         &
                                         x01, y01, z01,                            &
                                         xmaxCell, ymaxCell, zmaxCell,             &
                                         xminCell, yminCell, zminCell

  use Driver_interface,           ONLY : Driver_abortFlash

  use Interpolate_interface,      ONLY : Interpolate_cubic3DF,   &
                                         Interpolate_cubic3DFd1

  use RungeKutta_interface,       ONLY : RungeKutta_stepConfined

  use ed_interface,               ONLY : ed_maxConfinement3DRec,       &
                                         ed_minConfinement3DRec,       &
                                         ed_raytraceODEfunction3DRec,  &
                                         ed_time2FacesParabolicPath1D

  use ed_commInterface,           ONLY : ed_commHandleOffBlkRay,   &
                                         ed_commIncrementDeadRays, &
                                         ed_commProgressTransport

  implicit none

#include "constants.h"
#include "Flash.h"
#include "EnergyDeposition.h"

  real,    intent (in)    :: timeStep
  integer, intent (in)    :: rayFirst,  rayLast   
  integer, intent (in)    :: iminBlock, imaxBlock
  integer, intent (in)    :: jminBlock, jmaxBlock
  integer, intent (in)    :: kminBlock, kmaxBlock
  real,    intent (in)    :: xminBlock, xmaxBlock
  real,    intent (in)    :: yminBlock, ymaxBlock
  real,    intent (in)    :: zminBlock, zmaxBlock
  real,    intent (in)    :: deltaX,    deltaY,    deltaZ
  real,    intent (in)    :: deltaInvX, deltaInvY, deltaInvZ
  logical, intent (in)    :: blockReflectMinX
  logical, intent (in)    :: blockReflectMaxX
  logical, intent (in)    :: blockReflectMinY
  logical, intent (in)    :: blockReflectMaxY
  logical, intent (in)    :: blockReflectMinZ
  logical, intent (in)    :: blockReflectMaxZ
  real,    intent (inout) :: cellEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock,kminBlock:kmaxBlock)
  real,    intent (inout),OPTIONAL &
                          :: cellIntensityDepot (iminBlock:      ,jminBlock:         ,kminBlock:)

  logical :: badTimeStep
  logical :: blockFaceMinX, blockFaceMaxX
  logical :: blockFaceMinY, blockFaceMaxY
  logical :: blockFaceMinZ, blockFaceMaxZ
  logical :: cellFaceMinX,  cellFaceMaxX
  logical :: cellFaceMinY,  cellFaceMaxY
  logical :: cellFaceMinZ,  cellFaceMaxZ
  logical :: crossX, crossY, crossZ
  logical :: inDomain, inBlock
  logical :: newCell, outOfCell
  logical :: onBlockBoundaryCell
  logical :: rayOutOfBlock
  logical :: reflectX, reflectY, reflectZ
  logical :: stationaryRay
  logical :: velXeq0, velXgt0, velXlt0
  logical :: velYeq0, velYgt0, velYlt0
  logical :: velZeq0, velZgt0, velZlt0
  logical :: writeRay

  integer :: ip,jp,kp
  integer :: n
  integer :: nRayWritePositions
  integer :: numDeadRays
  integer :: ray
  integer :: rayTag
  integer :: rayWriteIndex

  real    :: c2div2nc
  real    :: cellDensity
  real    :: cellEnergy
  real    :: cellMass, cellMassInv
  real    :: cellPower
  real    :: cellRadEnergy, cellRadPower
  real    :: cellStepErrorX, cellStepErrorY, cellStepErrorZ
  real    :: cellVolume, cellVolumeInv
  real    :: cellWallThicknessHalf
  real    :: dist2minX, dist2minY, dist2minZ
  real    :: dist2maxX, dist2maxY, dist2maxZ
  real    :: minDistance
  real    :: nudgeX, nudgeY, nudgeZ
  real    :: rayErrorFrac
  real    :: rayPower
  real    :: rayEner
  real    :: rayX, rayY, rayZ
  real    :: stepTimeTry, stepTimeUsed, stepTimeNext
  real    :: tx, ty, tz
  real    :: velX, velY, velZ

  real    :: rayError     (1:8)
  real    :: rayErrorBase (1:8)
  real    :: rayIn        (1:8)
  real    :: rayOut       (1:8)
!
!
!     ...Define some variables. Fix the cell volume, which is independent of the cell location
!        inside the block. The error base for the ray power has to be set inside the deposition
!        loop, as it depends on the current ray power.
!
!
  numDeadRays = 0

  cellVolume            = deltaX * deltaY * deltaZ
  cellVolumeInv         = 1.0 / cellVolume
  cellWallThicknessHalf = 0.5 * ed_cellWallThickness
  cellStepErrorX        = ed_cellStepTolerance * deltaX
  cellStepErrorY        = ed_cellStepTolerance * deltaY
  cellStepErrorZ        = ed_cellStepTolerance * deltaZ

  rayErrorFrac = 1.0                        ! this means the error is controlled by the base only

  rayErrorBase (1) = cellStepErrorX         ! error bar on rayX
  rayErrorBase (2) = cellStepErrorY         ! error bar on rayY
  rayErrorBase (3) = cellStepErrorZ         ! error bar on rayZ
  rayErrorBase (4) = ed_infiniteSpeed       ! this in effect puts no error bars on velX yet
  rayErrorBase (5) = ed_infiniteSpeed       ! this in effect puts no error bars on velY yet
  rayErrorBase (6) = ed_infiniteSpeed       ! this in effect puts no error bars on velZ yet
!
!
!     ...Outer (threaded) loop over all rays associated with the current block.
!
!
!$omp do schedule (dynamic)
  do ray = rayFirst , rayLast

     call ed_commProgressTransport ()

     rayTag      = int (ed_rays (RAY_TAGS,ray))
     rayX        =      ed_rays (RAY_POSX,ray)
     rayY        =      ed_rays (RAY_POSY,ray)
     rayZ        =      ed_rays (RAY_POSZ,ray)
     velX        =      ed_rays (RAY_VELX,ray)
     velY        =      ed_rays (RAY_VELY,ray)
     velZ        =      ed_rays (RAY_VELZ,ray)
     rayPower    =      ed_rays (RAY_POWR,ray)
     rayCritDens =      ed_rays (RAY_DENC,ray)
     rayEner     =      0.0

     c2div2nc = ed_speedOfLightSquared / (rayCritDens + rayCritDens)

     accFactorX = - c2div2nc * deltaInvX
     accFactorY = - c2div2nc * deltaInvY
     accFactorZ = - c2div2nc * deltaInvZ
!
!
!     ...Decide, if we should write this ray out. If the case, start the writeout procedure.
!        If threading is done, the ray writing index must be protected from incrementation
!        by another thread and is saved in a local thread variable.
!
!
     writeRay = ed_laserIOWrite

     if (writeRay) then
         writeRay = mod (rayTag, ed_laserIORayFrequency) == 0
         if (writeRay) then
             !$omp critical (WriteRayIndex)
                   writeRay = ed_laserIONumberOfRaysWritten < ed_laserIOMaxNumberOfRays
                   if (writeRay) then
                       ed_laserIONumberOfRaysWritten = ed_laserIONumberOfRaysWritten + 1
                       rayWriteIndex = ed_laserIONumberOfRaysWritten
                   end if
             !$omp end critical (WriteRayIndex)
         end if
     end if

     if(writeRay) then
        nRayWritePositions = 0
        ed_laserIORayTags (rayWriteIndex) = rayTag
     end if
!
!
!     ...Find the indices (i,j,k) of the initial cell through which the ray will
!        enter the block. We know for sure that the ray enters the block, because
!        otherwise it would not be on the current block list. Check, on which of
!        the six possible faces the ray currently is.
!
!
     rayOutOfBlock =     (rayX < xminBlock) &
                    .or. (rayX > xmaxBlock) &
                    .or. (rayY < yminBlock) &
                    .or. (rayY > ymaxBlock) &
                    .or. (rayZ < zminBlock) &
                    .or. (rayZ > zmaxBlock)

     if (rayOutOfBlock) then
         call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: ray found out of block')
     end if

     dist2minX = rayX - xminBlock
     dist2maxX = xmaxBlock - rayX
     dist2minY = rayY - yminBlock
     dist2maxY = ymaxBlock - rayY
     dist2minZ = rayZ - zminBlock
     dist2maxZ = zmaxBlock - rayZ

     minDistance = min (dist2minX, dist2maxX, &
                        dist2minY, dist2maxY, &
                        dist2minZ, dist2maxZ)

     if (minDistance > ed_cellWallThickness) then
         call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: ray too far inside the block')
     end if

     i = iminBlock + int ( (rayX - xminBlock) * deltaInvX )
     j = jminBlock + int ( (rayY - yminBlock) * deltaInvY )
     k = kminBlock + int ( (rayZ - zminBlock) * deltaInvZ )

     onBlockBoundaryCell = (     (i == iminBlock) &
                            .or. (i == imaxBlock) &
                            .or. (j == jminBlock) &
                            .or. (j == jmaxBlock) &
                            .or. (k == kminBlock) &
                            .or. (k == kmaxBlock) )

     if (.not.onBlockBoundaryCell) then
          call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: ray not in a block boundary cell')
     end if

     xminCell = ed_cellEdges (i  ,1)
     xmaxCell = ed_cellEdges (i+1,1)
     yminCell = ed_cellEdges (j  ,2)
     ymaxCell = ed_cellEdges (j+1,2)
     zminCell = ed_cellEdges (k  ,3)
     zmaxCell = ed_cellEdges (k+1,3)

     dist2minX = abs (xminCell - rayX)
     dist2maxX = abs (xmaxCell - rayX)
     dist2minY = abs (yminCell - rayY)
     dist2maxY = abs (ymaxCell - rayY)
     dist2minZ = abs (zminCell - rayZ)
     dist2maxZ = abs (zmaxCell - rayZ)

     cellFaceMinX = (dist2minX <= cellWallThicknessHalf)
     cellFaceMaxX = (dist2maxX <= cellWallThicknessHalf)
     cellFaceMinY = (dist2minY <= cellWallThicknessHalf)
     cellFaceMaxY = (dist2maxY <= cellWallThicknessHalf)
     cellFaceMinZ = (dist2minZ <= cellWallThicknessHalf)
     cellFaceMaxZ = (dist2maxZ <= cellWallThicknessHalf)
!
!
!     ...Make sure the ray is also properly nudged into the corresponding cell.
!
!
     if (cellFaceMinX) rayX = xminCell + cellWallThicknessHalf
     if (cellFaceMaxX) rayX = xmaxCell - cellWallThicknessHalf
     if (cellFaceMinY) rayY = yminCell + cellWallThicknessHalf
     if (cellFaceMaxY) rayY = ymaxCell - cellWallThicknessHalf
     if (cellFaceMinZ) rayZ = zminCell + cellWallThicknessHalf
     if (cellFaceMaxZ) rayZ = zmaxCell - cellWallThicknessHalf

     velXeq0 = (velX == 0.0)
     velYeq0 = (velY == 0.0)
     velZeq0 = (velZ == 0.0)

     stationaryRay = (velXeq0 .and. velYeq0 .and. velZeq0)

     if (stationaryRay) then
         call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: stationary ray at a block face boundary')
     end if
!
!
!     ...Get extra needed info about the initial cell (i,j,k). The calculation of the inverse
!        length of the cell edges have to be based on the actual min/max cell values. If
!        based on the (global per block) cell delta values, roundoff errors introduced during
!        evaluation of the min/max cell values can lead to abortion of the cubic interpolation
!        routines, which strictly! require [0,1] rescaled coordinates.
!
!
     cellZbar     = ed_cellZbar    (i,j,k)
     cellDensity  = ed_cellDensity (i,j,k)
     cellMass     = cellDensity * cellVolume
     cellMassInv  = 1.0 / cellMass
     cellEdgeInvX = 1.0 / (xmaxCell - xminCell)
     cellEdgeInvY = 1.0 / (ymaxCell - yminCell)
     cellEdgeInvZ = 1.0 / (zmaxCell - zminCell)
!
!
!     ...We are ready to follow the ray's path through all the cells of the current
!        block. At this point, the ray is in the initial cell (i,j,k) and is ready to move
!        through the block. In case a laser IO is performed on the rays, store the initial
!        ray IO data.
!
!
     if (writeRay) then
         nRayWritePositions = nRayWritePositions + 1
         if(nRayWritePositions <= ed_laserIOMaxNumberOfPositions) then
            ed_laserIONumberOfPositions (rayWriteIndex                           ) = nRayWritePositions
            ed_laserIORayPower          (rayWriteIndex, nRayWritePositions       ) = rayPower
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, IAXIS) = rayX
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = rayY
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = rayZ
         end if
     end if
!
!
!-------------------- Loop following ray through cells in block --------------------------------------------
!
!
     stepTimeNext = ed_infiniteTime

     do                                ! indefinite loop through the block cells
                                       ! will be broken (exit) by the various conditions
                                       ! of the ray (no power, out of domain, etc)
!
!
!     ...From the current position, velocity and accelleration of the ray, we determine
!        the initial stepping time to the closest cell wall by assuming constant acceleration.
!        This will be the maximal step time used for the Runge Kutta stepper.
!
!
        x01 = (rayX - xminCell) * cellEdgeInvX                          ! rescaled [0,1] ray x coordinate
        y01 = (rayY - yminCell) * cellEdgeInvY                          ! rescaled [0,1] ray y coordinate
        z01 = (rayZ - zminCell) * cellEdgeInvZ                          ! rescaled [0,1] ray z coordinate

        Nele (1:4) = Interpolate_cubic3DFd1 (ed_cellCubicNele (1:64,i,j,k), x01,y01,z01)

        accX = accFactorX * Nele (2)                                    ! acceleration in x-direction
        accY = accFactorY * Nele (3)                                    ! acceleration in y-direction
        accZ = accFactorZ * Nele (4)                                    ! acceleration in z-direction

        tx = ed_time2FacesParabolicPath1D (rayX, velX, accX, xminCell, xmaxCell, ed_infiniteTime)
        ty = ed_time2FacesParabolicPath1D (rayY, velY, accY, yminCell, ymaxCell, ed_infiniteTime)
        tz = ed_time2FacesParabolicPath1D (rayZ, velZ, accZ, zminCell, zmaxCell, ed_infiniteTime)

        stepTimeTry = min (tx, ty, tz, stepTimeNext)                    ! initial time step for RK

        if (stepTimeTry == ed_infiniteTime .or. stepTimeTry == 0.0) then
            call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: infinite/zero stepping time for a cell')
        end if
!
!
!     ...Assemble the vector of dependent variables to be passed to the RK stepper. Since this
!        will be a confined RK step, we also pass the minimum/maximum cell wall coordinates.
!        The RK stepper has to stay within the cell boundaries until one of the cell walls is
!        hit.
!
!        Before calling the RK stepper, set the save computations keyword to true. This signals
!        the ODE function routine to use the current acc(X,Y,Z) and (x,y,z)01 values for evaluation
!        of the first RK point. It avoids the rather costly recalculation of the Nele values and
!        their derivatives at the current ray location. The keyword is set imediately to false
!        after that (inside the ODE function routine), because the RK stepper saves an internal
!        copy of the first ODE function vector.
!
!
        rayIn (1) = rayX
        rayIn (2) = rayY
        rayIn (3) = rayZ
        rayIn (4) = velX
        rayIn (5) = velY
        rayIn (6) = velZ
        rayIn (7) = rayPower
        rayIn (8) = 0

        rayErrorBase (7) = rayPower * ed_powerStepTolerance   ! power error base depends on current power
        rayErrorBase (0) = 0.0

        saveComputations = .true.

        if (ed_irradVar > 0) then
           rayIn (8) = rayEner
           rayErrorBase (8) = max(rayPower*stepTimeTry,rayEner) * ed_powerStepTolerance   ! ener error base depends on current ener
           call  RungeKutta_stepConfined (ed_RungeKuttaMethod,         &
                                       ed_raytraceODEfunction3DRec, &
                                       3,                           &  ! # of confined variables
                                       0.0,                         &  ! dummy time -> ODE independent of time
                                       rayIn        (1:8),          &
                                       ed_minConfinement3DRec,      &  ! lower limit confinement function
                                       ed_maxConfinement3DRec,      &  ! upper limit confinement function
                                       rayErrorFrac,                &
                                       rayErrorBase (1:8),          &
                                       stepTimeTry,                 &
                                       stepTimeUsed,                &  ! actual stepping time used
                                       stepTimeNext,                &  ! estimates of next stepping time
                                       rayOut       (1:8),          &
                                       rayError     (1:8)           )  ! can be +ve or -ve

        else
           rayError (8) = 0.0
           call  RungeKutta_stepConfined (ed_RungeKuttaMethod,         &
                                       ed_raytraceODEfunction3DRec, &
                                       3,                           &  ! # of confined variables
                                       0.0,                         &  ! dummy time -> ODE independent of time
                                       rayIn        (1:7),          &
                                       ed_minConfinement3DRec,      &  ! lower limit confinement function
                                       ed_maxConfinement3DRec,      &  ! upper limit confinement function
                                       rayErrorFrac,                &
                                       rayErrorBase (1:7),          &
                                       stepTimeTry,                 &
                                       stepTimeUsed,                &  ! actual stepping time used
                                       stepTimeNext,                &  ! estimates of next stepping time
                                       rayOut       (1:7),          &
                                       rayError     (1:7)           )  ! can be +ve or -ve
        end if
!
!
!     ...The confined RK step has been taken. Check, if the errors are within the error bars and check
!        where the ray is currently located. Add the lost ray power to the cell energy and update the
!        (diminished) power of the ray.
!
!
        if (any (abs (rayError (1:8)) > rayErrorFrac * rayErrorBase (1:8))) then
            call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: error(s) in RK step out of bounds!')
        end if

        rayX = rayOut (1)
        rayY = rayOut (2)
        rayZ = rayOut (3)
        velX = rayOut (4)
        velY = rayOut (5)
        velZ = rayOut (6)

        outOfCell =     (rayX < xminCell - cellWallThicknessHalf) &    ! for debugging purposes
                   .or. (rayX > xmaxCell + cellWallThicknessHalf) &    ! will be removed once the
                   .or. (rayY < yminCell - cellWallThicknessHalf) &    ! code (the confined RK stepper)
                   .or. (rayY > ymaxCell + cellWallThicknessHalf) &    ! is running properly
                   .or. (rayZ < zminCell - cellWallThicknessHalf) &
                   .or. (rayZ > zmaxCell + cellWallThicknessHalf)

        if (outOfCell) then
            call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: RK stepped out of cell confinement!')
        end if

        newCell =     (rayX < xminCell + cellWallThicknessHalf) &      ! checks, if the current
                 .or. (rayX > xmaxCell - cellWallThicknessHalf) &      ! confined RK step has hit
                 .or. (rayY < yminCell + cellWallThicknessHalf) &      ! one of the cell walls
                 .or. (rayY > ymaxCell - cellWallThicknessHalf) &
                 .or. (rayZ < zminCell + cellWallThicknessHalf) &
                 .or. (rayZ > zmaxCell - cellWallThicknessHalf)

        velXeq0 = (velX == 0.0)
        velYeq0 = (velY == 0.0)
        velZeq0 = (velZ == 0.0)
 
        stationaryRay = velXeq0 .and. velYeq0 .and. velZeq0

        cellPower  = rayPower - rayOut (7)
        cellEnergy = cellPower * timeStep

        if (ed_depoVarIsPerMass) then
            cellEnergyDepot (i,j,k) = cellEnergyDepot (i,j,k) + cellEnergy * cellMassInv
        else
            cellEnergyDepot (i,j,k) = cellEnergyDepot (i,j,k) + cellEnergy * cellVolumeInv
        end if

        rayPower = rayOut (7)

        if (ed_irradVar > 0) then
           rayEner = rayOut (8)

           if (newCell .OR. (rayPower <= ed_rayZeroPower) .OR. stationaryRay) then
              cellRadEnergy   = rayEner
              cellIntensityDepot (i,j,k) = cellIntensityDepot (i,j,k) + cellRadEnergy * cellVolumeInv
              rayEner = 0.0
           end if
        end if

        if (rayPower <= ed_rayZeroPower) then
            ed_rays (RAY_BLCK,ray) = real (NONEXISTENT)
            exit
        end if
!
!
!     ...If the ray is stationary (no movement), mark the ray as nonexistent and exit the
!        block loop. In case a laser IO is performed on the rays, store the current ray IO data.
!
!

        if (stationaryRay) then
            write (*,*) ' stationary ray detected! Removing it from list ... '
            ed_rays (RAY_BLCK,ray) = real (NONEXISTENT)
            exit
        end if

        if (writeRay) then
            nRayWritePositions = nRayWritePositions + 1
            if(nRayWritePositions <= ed_laserIOMaxNumberOfPositions) then
               ed_laserIONumberOfPositions (rayWriteIndex                           ) = nRayWritePositions
               ed_laserIORayPower          (rayWriteIndex, nRayWritePositions       ) = rayPower
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, IAXIS) = rayX
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = rayY
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = rayZ
            end if
        end if
!
!
!     ...If, at the current stage, the ray enters a new cell, we have to determine: 1) which new
!        cell (i,j,k) it is and 2) the appropriate nudging values on the ray's position. Due to
!        possible reflective boundary conditions on the block faces, it can happen that the ray
!        stays in the original cell. After handling the logistics inside the following 'if'
!        statement, the new cell indices i,j,k are either the old ones or new ones.
!
!
        if (newCell) then

            dist2minX = abs (xminCell - rayX)
            dist2maxX = abs (xmaxCell - rayX)
            dist2minY = abs (yminCell - rayY)
            dist2maxY = abs (ymaxCell - rayY)
            dist2minZ = abs (zminCell - rayZ)
            dist2maxZ = abs (zmaxCell - rayZ)

            minDistance = min (dist2minX, dist2maxX, &
                               dist2minY, dist2maxY, &
                               dist2minZ, dist2maxZ)

            if (minDistance > cellWallThicknessHalf) then
                call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: ray to far away from cell face')
            end if

            cellFaceMinX = (dist2minX <= cellWallThicknessHalf)
            cellFaceMaxX = (dist2maxX <= cellWallThicknessHalf)
            cellFaceMinY = (dist2minY <= cellWallThicknessHalf)
            cellFaceMaxY = (dist2maxY <= cellWallThicknessHalf)
            cellFaceMinZ = (dist2minZ <= cellWallThicknessHalf)
            cellFaceMaxZ = (dist2maxZ <= cellWallThicknessHalf)

            velXgt0 = (velX  > 0.0)
            velXlt0 = (velX  < 0.0)
            velYgt0 = (velY  > 0.0)
            velYlt0 = (velY  < 0.0)
            velZgt0 = (velZ  > 0.0)
            velZlt0 = (velZ  < 0.0)

            crossX = .false.
            crossY = .false.
            crossZ = .false.

            nudgeX = 0.0
            nudgeY = 0.0
            nudgeZ = 0.0

            ip = i
            jp = j
            kp = k

            if (cellFaceMinX) then

                rayX   = xminCell
                nudgeX = + cellWallThicknessHalf

                if (velXlt0) then
                    i = i - 1
                    crossX = .true.
                end if

            else if (cellFaceMaxX) then

                rayX   = xmaxCell
                nudgeX = - cellWallThicknessHalf

                if (velXgt0) then
                    i = i + 1
                    crossX = .true.
                end if

            end if

            if (cellFaceMinY) then

                rayY   = yminCell
                nudgeY = + cellWallThicknessHalf

                if (velYlt0) then
                    j = j - 1
                    crossY = .true.
                end if

            else if (cellFaceMaxY) then

                rayY   = ymaxCell
                nudgeY = - cellWallThicknessHalf

                if (velYgt0) then
                    j = j + 1
                    crossY = .true.
                end if

            end if

            if (cellFaceMinZ) then

                rayZ   = zminCell
                nudgeZ = + cellWallThicknessHalf

                if (velZlt0) then
                    k = k - 1
                    crossZ = .true.
                end if

            else if (cellFaceMaxZ) then

                rayZ   = zmaxCell
                nudgeZ = - cellWallThicknessHalf

                if (velZgt0) then
                    k = k + 1
                    crossZ = .true.
                end if

            end if

            blockFaceMinX = (rayX == xminBlock)
            blockFaceMaxX = (rayX == xmaxBlock)
            blockFaceMinY = (rayY == yminBlock)
            blockFaceMaxY = (rayY == ymaxBlock)
            blockFaceMinZ = (rayZ == zminBlock)
            blockFaceMaxZ = (rayZ == zmaxBlock)

            reflectX =     (blockFaceMinX .and. blockReflectMinX .and. velXlt0) &
                      .or. (blockFaceMaxX .and. blockReflectMaxX .and. velXgt0)
            reflectY =     (blockFaceMinY .and. blockReflectMinY .and. velYlt0) &
                      .or. (blockFaceMaxY .and. blockReflectMaxY .and. velYgt0)
            reflectZ =     (blockFaceMinZ .and. blockReflectMinZ .and. velZlt0) &
                      .or. (blockFaceMaxZ .and. blockReflectMaxZ .and. velZgt0)

            if (reflectX) then
                i = ip
                velX = - velX
                crossX = .false.
            end if

            if (reflectY) then
                j = jp
                velY = - velY
                crossY = .false.
            end if

            if (reflectZ) then
                k = kp
                velZ = - velZ
                crossZ = .false.
            end if

            if (crossX) then
                nudgeX = (i - ip) * cellWallThicknessHalf
            end if

            if (crossY) then
                nudgeY = (j - jp) * cellWallThicknessHalf
            end if

            if (crossZ) then
                nudgeZ = (k - kp) * cellWallThicknessHalf
            end if

            rayX = rayX + nudgeX
            rayY = rayY + nudgeY
            rayZ = rayZ + nudgeZ

            newCell = crossX .or. crossY .or. crossZ

        end if
!
!
!     ...We are now sure about the target cell. Check, if the target cell (i,j,k) is still within the block.
!        If it is, we check if this is a new cell, in which case we update the cell info. If the target cell
!        is not within the block, check if the ray coordinates are still within the defined domain. If not,
!        store its latest data and mark it as nonexistent. If the ray is still within the domain boundaries,
!        exit the current block loop.
!
!
        inBlock =      (i >= iminBlock) &
                 .and. (i <= imaxBlock) &
                 .and. (j >= jminBlock) &
                 .and. (j <= jmaxBlock) &
                 .and. (k >= kminBlock) &
                 .and. (k <= kmaxBlock)

        if (inBlock) then

            if (newCell) then

                xminCell    = ed_cellEdges   (i  ,1)
                xmaxCell    = ed_cellEdges   (i+1,1)
                yminCell    = ed_cellEdges   (j  ,2)
                ymaxCell    = ed_cellEdges   (j+1,2)
                zminCell    = ed_cellEdges   (k  ,3)
                zmaxCell    = ed_cellEdges   (k+1,3)

                cellZbar     = ed_cellZbar    (i,j,k)
                cellDensity  = ed_cellDensity (i,j,k)
                cellMass     = cellDensity * cellVolume
                cellMassInv  = 1.0 / cellMass
                cellEdgeInvX = 1.0 / (xmaxCell - xminCell)
                cellEdgeInvY = 1.0 / (ymaxCell - yminCell)
                cellEdgeInvZ = 1.0 / (zmaxCell - zminCell)

            end if

        else

            ed_rays (RAY_POSX,ray) = rayX
            ed_rays (RAY_POSY,ray) = rayY
            ed_rays (RAY_POSZ,ray) = rayZ
            ed_rays (RAY_VELX,ray) = velX
            ed_rays (RAY_VELY,ray) = velY
            ed_rays (RAY_VELZ,ray) = velZ
            ed_rays (RAY_POWR,ray) = rayPower

            inDomain =      (rayX > ed_xminDomain) &
                      .and. (rayX < ed_xmaxDomain) &
                      .and. (rayY > ed_yminDomain) &
                      .and. (rayY < ed_ymaxDomain) &
                      .and. (rayZ > ed_zminDomain) &
                      .and. (rayZ < ed_zmaxDomain)

            if (.not.inDomain) then
                 ed_rays (RAY_POSX,ray) = rayX - nudgeX             ! undo the nudging
                 ed_rays (RAY_POSY,ray) = rayY - nudgeY             ! since it is not
                 ed_rays (RAY_POSZ,ray) = rayZ - nudgeZ             ! needed anymore
                 ed_rays (RAY_BLCK,ray) = real (RAY_OUTDOMAIN)
                 ed_energyOutTimeStep   = ed_energyOutTimeStep + rayPower * timeStep
                 numDeadRays            = numDeadRays + 1
            else
                 call ed_commHandleOffBlkRay (ray)
            end if

            exit

        end if
!
!
!-------------------- End loop following ray through cells in block --------------------------------------------
!
!
     end do
!
!
!     ...Check to see if we ran out of laser IO buffer space
!
!
     if(writeRay .and. (nRayWritePositions > ed_laserIOMaxNumberOfPositions) ) then
        print *, "[ed_traceBlockRays3DRec] Ray ", ray, &
                 " ran out of IO buffer space (", nRayWritePositions, ")"
     end if
!
!
!     ...Consider next ray.
!
!
  end do
!$omp end do nowait

  if (numDeadRays > 0) then
      call ed_commIncrementDeadRays (numDeadRays)
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_traceBlockRays3DRec
