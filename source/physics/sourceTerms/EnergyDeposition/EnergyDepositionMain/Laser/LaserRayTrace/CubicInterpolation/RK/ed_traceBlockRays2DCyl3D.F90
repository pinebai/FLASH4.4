!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/RK/ed_traceBlockRays2DCyl3D
!!
!! NAME
!!
!!  ed_traceBlockRays2DCyl3D
!!
!! SYNOPSIS
!!
!!  call ed_traceBlockRays2DCyl3D (real    (in)    :: timeStep,
!!                                 integer (in)    :: rayFirst
!!                                 integer (in)    :: rayLast,
!!                                 integer (in)    :: iminBlock,
!!                                 integer (in)    :: imaxBlock,
!!                                 integer (in)    :: jminBlock,
!!                                 integer (in)    :: jmaxBlock,
!!                                 real    (in)    :: xminBlock,
!!                                 real    (in)    :: xmaxBlock,
!!                                 real    (in)    :: zminBlock,
!!                                 real    (in)    :: zmaxBlock,
!!                                 real    (in)    :: deltaX,
!!                                 real    (in)    :: deltaZ,
!!                                 real    (in)    :: deltaInvX,
!!                                 real    (in)    :: deltaInvZ,
!!                                 logical (in)    :: blockReflectMinX,
!!                                 logical (in)    :: blockReflectMaxX,
!!                                 logical (in)    :: blockReflectMinZ,
!!                                 logical (in)    :: blockReflectMaxZ,
!!                                 real    (inout) :: wedgeEnergyDepot (:,:))
!!
!! DESCRIPTION
!!
!!  Traces the movement of the current collection of active 3D cartesian rays through
!!  one block for 2D cylindrical geometries. On exit, each ray has either:
!!
!!            i)  reached a different (yet unknown) 2D cylindrical block
!!           ii)  has been absorbed -> is nonexistent
!!          iii)  has reached the domain boundary and exited -> is nonexistent.
!!
!!  The 3-dimensional shape of the 2D cylindrical geometry is approximated as a 3D wedge,
!!  whose width extends in the y-direction. The y-direction constitutes the linear angular
!!  approximation. The smaller the width of the wedge, the more accurate the cylindrical
!!  representation becomes. Rays hitting on one of the wedge's y-directional boundaries
!!  stay in the same 2D cylindrical block and automatically jump to the other y-directional
!!  boundary, thus mimicking travel inside a cylindrical shell. The y-directional boundaries
!!  are defined by two lines 'y = m * x' and 'y = - m * x' with opposite slope, whose value
!!  depends on the wedges angle of aperture. The following picture shows the wedge with
!!  origin at the 2D cylindrical domain origin:
!!
!!
!!     y axis
!!       |                                                     *  <--- y = + m * x wedge boundary
!!       |                                            *        |
!!       |                                   *        |        |
!!       |                          *        |        |        |
!!       |                 *        |        |        |        |
!!       |        *        |        |        |        |         
!!       O--------|--------|--------|--------|--------|------- X ---------------------> R or x
!!       |        *        |        |        |        |         
!!       |                 *        |        |        |        |  X = 2D cylindrical
!!       |                          *        |        |        |      domain boundary
!!       |                                   *        |        |
!!       |                                            *        |
!!       |                                                     *  <--- y = - m * x wedge boundary
!!
!!
!!
!!  The z-direction is the same for either the 2D cylindrical or the 3D cartesian picture.
!!  The ray tracing is performed using Runge Kutta integration within each cell.
!!
!! ARGUMENTS
!!
!!  timeStep         : current timestep value
!!  rayFirst         : first ray index to be considered
!!  rayLast          : last ray index to be considered
!!  iminBlock        : minimum wedge i-index limit defining the interior block
!!  imaxBlock        : maximum wedge i-index limit defining the interior block
!!  jminBlock        : minimum wedge j-index limit defining the interior block
!!  jmaxBlock        : maximum wedge j-index limit defining the interior block
!!  xminBlock        : minimum x-coordinate limit of the block
!!  xmaxBlock        : maximum x-coordinate limit of the block
!!  zminBlock        : minimum z-coordinate limit of the block
!!  zmaxBlock        : maximum z-coordinate limit of the block
!!  deltaX           : the wedge's x-dimension
!!  deltaZ           : the wedge's z-dimension
!!  deltaInvX        : inverse of the wedge's x-dimension
!!  deltaInvZ        : inverse of the wedge's z-dimension
!!  blockReflectMinX : is the block boundary on the minimum x-side reflective ?
!!  blockReflectMaxX : is the block boundary on the maximum x-side reflective ?
!!  blockReflectMinZ : is the block boundary on the minimum z-side reflective ?
!!  blockReflectMaxZ : is the block boundary on the maximum z-side reflective ?
!!  wedgeEnergyDepot : array collecting the ray energy deposition for each wedge
!!
!! NOTES
!!
!!  The code allows for threading to be used on the outer ray trace loop.
!!  The paths of the rays are computed using the geometric optics approximation. 
!!  Inside each wedge, the paths of the rays are evaluated stepwise, using the
!!  bicubic expansions of the number electron density and electron temperature
!!  grid.
!!
!!***

subroutine ed_traceBlockRays2DCyl3D (timeStep,                          &
                                     rayFirst,  rayLast,                &
                                     iminBlock, imaxBlock,              &
                                     jminBlock, jmaxBlock,              &
                                     xminBlock, xmaxBlock,              &
                                     zminBlock, zmaxBlock,              &
                                     deltaX,    deltaZ,                 &
                                     deltaInvX, deltaInvZ,              &
                                     blockReflectMinX,                  &
                                     blockReflectMaxX,                  &
                                     blockReflectMinZ,                  &
                                     blockReflectMaxZ,                  &
                                                       wedgeEnergyDepot ) 

  use EnergyDeposition_data,      ONLY : ed_Boltzmann,                   &
                                         ed_cellCubicNele,               &
                                         ed_cellCubicTele,               &
                                         ed_cellDensity,                 &
                                         ed_cellEdges,                   &
                                         ed_cellStepTolerance,           &
                                         ed_cellVolume,                  &
                                         ed_cellWallThickness,           &
                                         ed_cellZbar,                    &
                                         ed_depoVarIsPerMass,            &
                                         ed_electronMass,                &
                                         ed_electronCharge,              &
                                         ed_energyOutTimestep,           &
                                         ed_infiniteSpeed,               &
                                         ed_infiniteTime,                &
                                         ed_largestPositiveReal,         &
                                         ed_laser3Din2DwedgeCosine,      &
                                         ed_laser3Din2DwedgeSine,        &
                                         ed_laser3Din2DwedgeSlope,       &
                                         ed_laserIOMaxNumberOfPositions, &
                                         ed_laserIOMaxNumberOfRays,      &
                                         ed_laserIONumberOfPositions,    &
                                         ed_laserIONumberOfRaysWritten,  &
                                         ed_laserIORayPositions,         &
                                         ed_laserIORayPower,             &
                                         ed_laserIORayFrequency,         &
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
                                         ed_ymaxDomain

  use ed_raytraceODEfunctionData, ONLY : accX, accZ,                   &
                                         accFactorX, accFactorZ,       &
                                         i, j,                         &
                                         Nele,                         &
                                         rayCritDens,                  &
                                         saveComputations,             &
                                         wedgeEdgeInvX, wedgeEdgeInvZ, &
                                         wedgeSlope,                   &
                                         wedgeZbar,                    &
                                         x01, z01,                     &
                                         xmaxWedge, zmaxWedge,         &
                                         xminWedge, zminWedge

  use Driver_interface,           ONLY : Driver_abortFlash

  use Interpolate_interface,      ONLY : Interpolate_cubic2DF,   &
                                         Interpolate_cubic2DFd1

  use RungeKutta_interface,       ONLY : RungeKutta_stepConfined

  use ed_interface,               ONLY : ed_maxConfinement2DCyl3D,      &
                                         ed_minConfinement2DCyl3D,      &
                                         ed_raytraceODEfunction2DCyl3D, &
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
  real,    intent (in)    :: xminBlock, xmaxBlock
  real,    intent (in)    :: zminBlock, zmaxBlock
  real,    intent (in)    :: deltaX,    deltaZ
  real,    intent (in)    :: deltaInvX, deltaInvZ
  logical, intent (in)    :: blockReflectMinX
  logical, intent (in)    :: blockReflectMaxX
  logical, intent (in)    :: blockReflectMinZ
  logical, intent (in)    :: blockReflectMaxZ
  real,    intent (inout) :: wedgeEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock)

  logical :: badTimeStep
  logical :: blockFaceMinX, blockFaceMaxX
  logical :: blockFaceMinZ, blockFaceMaxZ
  logical :: crossX, crossZ
  logical :: inDomain, inBlock
  logical :: newWedge, newWedgeIJ, outOfWedge
  logical :: onBlockBoundaryWedge
  logical :: rayOutOfBlock
  logical :: reflectX, reflectZ
  logical :: stationaryRay
  logical :: velXeq0, velXgt0, velXlt0
  logical :: velYeq0, velYgt0, velYlt0
  logical :: velZeq0, velZgt0, velZlt0
  logical :: wedgeFaceMinX, wedgeFaceMaxX
  logical :: wedgeFaceMinY, wedgeFaceMaxY
  logical :: wedgeFaceMinZ, wedgeFaceMaxZ
  logical :: writeRay

  integer :: ip,jp
  integer :: n
  integer :: nRayWritePositions
  integer :: numDeadRays
  integer :: ray
  integer :: rayTag
  integer :: rayWriteIndex

  integer :: rayBlock

  real    :: c2div2nc
  real    :: dist2minX, dist2maxX
  real    :: dist2minY, dist2maxY
  real    :: dist2minZ, dist2maxZ
  real    :: minDistance
  real    :: nudgeX, nudgeZ
  real    :: powerLossFactor
  real    :: rayErrorFrac
  real    :: rayPower
  real    :: rayX, rayY, rayZ
  real    :: stepTimeTry, stepTimeUsed, stepTimeNext
  real    :: tw, tx, ty, tz
  real    :: velX, velY, velZ
  real    :: wedgeCosine, wedgeSine
  real    :: wedgeDensity
  real    :: wedgeEnergy
  real    :: wedgeMass, wedgeMassInv
  real    :: wedgePower
  real    :: wedgeStepErrorX, wedgeStepErrorZ
  real    :: wedgeVolume, wedgeVolumeInv
  real    :: wedgeWallThicknessHalf
  real    :: yminWedge, ymaxWedge

  real    :: rayError     (1:7)
  real    :: rayErrorBase (1:7)
  real    :: rayIn        (1:7)
  real    :: rayOut       (1:7)
!
!
!     ...Define some variables. The error base for the ray power has to be set
!        inside the deposition loop, as it depends on the current ray power.
!
!
  numDeadRays = 0

  wedgeCosine            = ed_laser3Din2DwedgeCosine
  wedgeSine              = ed_laser3Din2DwedgeSine
  wedgeSlope             = ed_laser3Din2DwedgeSlope
  wedgeWallThicknessHalf = ed_cellWallThickness * 0.5
  wedgeStepErrorX        = ed_cellStepTolerance * deltaX
  wedgeStepErrorZ        = ed_cellStepTolerance * deltaZ

  rayErrorFrac = 1.0                          ! this means the error is controlled by the base only

  rayErrorBase (1) = wedgeStepErrorX          ! error bar on rayX
  rayErrorBase (2) = ed_largestPositiveReal   ! this in effect puts no error bars on rayY
  rayErrorBase (3) = wedgeStepErrorZ          ! error bar on rayZ
  rayErrorBase (4) = ed_infiniteSpeed         ! this in effect puts no error bars on velX yet
  rayErrorBase (5) = ed_infiniteSpeed         ! this in effect puts no error bars on velY yet
  rayErrorBase (6) = ed_infiniteSpeed         ! this in effect puts no error bars on velZ yet
!
!
!     ...Outer (threaded) loop over all rays associated with the current 2D cylindrical block.
!
!
!$omp do schedule (dynamic)
  do ray = rayFirst , rayLast

     call ed_commProgressTransport ()

     rayBlock    = int (ed_rays (RAY_BLCK,ray))

     rayTag      = int (ed_rays (RAY_TAGS,ray))
     rayX        =      ed_rays (RAY_POSX,ray)
     rayY        =      ed_rays (RAY_POSY,ray)
     rayZ        =      ed_rays (RAY_POSZ,ray)
     velX        =      ed_rays (RAY_VELX,ray)
     velY        =      ed_rays (RAY_VELY,ray)
     velZ        =      ed_rays (RAY_VELZ,ray)
     rayPower    =      ed_rays (RAY_POWR,ray)
     rayCritDens =      ed_rays (RAY_DENC,ray)

     c2div2nc = ed_speedOfLightSquared / (rayCritDens + rayCritDens)

     accFactorX = - c2div2nc * deltaInvX
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
                   ed_laserIONumberOfRaysWritten = ed_laserIONumberOfRaysWritten + 1
                   rayWriteIndex = ed_laserIONumberOfRaysWritten
             !$omp end critical (WriteRayIndex)
         end if
     end if

     if(writeRay) then
        nRayWritePositions = 0
        ed_laserIORayTags (rayWriteIndex) = rayTag
     end if
!
!
!     ...Find the indices (i,j) of the initial wedge through which the ray will
!        enter the block. We know for sure that the ray enters the block, because
!        otherwise it would not be on the current block list. Check, on which of
!        the four possible 2D cylindrical wedge faces the ray currently is.
!
!
     rayOutOfBlock =     (rayX < xminBlock) &
                    .or. (rayX > xmaxBlock) &
                    .or. (rayZ < zminBlock) &
                    .or. (rayZ > zmaxBlock)

     if (rayOutOfBlock) then
         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: ray found out of block')
     end if

     dist2minX = rayX - xminBlock
     dist2maxX = xmaxBlock - rayX
     dist2minZ = rayZ - zminBlock
     dist2maxZ = zmaxBlock - rayZ

     minDistance = min (dist2minX, dist2maxX, &
                        dist2minZ, dist2maxZ)

     if (minDistance > ed_cellWallThickness) then
         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: ray too far inside the block')
     end if

     i = iminBlock + int ( (rayX - xminBlock) * deltaInvX )
     j = jminBlock + int ( (rayZ - zminBlock) * deltaInvZ )

     onBlockBoundaryWedge = (     (i == iminBlock) &
                             .or. (i == imaxBlock) &
                             .or. (j == jminBlock) &
                             .or. (j == jmaxBlock) )

     if (.not.onBlockBoundaryWedge) then
          call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: ray not in a block boundary wedge')
     end if

     xminWedge = ed_cellEdges (i  ,1)
     xmaxWedge = ed_cellEdges (i+1,1)
     zminWedge = ed_cellEdges (j  ,2)
     zmaxWedge = ed_cellEdges (j+1,2)

     dist2minX = abs (xminWedge - rayX)
     dist2maxX = abs (xmaxWedge - rayX)
     dist2minZ = abs (zminWedge - rayZ)
     dist2maxZ = abs (zmaxWedge - rayZ)

     wedgeFaceMinX = (dist2minX <= wedgeWallThicknessHalf)
     wedgeFaceMaxX = (dist2maxX <= wedgeWallThicknessHalf)
     wedgeFaceMinZ = (dist2minZ <= wedgeWallThicknessHalf)
     wedgeFaceMaxZ = (dist2maxZ <= wedgeWallThicknessHalf)
!
!
!     ...Make sure the ray is also properly nudged into the corresponding wedge.
!
!
     if (wedgeFaceMinX) rayX = xminWedge + wedgeWallThicknessHalf
     if (wedgeFaceMaxX) rayX = xmaxWedge - wedgeWallThicknessHalf
     if (wedgeFaceMinZ) rayZ = zminWedge + wedgeWallThicknessHalf
     if (wedgeFaceMaxZ) rayZ = zmaxWedge - wedgeWallThicknessHalf

     velXeq0 = (velX == 0.0)
     velYeq0 = (velY == 0.0)
     velZeq0 = (velZ == 0.0)

     stationaryRay = (velXeq0 .and. velYeq0 .and. velZeq0)

     if (stationaryRay) then
         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: stationary ray at a block face boundary')
     end if
!
!
!     ...Get extra needed info about the initial wedge (i,j). The calculation of the inverse
!        length of the wedge edges have to be based on the actual min/max wedge values. If
!        based on the (global per block) wedge delta values, roundoff errors introduced during
!        evaluation of the min/max wedge values can lead to abortion of the cubic interpolation
!        routines, which strictly! require [0,1] rescaled coordinates.
!
!
     wedgeZbar      = ed_cellZbar    (i,j,1)
     wedgeDensity   = ed_cellDensity (i,j,1)
     wedgeVolume    = ed_cellVolume  (i,j,1)
     wedgeVolumeInv = 1.0 / wedgeVolume
     wedgeMass      = wedgeDensity * wedgeVolume
     wedgeMassInv   = 1.0 / wedgeMass
     wedgeEdgeInvX  = 1.0 / (xmaxWedge - xminWedge)
     wedgeEdgeInvZ  = 1.0 / (zmaxWedge - zminWedge)
!
!
!     ...We are ready to follow the ray's path through all the wedges of the current
!        block. The current wedge indices (i,j) and the previous wedge indices (ip,jp)
!        will be updated as we move through the block. In case a laser IO is performed
!        on the rays, store the initial ray IO data.
!
!
     if (writeRay) then
         nRayWritePositions = nRayWritePositions + 1
         if(nRayWritePositions <= ed_laserIOMaxNumberOfPositions) then
            ed_laserIONumberOfPositions (rayWriteIndex                           ) = nRayWritePositions
            ed_laserIORayPower          (rayWriteIndex, nRayWritePositions       ) = rayPower
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, IAXIS) = rayX
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = rayZ
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = 0.0
         end if
     end if
!
!
!-------------------- Loop following ray through wedges in block --------------------------------------------
!
!
     stepTimeNext = ed_infiniteTime

     do                                ! indefinite loop through the block wedges
                                       ! will be broken (exit) by the various conditions
                                       ! of the ray (no power, out of domain, etc)
!
!
!     ...From the current position, velocity and accelleration of the ray, we determine
!        the initial stepping time to the closest cell wall by assuming constant acceleration.
!        This will be the maximal step time used for the Runge Kutta stepper.
!
!
        x01 = (rayX - xminWedge) * wedgeEdgeInvX                        ! rescaled [0,1] ray x coordinate
        z01 = (rayZ - zminWedge) * wedgeEdgeInvZ                        ! rescaled [0,1] ray z coordinate

        Nele (1:3) = Interpolate_cubic2DFd1 (ed_cellCubicNele (1:16,i,j,1), x01,z01)

        accX = accFactorX * Nele (2)                                    ! acceleration in x-direction
        accZ = accFactorZ * Nele (3)                                    ! acceleration in z-direction

        ymaxWedge = wedgeSlope * rayX
        yminWedge = - ymaxWedge

        tx = ed_time2FacesParabolicPath1D (rayX, velX, accX, xminWedge, xmaxWedge, ed_infiniteTime)
        tz = ed_time2FacesParabolicPath1D (rayZ, velZ, accZ, zminWedge, zmaxWedge, ed_infiniteTime)

        ty = ed_time2FacesParabolicPath1D (rayY,                     &
                                           velY + wedgeSlope * velX, &   ! normal vel to yminWedge
                                                  wedgeSlope * accX, &   ! normal acc to yminWedge
                                           yminWedge,                &   ! just one face here
                                           yminWedge,                &   ! (ymin/ymax not parallel)
                                           ed_infiniteTime           )

        tw = ed_time2FacesParabolicPath1D (rayY,                     &
                                           velY - wedgeSlope * velX, &   ! normal vel to ymaxWedge
                                                - wedgeSlope * accX, &   ! normal acc to ymaxWedge
                                           ymaxWedge,                &   ! just one face here
                                           ymaxWedge,                &   ! (ymin/ymax not parallel)
                                           ed_infiniteTime           )

        stepTimeTry = min (tw, tx, ty, tz, stepTimeNext)                 ! initial time step for RK

        if (stepTimeTry == ed_infiniteTime .or. stepTimeTry == 0.0) then
            call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: infinite/zero stepping time for a cell')
        end if
!
!
!     ...Assemble the vector of dependent variables to be passed to the RK stepper. Since this
!        will be a confined RK step, we also pass the minimum/maximum cell wall coordinates.
!        The RK stepper has to stay within the cell boundaries until one of the cell walls is
!        hit.
!
!        Before calling the RK stepper, set the save computations keyword to true. This signals
!        the ODE function routine to use the current acc(X,Z) and (x,z)01 values for evaluation
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

        rayErrorBase (7) = rayPower * ed_powerStepTolerance   ! power error base depends on current power

        saveComputations = .true.

        call  RungeKutta_stepConfined (ed_RungeKuttaMethod,           &
                                       ed_raytraceODEfunction2DCyl3D, &
                                       3,                             &  ! # of confined variables
                                       0.0,                           &  ! dummy time -> ODE independent of time
                                       rayIn        (1:7),            &
                                       ed_minConfinement2DCyl3D,      &  ! lower limit confinement function
                                       ed_maxConfinement2DCyl3D,      &  ! upper limit confinement function
                                       rayErrorFrac,                  &
                                       rayErrorBase (1:7),            &
                                       stepTimeTry,                   &
                                       stepTimeUsed,                  &  ! actual stepping time used
                                       stepTimeNext,                  &  ! next stepping time
                                       rayOut       (1:7),            &
                                       rayError     (1:7)             )  ! can be +ve or -ve
!
!
!     ...The confined RK step has been taken. Check, if the errors are within the error bars and check
!        where the ray is currently located. Add the lost ray power to the cell energy and update the
!        (diminished) power of the ray.
!
!
        if (any (abs (rayError (1:7)) > rayErrorFrac * rayErrorBase (1:7))) then
            call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: error(s) in RK step out of bounds!')
        end if

        rayX = rayOut (1)
        rayY = rayOut (2)
        rayZ = rayOut (3)
        velX = rayOut (4)
        velY = rayOut (5)
        velZ = rayOut (6)

        ymaxWedge = wedgeSlope * rayX
        yminWedge = - ymaxWedge

        outOfWedge =     (rayX < xminWedge - wedgeWallThicknessHalf) &    ! for debugging purposes
                    .or. (rayX > xmaxWedge + wedgeWallThicknessHalf) &    ! will be removed once the
                    .or. (rayY < yminWedge - wedgeWallThicknessHalf) &    ! code (the confined RK stepper)
                    .or. (rayY > ymaxWedge + wedgeWallThicknessHalf) &    ! is running properly
                    .or. (rayZ < zminWedge - wedgeWallThicknessHalf) &
                    .or. (rayZ > zmaxWedge + wedgeWallThicknessHalf)

        if (outOfWedge) then
            call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: RK stepped out of wedge confinement!')
        end if

        newWedge =     (rayX < xminWedge + wedgeWallThicknessHalf) &      ! checks, if the current
                  .or. (rayX > xmaxWedge - wedgeWallThicknessHalf) &      ! confined RK step has hit
                  .or. (rayY < yminWedge + wedgeWallThicknessHalf) &      ! one of the cell walls
                  .or. (rayY > ymaxWedge - wedgeWallThicknessHalf) &
                  .or. (rayZ < zminWedge + wedgeWallThicknessHalf) &
                  .or. (rayZ > zmaxWedge - wedgeWallThicknessHalf)

        wedgePower  = rayPower - rayOut (7)
        wedgeEnergy = wedgePower * timeStep

        if (ed_depoVarIsPerMass) then
            wedgeEnergyDepot (i,j) = wedgeEnergyDepot (i,j) + wedgeEnergy * wedgeMassInv
        else
            wedgeEnergyDepot (i,j) = wedgeEnergyDepot (i,j) + wedgeEnergy * wedgeVolumeInv
        end if

        rayPower = rayOut (7)

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
        velXeq0 = (velX == 0.0)
        velYeq0 = (velY == 0.0)
        velZeq0 = (velZ == 0.0)

        stationaryRay  = velXeq0 .and. velYeq0 .and. velZeq0

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
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = rayZ
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = 0.0
            end if
        end if
!
!
!     ...If, at the current stage, the ray enters a new wedge, we have to determine: 1) which new
!        wedge (i,j) it is and 2) the appropriate nudging values on the ray's position. Due to
!        possible reflective boundary conditions on the block faces, it can happen that the ray
!        stays in the original wedge. After handling the logistics inside the following 'if'
!        statement, the new wedge indices i,j are either the old ones or new ones. The keyword
!        'newWedgeIJ' will indicate, if the indices i,j will change.
!
!
        newWedgeIJ = .false.

        if (newWedge) then

            ymaxWedge = + wedgeSlope * rayX
            yminWedge = - ymaxWedge

            dist2minX = abs (xminWedge - rayX)
            dist2maxX = abs (xmaxWedge - rayX)
            dist2minY = abs (yminWedge - rayY)
            dist2maxY = abs (ymaxWedge - rayY)
            dist2minZ = abs (zminWedge - rayZ)
            dist2maxZ = abs (zmaxWedge - rayZ)

            minDistance = min (dist2minX, dist2maxX, &
                               dist2minY, dist2maxY, &
                               dist2minZ, dist2maxZ)

            if (minDistance > wedgeWallThicknessHalf) then
                call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: ray too far away from wedge face')
            end if

            wedgeFaceMinX = (dist2minX <= wedgeWallThicknessHalf)
            wedgeFaceMaxX = (dist2maxX <= wedgeWallThicknessHalf)
            wedgeFaceMinY = (dist2minY <= wedgeWallThicknessHalf) .and. (rayY < 0.0)
            wedgeFaceMaxY = (dist2maxY <= wedgeWallThicknessHalf) .and. (rayY > 0.0)
            wedgeFaceMinZ = (dist2minZ <= wedgeWallThicknessHalf)
            wedgeFaceMaxZ = (dist2maxZ <= wedgeWallThicknessHalf)

            velXgt0 = (velX  > 0.0)
            velXlt0 = (velX  < 0.0)
            velZgt0 = (velZ  > 0.0)
            velZlt0 = (velZ  < 0.0)

            crossX = .false.
            crossZ = .false.

            nudgeX = 0.0
            nudgeZ = 0.0

            ip = i
            jp = j

            if (wedgeFaceMinX) then

                rayX   = xminWedge
                nudgeX = + wedgeWallThicknessHalf

                if (velXlt0) then
                    i = i - 1
                    crossX = .true.
                end if

            else if (wedgeFaceMaxX) then

                rayX   = xmaxWedge
                nudgeX = - wedgeWallThicknessHalf

                if (velXgt0) then
                    i = i + 1
                    crossX = .true.
                end if

            end if

            if (wedgeFaceMinY) then

                if (wedgeSlope * velX < - velY) then
                    rayY = ymaxWedge
                    velX = velX * wedgeCosine - velY * wedgeSine
                    velY = velX * wedgeSine   + velY * wedgeCosine
                else
                    rayY = yminWedge
                end if

            else if (wedgeFaceMaxY) then

                if (wedgeSlope * velX < velY) then
                    rayY = yminWedge
                    velX = velY * wedgeSine   + velX * wedgeCosine
                    velY = velY * wedgeCosine - velX * wedgeSine
                else
                    rayY = ymaxWedge
                end if

            end if

            if (wedgeFaceMinZ) then

                rayZ   = zminWedge
                nudgeZ = + wedgeWallThicknessHalf

                if (velZlt0) then
                    j = j - 1
                    crossZ = .true.
                end if

            else if (wedgeFaceMaxZ) then

                rayZ   = zmaxWedge
                nudgeZ = - wedgeWallThicknessHalf

                if (velZgt0) then
                    j = j + 1
                    crossZ = .true.
                end if

            end if

            blockFaceMinX = (rayX == xminBlock)
            blockFaceMaxX = (rayX == xmaxBlock)
            blockFaceMinZ = (rayZ == zminBlock)
            blockFaceMaxZ = (rayZ == zmaxBlock)

            reflectX =     (blockFaceMinX .and. blockReflectMinX .and. velXlt0) &
                      .or. (blockFaceMaxX .and. blockReflectMaxX .and. velXgt0)
            reflectZ =     (blockFaceMinZ .and. blockReflectMinZ .and. velZlt0) &
                      .or. (blockFaceMaxZ .and. blockReflectMaxZ .and. velZgt0)

            if (reflectX) then
                i = ip
                velX = - velX
                crossX = .false.
            end if

            if (reflectZ) then
                j = jp
                velZ = - velZ
                crossZ = .false.
            end if

            if (crossX) then
                nudgeX = (i - ip) * wedgeWallThicknessHalf
            end if

            if (crossZ) then
                nudgeZ = (j - jp) * wedgeWallThicknessHalf
            end if

            rayX = rayX + nudgeX
            rayZ = rayZ + nudgeZ

            newWedgeIJ = crossX .or. crossZ

        end if
!
!
!     ...We are now sure about the target wedge. Check, if the target wedge (i,j) is still within the block.
!        If it is, we check if this is a wedge with new i,j indices, in which case we update the i,j info.
!        If the target wedge is not within the block, check if the ray coordinates are still within the
!        defined domain. If not, store its latest data and mark it as nonexistent. If the ray is still
!        within the domain boundaries, exit the current block loop.
!
!
        inBlock =      (i >= iminBlock) &
                 .and. (i <= imaxBlock) &
                 .and. (j >= jminBlock) &
                 .and. (j <= jmaxBlock)

        if (inBlock) then

            if (newWedgeIJ) then

                xminWedge    = ed_cellEdges (i  ,1)
                xmaxWedge    = ed_cellEdges (i+1,1)
                zminWedge    = ed_cellEdges (j  ,2)
                zmaxWedge    = ed_cellEdges (j+1,2)

                wedgeZbar      = ed_cellZbar    (i,j,1)
                wedgeDensity   = ed_cellDensity (i,j,1)
                wedgeVolume    = ed_cellVolume  (i,j,1)
                wedgeVolumeInv = 1.0 / wedgeVolume
                wedgeMass      = wedgeDensity * wedgeVolume
                wedgeMassInv   = 1.0 / wedgeMass
                wedgeEdgeInvX  = 1.0 / (xmaxWedge - xminWedge)
                wedgeEdgeInvZ  = 1.0 / (zmaxWedge - zminWedge)

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
                      .and. (rayZ > ed_yminDomain) &
                      .and. (rayZ < ed_ymaxDomain)

            if (.not.inDomain) then
                 ed_rays (RAY_POSX,ray) = rayX - nudgeX             ! undo the nudging, since
                 ed_rays (RAY_POSZ,ray) = rayZ - nudgeZ             ! it is not needed anymore
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
!-------------------- End loop following ray through wedges in block --------------------------------------------
!
!
     end do
!
!
!     ...Check to see if we ran out of laser IO buffer space.
!
!
     if (writeRay .and. (nRayWritePositions > ed_laserIOMaxNumberOfPositions) ) then
         print *, "[ed_traceBlockRays2DCyl3D] Ray ", ray, &
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
end subroutine ed_traceBlockRays2DCyl3D
