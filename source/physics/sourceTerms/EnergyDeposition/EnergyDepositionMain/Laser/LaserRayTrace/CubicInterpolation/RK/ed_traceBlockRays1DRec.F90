!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/RK/ed_traceBlockRays1DRec
!!
!! NAME
!!
!!  ed_traceBlockRays1DRec
!!
!! SYNOPSIS
!!
!!  call ed_traceBlockRays1DRec (real    (in)    :: timeStep,
!!                               integer (in)    :: rayFirst
!!                               integer (in)    :: rayLast,
!!                               integer (in)    :: iminBlock,
!!                               integer (in)    :: imaxBlock,
!!                               real    (in)    :: xminBlock,
!!                               real    (in)    :: xmaxBlock,
!!                               real    (in)    :: deltaX,
!!                               real    (in)    :: deltaInvX,
!!                               logical (in)    :: blockReflectMinX,
!!                               logical (in)    :: blockReflectMaxX,
!!                               real    (inout) :: cellEnergyDepot (:))
!!
!! DESCRIPTION
!!
!!  Traces the movement of the current collection of active rays through one block for
!!  those geometries consisting formally of 1D rectangular grids (cartesian + spherical).
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
!!  xminBlock        : minimum x-coordinate limit of the block
!!  xmaxBlock        : maximum x-coordinate limit of the block
!!  deltaX           : the cell's x-dimension
!!  deltaInvX        : inverse of the cell's x-dimension
!!  blockReflectMinX : is the block boundary on the minimum x-side reflective ?
!!  blockReflectMaxX : is the block boundary on the maximum x-side reflective ?
!!  cellEnergyDepot  : array collecting the ray energy deposition for each cell
!!
!! NOTES
!!
!!  The code allows for threading to be used on the outer ray trace loop.
!!  The paths of the rays are computed using the geometric optics approximation. 
!!  Inside each cell, the paths of the rays are evaluated stepwise, using the
!!  monocubic expansions of the number electron density and electron temperature
!!  grid.
!!
!!***

subroutine ed_traceBlockRays1DRec (timeStep,                          &
                                   rayFirst,  rayLast,                &
                                   iminBlock, imaxBlock,              &
                                   xminBlock, xmaxBlock,              &
                                   deltaX,                            &
                                   deltaInvX,                         &
                                   blockReflectMinX,                  &
                                   blockReflectMaxX,                  &
                                                      cellEnergyDepot ) 

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
                                         ed_xmaxDomain

  use ed_raytraceODEfunctionData, ONLY : accX,             &
                                         accFactorX,       &
                                         cellEdgeInvX,     &
                                         cellZbar,         &
                                         i,                &
                                         Nele,             &
                                         rayCritDens,      &
                                         saveComputations, &
                                         x01,              &
                                         xmaxCell,         &
                                         xminCell

  use Driver_interface,           ONLY : Driver_abortFlash

  use Interpolate_interface,      ONLY : Interpolate_cubic1DF,  &
                                         Interpolate_cubic1DFd1

  use RungeKutta_interface,       ONLY : RungeKutta_stepConfined

  use ed_interface,               ONLY : ed_maxConfinement1DRec,       &
                                         ed_minConfinement1DRec,       &
                                         ed_raytraceODEfunction1Drec,  &
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
  real,    intent (in)    :: xminBlock, xmaxBlock
  real,    intent (in)    :: deltaX
  real,    intent (in)    :: deltaInvX
  logical, intent (in)    :: blockReflectMinX
  logical, intent (in)    :: blockReflectMaxX
  real,    intent (inout) :: cellEnergyDepot (iminBlock:imaxBlock)

  logical :: badTimeStep
  logical :: blockFaceMinX, blockFaceMaxX
  logical :: cellFaceMinX, cellFaceMaxX
  logical :: crossX
  logical :: inDomain, inBlock
  logical :: newCell, outOfCell
  logical :: onBlockBoundaryCell
  logical :: rayOutOfBlock
  logical :: reflectX
  logical :: stationaryRay
  logical :: velXeq0, velXgt0, velXlt0
  logical :: writeRay

  integer :: ip
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
  real    :: cellStepErrorX
  real    :: cellVolume, cellVolumeInv
  real    :: cellWallThicknessHalf
  real    :: dist2minX, dist2maxX
  real    :: minDistance
  real    :: nudgeX
  real    :: rayErrorFrac
  real    :: rayPower
  real    :: rayX
  real    :: stepTimeTry, stepTimeUsed, stepTimeNext
  real    :: tx
  real    :: velX

  real    :: rayError     (1:3)
  real    :: rayErrorBase (1:3)
  real    :: rayIn        (1:3)
  real    :: rayOut       (1:3)
!
!
!     ...Define some variables. The error base for the ray power has to be set
!        inside the deposition loop, as it depends on the current ray power.
!
!
  numDeadRays = 0

  cellWallThicknessHalf = 0.5 * ed_cellWallThickness
  cellStepErrorX        = ed_cellStepTolerance * deltaX

  rayErrorFrac = 1.0                       ! this means the error is controlled by the base only

  rayErrorBase (1) = cellStepErrorX        ! error bar on rayX
  rayErrorBase (2) = ed_infiniteSpeed      ! this in effect puts no error bars on velX yet
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
     velX        =      ed_rays (RAY_VELX,ray)
     rayPower    =      ed_rays (RAY_POWR,ray)
     rayCritDens =      ed_rays (RAY_DENC,ray)

     c2div2nc = ed_speedOfLightSquared / (rayCritDens + rayCritDens)

     accFactorX = - c2div2nc * deltaInvX
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
!     ...Find the index (i) of the initial cell through which the ray will
!        enter the block. We know for sure that the ray enters the block, because
!        otherwise it would not be on the current block list. Check, on which of
!        the two possible faces the ray currently is.
!
!
     rayOutOfBlock =  (rayX < xminBlock) .or. (rayX > xmaxBlock)

     if (rayOutOfBlock) then
         call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: ray found out of block')
     end if

     dist2minX = rayX - xminBlock
     dist2maxX = xmaxBlock - rayX

     minDistance = min (dist2minX, dist2maxX)

     if (minDistance > ed_cellWallThickness) then
         call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: ray too far inside the block')
     end if

     i = iminBlock + int ( (rayX - xminBlock) * deltaInvX )

     onBlockBoundaryCell = (i == iminBlock) .or. (i == imaxBlock)

     if (.not.onBlockBoundaryCell) then
          call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: ray not in a block boundary cell')
     end if

     xminCell = ed_cellEdges (i  ,1)
     xmaxCell = ed_cellEdges (i+1,1)

     dist2minX = abs (xminCell - rayX)
     dist2maxX = abs (xmaxCell - rayX)

     cellFaceMinX = (dist2minX <= cellWallThicknessHalf)
     cellFaceMaxX = (dist2maxX <= cellWallThicknessHalf)
!
!
!     ...Make sure the ray is also properly nudged into the corresponding cell.
!
!
     if (cellFaceMinX) rayX = xminCell + cellWallThicknessHalf
     if (cellFaceMaxX) rayX = xmaxCell - cellWallThicknessHalf

     velXeq0 = (velX == 0.0)

     stationaryRay = velXeq0

     if (stationaryRay) then
         call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: stationary ray at a block face boundary')
     end if
!
!
!     ...Get extra needed info about the initial cell (i). The calculation of the inverse
!        length of the cell edges have to be based on the actual min/max cell values. If
!        based on the (global per block) cell delta values, roundoff errors introduced during
!        evaluation of the min/max cell values can lead to abortion of the cubic interpolation
!        routines, which strictly! require [0,1] rescaled coordinates.
!
!
     cellZbar      = ed_cellZbar    (i,1,1)
     cellDensity   = ed_cellDensity (i,1,1)
     cellVolume    = ed_cellVolume  (i,1,1)
     cellVolumeInv = 1.0 / cellVolume
     cellMass      = cellDensity * cellVolume
     cellMassInv   = 1.0 / cellMass
     cellEdgeInvX  = 1.0 / (xmaxCell - xminCell)
!
!
!     ...We are ready to follow the ray's path through all the cells of the current
!        block. The current cell index (i) and the previous cell index (ip) will be
!        updated as we move through the block. In case a laser IO is performed on the
!        rays, store the initial ray IO data.
!
!
     if (writeRay) then
         nRayWritePositions = nRayWritePositions + 1
         if(nRayWritePositions <= ed_laserIOMaxNumberOfPositions) then
            ed_laserIONumberOfPositions (rayWriteIndex                           ) = nRayWritePositions
            ed_laserIORayPower          (rayWriteIndex, nRayWritePositions       ) = rayPower
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, IAXIS) = rayX
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = 0.0
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = 0.0
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

        Nele (1:2) = Interpolate_cubic1DFd1 (ed_cellCubicNele (1:4,i,1,1), x01)

        accX = accFactorX * Nele (2)                                    ! acceleration in x-direction

        tx = ed_time2FacesParabolicPath1D (rayX, velX, accX, xminCell, xmaxCell, ed_infiniteTime)

        stepTimeTry = min (tx, stepTimeNext)                            ! initial time step for RK

        if (stepTimeTry == ed_infiniteTime .or. stepTimeTry == 0.0) then
            call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: infinite/zero stepping time for a cell')
        end if
!
!
!     ...Assemble the vector of dependent variables to be passed to the RK stepper. Since this
!        will be a confined RK step, we also pass the minimum/maximum cell wall coordinates.
!        The RK stepper has to stay within the cell boundaries until one of the cell walls is
!        hit.
!
!        Before calling the RK stepper, set the save computations keyword to true. This signals
!        the ODE function routine to use the current accX and x01 values for evaluation
!        of the first RK point. It avoids the rather costly recalculation of the Nele values and
!        their derivatives at the current ray location. The keyword is set imediately to false
!        after that (inside the ODE function routine), because the RK stepper saves an internal
!        copy of the first ODE function vector.
!
!
        rayIn (1) = rayX
        rayIn (2) = velX
        rayIn (3) = rayPower

        rayErrorBase (3) = rayPower * ed_powerStepTolerance   ! power error base depends on current power

        saveComputations = .true.

        call  RungeKutta_stepConfined (ed_RungeKuttaMethod,         &
                                       ed_raytraceODEfunction1Drec, &
                                       1,                           &  ! # of confined variables
                                       0.0,                         &  ! dummy time -> ODE independent of time
                                       rayIn        (1:3),          &
                                       ed_minConfinement1DRec,      &  ! lower limit confinement function
                                       ed_maxConfinement1DRec,      &  ! upper limit confinement function
                                       rayErrorFrac,                &
                                       rayErrorBase (1:3),          &
                                       stepTimeTry,                 &
                                       stepTimeUsed,                &  ! actual stepping time used
                                       stepTimeNext,                &  ! estimates of next stepping time
                                       rayOut       (1:3),          &
                                       rayError     (1:3)           )  ! can be +ve or -ve
!
!
!     ...The confined RK step has been taken. Check, if the errors are within the error bars and check
!        where the ray is currently located. Add the lost ray power to the cell energy and update the
!        (diminished) power of the ray.
!
!
        if (any (abs (rayError (1:3)) > rayErrorFrac * rayErrorBase (1:3))) then
            call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: error(s) in RK step out of bounds!')
        end if

        rayX = rayOut (1)
        velX = rayOut (2)

        outOfCell =     (rayX < xminCell - cellWallThicknessHalf) &    ! for debugging purposes
                   .or. (rayX > xmaxCell + cellWallThicknessHalf)      ! (will be removed later)

        if (outOfCell) then
            call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: RK stepped out of cell confinement!')
        end if

        newCell =     (rayX < xminCell + cellWallThicknessHalf) &      ! checks, if the current confined RK
                 .or. (rayX > xmaxCell - cellWallThicknessHalf)        ! step has hit one of the cell walls

        cellPower  = rayPower - rayOut (3)
        cellEnergy = cellPower * timeStep

        if (ed_depoVarIsPerMass) then
            cellEnergyDepot (i) = cellEnergyDepot (i) + cellEnergy * cellMassInv
        else
            cellEnergyDepot (i) = cellEnergyDepot (i) + cellEnergy * cellVolumeInv
        end if

        rayPower = rayOut (3)

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
 
        stationaryRay = velXeq0

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
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = 0.0
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = 0.0
            end if
        end if
!
!
!     ...If, at the current stage, the ray enters a new cell, we have to determine: 1) which new
!        cell (i) it is and 2) the appropriate nudging value on the ray's position. Due to
!        possible reflective boundary conditions on the block faces, it can happen that the ray
!        stays in the original cell. After handling the logistics inside the following 'if'
!        statement, the new cell index i is either the old one or a new one.
!
!
        if (newCell) then

            dist2minX = abs (xminCell - rayX)
            dist2maxX = abs (xmaxCell - rayX)

            minDistance = min (dist2minX, dist2maxX)

            if (minDistance > cellWallThicknessHalf) then
                call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: ray too far away from cell face')
            end if

            cellFaceMinX = (dist2minX <= cellWallThicknessHalf)
            cellFaceMaxX = (dist2maxX <= cellWallThicknessHalf)

            velXgt0 = (velX  > 0.0)
            velXlt0 = (velX  < 0.0)

            crossX = .false.
            nudgeX = 0.0

            ip = i

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

            blockFaceMinX = (rayX == xminBlock)
            blockFaceMaxX = (rayX == xmaxBlock)

            reflectX =     (blockFaceMinX .and. blockReflectMinX .and. velXlt0) &
                      .or. (blockFaceMaxX .and. blockReflectMaxX .and. velXgt0)

            if (reflectX) then
                i = ip
                velX = - velX
                crossX = .false.
            end if

            if (crossX) then
                nudgeX = (i - ip) * cellWallThicknessHalf
            end if

            rayX = rayX + nudgeX

            newCell = crossX

        end if
!
!
!     ...We are now sure about the target cell. Check, if the target cell (i) is still within the block.
!        If it is, we check if this is a new cell, in which case we update the cell info. If the target
!        cell is not within the block, check if the ray coordinates are still within the defined domain.
!        If not, store its latest data and mark it as nonexistent. If the ray is still within the domain
!        boundaries, exit the current block loop.
!
!
        inBlock = (i >= iminBlock) .and. (i <= imaxBlock)

        if (inBlock) then

            if (newCell) then

                xminCell = ed_cellEdges   (i  ,1)
                xmaxCell = ed_cellEdges   (i+1,1)

                cellZbar      = ed_cellZbar    (i,1,1)
                cellDensity   = ed_cellDensity (i,1,1)
                cellVolume    = ed_cellVolume  (i,1,1)
                cellVolumeInv = 1.0 / cellVolume
                cellMass      = cellDensity * cellVolume
                cellMassInv   = 1.0 / cellMass
                cellEdgeInvX  = 1.0 / (xmaxCell - xminCell)

            end if

        else

            ed_rays (RAY_POSX,ray) = rayX
            ed_rays (RAY_VELX,ray) = velX
            ed_rays (RAY_POWR,ray) = rayPower

            inDomain = (rayX > ed_xminDomain) .and. (rayX < ed_xmaxDomain)

            if (.not.inDomain) then
                 ed_rays (RAY_POSX,ray) = rayX - nudgeX          ! undo the nudging (not needed anymore)
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
        print *, "[ed_traceBlockRays1DRec] Ray ", ray, &
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
end subroutine ed_traceBlockRays1DRec
