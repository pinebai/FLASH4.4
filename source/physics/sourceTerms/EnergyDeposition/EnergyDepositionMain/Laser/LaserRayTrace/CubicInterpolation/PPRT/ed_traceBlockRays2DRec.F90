!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/PPRT/ed_traceBlockRays2DRec
!!
!! NAME
!!
!!  ed_traceBlockRays2DRec
!!
!! SYNOPSIS
!!
!!  call ed_traceBlockRays2DRec (real    (in)    :: timeStep,
!!                               integer (in)    :: rayFirst
!!                               integer (in)    :: rayLast,
!!                               integer (in)    :: iminBlock,
!!                               integer (in)    :: imaxBlock,
!!                               integer (in)    :: jminBlock,
!!                               integer (in)    :: jmaxBlock,
!!                               real    (in)    :: xminBlock,
!!                               real    (in)    :: xmaxBlock,
!!                               real    (in)    :: yminBlock,
!!                               real    (in)    :: ymaxBlock,
!!                               real    (in)    :: deltaX,
!!                               real    (in)    :: deltaY,
!!                               real    (in)    :: deltaInvX,
!!                               real    (in)    :: deltaInvY,
!!                               logical (in)    :: blockReflectMinX,
!!                               logical (in)    :: blockReflectMaxX,
!!                               logical (in)    :: blockReflectMinY,
!!                               logical (in)    :: blockReflectMaxY,
!!                               real    (inout) :: cellEnergyDepot (:,:))
!!
!! DESCRIPTION
!!
!!  Traces the movement of the current collection of active rays through one block for
!!  those geometries consisting formally of 2D rectangular grids (cartesian + cylindrical).
!!  On exit, each ray has either:
!!
!!            i)  reached a different (yet unknown) block
!!           ii)  has been absorbed -> is nonexistent
!!          iii)  has reached the domain boundary and exited -> is nonexistent.
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
!!  xminBlock        : minimum x-coordinate limit of the block
!!  xmaxBlock        : maximum x-coordinate limit of the block
!!  yminBlock        : minimum y-coordinate limit of the block
!!  ymaxBlock        : maximum y-coordinate limit of the block
!!  deltaX           : the cell's x-dimension
!!  deltaY           : the cell's y-dimension
!!  deltaInvX        : inverse of the cell's x-dimension
!!  deltaInvY        : inverse of the cell's y-dimension
!!  blockReflectMinX : is the block boundary on the minimum x-side reflective ?
!!  blockReflectMaxX : is the block boundary on the maximum x-side reflective ?
!!  blockReflectMinY : is the block boundary on the minimum y-side reflective ?
!!  blockReflectMaxY : is the block boundary on the maximum y-side reflective ?
!!  cellEnergyDepot  : array collecting the ray energy deposition for each cell
!!
!! NOTES
!!        
!!  The code allows for threading to be used on the outer ray trace loop.
!!  The paths of the rays are computed using the geometric optics approximation.
!!  Inside each cell, the paths of the rays are evaluated stepwise, using the
!!  bicubic expansions of the number electron density and electron temperature
!!  grid.
!!
!!***

subroutine ed_traceBlockRays2DRec (timeStep,                          &
                                   rayFirst,  rayLast,                &
                                   iminBlock, imaxBlock,              &
                                   jminBlock, jmaxBlock,              &
                                   xminBlock, xmaxBlock,              &
                                   yminBlock, ymaxBlock,              &
                                   deltaX,    deltaY,                 &
                                   deltaInvX, deltaInvY,              &
                                   blockReflectMinX,                  &
                                   blockReflectMaxX,                  &
                                   blockReflectMinY,                  &
                                   blockReflectMaxY,                  &
                                                      cellEnergyDepot ) 

  use EnergyDeposition_data,  ONLY : ed_Boltzmann,                   &
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
                                     ed_rays,                        &
                                     ed_raysMovedIntoDomain,         &
                                     ed_rayZeroPower,                &
                                     ed_speedOfLightSquared,         &
                                     ed_unitRoundoff,                &
                                     ed_xminDomain,                  &
                                     ed_xmaxDomain,                  &
                                     ed_yminDomain,                  &
                                     ed_ymaxDomain

  use Driver_interface,       ONLY : Driver_abortFlash

  use Interpolate_interface,  ONLY : Interpolate_cubic2DF,  &
                                     Interpolate_cubic2DFd1

  use ed_interface,           ONLY : ed_CoulombFactor,             &
                                     ed_inverseBremsstrahlungRate, &
                                     ed_time2FacesParabolicPath1D

  use ed_commInterface,       ONLY : ed_commHandleOffBlkRay,   &
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
  real,    intent (in)    :: yminBlock, ymaxBlock
  real,    intent (in)    :: deltaX,    deltaY
  real,    intent (in)    :: deltaInvX, deltaInvY
  logical, intent (in)    :: blockReflectMinX
  logical, intent (in)    :: blockReflectMaxX
  logical, intent (in)    :: blockReflectMinY
  logical, intent (in)    :: blockReflectMaxY
  real,    intent (inout) :: cellEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock)

  logical :: badTimeStep
  logical :: blockFaceMinX, blockFaceMaxX
  logical :: blockFaceMinY, blockFaceMaxY
  logical :: cellFaceMinX,  cellFaceMaxX
  logical :: cellFaceMinY,  cellFaceMaxY
  logical :: crossX, crossY
  logical :: inDomain, inBlock
  logical :: newCell
  logical :: onBlockBoundaryCell
  logical :: rayOutOfBlock
  logical :: reflectX, reflectY
  logical :: stationaryRay
  logical :: velXeq0, velXgt0, velXlt0
  logical :: velYeq0, velYgt0, velYlt0
  logical :: writeRay

  integer :: i,j,n
  integer :: ip,jp
  integer :: nRayWritePositions
  integer :: numDeadRays
  integer :: ray
  integer :: rayTag
  integer :: rayWriteIndex

  real    :: a1X, a1Y
  real    :: accX, accY
  real    :: c2div2nc
  real    :: cellDensity
  real    :: cellEnergy
  real    :: cellMass, cellMassInv
  real    :: cellPower
  real    :: cellStepX, cellStepY
  real    :: cellVolume, cellVolumeInv
  real    :: cellWallThicknessHalf
  real    :: cellZbar
  real    :: dist2minX, dist2minY
  real    :: dist2maxX, dist2maxY
  real    :: dX, dY
  real    :: gaussTime, gaussTimeHalf
  real    :: integral
  real    :: lnLambda
  real    :: minDistance
  real    :: Nele, Tele
  real    :: nu
  real    :: nudgeX, nudgeY
  real    :: powerLossFactor
  real    :: r0X, r0Y, r1X, r1Y, r2X, r2Y
  real    :: rayCritDens
  real    :: rayPower
  real    :: rayX, rayY
  real    :: stepTime, stepTimeHalf, stepTimeFourth
  real    :: v0X, v0Y, v1X, v1Y
  real    :: velX, velY
  real    :: x, y
  real    :: xminCell, yminCell
  real    :: xmaxCell, ymaxCell

  real    :: NeleDerv (1:3)

  real, parameter :: GaussianRoot (1:2) = (/1.577350269189626,    &  ! is 1 + 1/sqrt(3) for int limits [0,1]
                                            4.226497308103742e-01/)  ! is 1 - 1/sqrt(3) for int limits [0,1]
!
!
!     ...Define some variables.
!
!
  numDeadRays = 0

  cellWallThicknessHalf = 0.5 * ed_cellWallThickness
  cellStepX             = ed_cellStepTolerance * deltaX
  cellStepY             = ed_cellStepTolerance * deltaY
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
     velX        =      ed_rays (RAY_VELX,ray)
     velY        =      ed_rays (RAY_VELY,ray)
     rayPower    =      ed_rays (RAY_POWR,ray)
     rayCritDens =      ed_rays (RAY_DENC,ray)

     c2div2nc = ed_speedOfLightSquared / (rayCritDens + rayCritDens)
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
!     ...Find the indices (i,j) of the initial cell through which the ray will
!        enter the block. We know for sure that the ray enters the block, because
!        otherwise it would not be on the current block list. Check, on which of
!        the four possible faces the ray currently is.
!
!
     rayOutOfBlock =     (rayX < xminBlock) &
                    .or. (rayX > xmaxBlock) &
                    .or. (rayY < yminBlock) &
                    .or. (rayY > ymaxBlock)

     if (rayOutOfBlock) then
         call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: ray found out of block')
     end if

     dist2minX = rayX - xminBlock
     dist2maxX = xmaxBlock - rayX
     dist2minY = rayY - yminBlock
     dist2maxY = ymaxBlock - rayY

     minDistance = min (dist2minX, dist2maxX, &
                        dist2minY, dist2maxY)

     if (minDistance > ed_cellWallThickness) then
         call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: ray too far inside the block')
     end if

     i = iminBlock + int ( (rayX - xminBlock) * deltaInvX )
     j = jminBlock + int ( (rayY - yminBlock) * deltaInvY )

     onBlockBoundaryCell = (     (i == iminBlock) &
                            .or. (i == imaxBlock) &
                            .or. (j == jminBlock) &
                            .or. (j == jmaxBlock) )

     if (.not.onBlockBoundaryCell) then
          call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: ray not in a block boundary cell')
     end if

     xminCell = ed_cellEdges (i  ,1)
     xmaxCell = ed_cellEdges (i+1,1)
     yminCell = ed_cellEdges (j  ,2)
     ymaxCell = ed_cellEdges (j+1,2)

     dist2minX = abs (xminCell - rayX)
     dist2maxX = abs (xmaxCell - rayX)
     dist2minY = abs (yminCell - rayY)
     dist2maxY = abs (ymaxCell - rayY)

     cellFaceMinX = (dist2minX <= cellWallThicknessHalf)
     cellFaceMaxX = (dist2maxX <= cellWallThicknessHalf)
     cellFaceMinY = (dist2minY <= cellWallThicknessHalf)
     cellFaceMaxY = (dist2maxY <= cellWallThicknessHalf)
!
!
!     ...Make sure the ray is also properly nudged into the corresponding cell.
!
!
     if (cellFaceMinX) rayX = xminCell + cellWallThicknessHalf
     if (cellFaceMaxX) rayX = xmaxCell - cellWallThicknessHalf
     if (cellFaceMinY) rayY = yminCell + cellWallThicknessHalf
     if (cellFaceMaxY) rayY = ymaxCell - cellWallThicknessHalf

     velXeq0 = (velX == 0.0)
     velYeq0 = (velY == 0.0)

     stationaryRay = (velXeq0 .and. velYeq0)

     if (stationaryRay) then
         call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: stationary ray at a block face boundary')
     end if
!
!
!     ...Get extra needed info about the initial cell (i,j).
!
!
     cellZbar      = ed_cellZbar    (i,j,1)
     cellDensity   = ed_cellDensity (i,j,1)
     cellVolume    = ed_cellVolume  (i,j,1)
     cellVolumeInv = 1.0 / cellVolume
     cellMass      = cellDensity * cellVolume
     cellMassInv   = 1.0 / cellMass
!
!
!     ...We are ready to follow the ray's path through all the cells of the current
!        block. At this point, the ray is in the initial cell (i,j) and is ready to move
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
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = 0.0
         end if
     end if
!
!
!-------------------- Loop following ray through cells in block --------------------------------------------
!
!
     do                                ! indefinite loop through the block cells
                                       ! will be broken (exit) by the various conditions
                                       ! of the ray (no power, out of domain, etc)
!
!
!     ...From the current position, velocity and accelleration of the ray, we determine
!        the stepping time, such that the ray's parabolic path is accurate to within
!        the demanded position error tolerance. After this section, we know the exact
!        analytic parabolic ray path during the stepping time. We also know that this
!        path is entirely within the current cell.
!
!
        x = (rayX - xminCell) * deltaInvX                               ! rescaled [0,1] ray x coordinate
        y = (rayY - yminCell) * deltaInvY                               ! rescaled [0,1] ray y coordinate

        NeleDerv (1:3) = Interpolate_cubic2DFd1 (ed_cellCubicNele (1:16,i,j,1), x,y)

        dX  = NeleDerv (2) * deltaInvX                                  ! d/dx -> d/dX
        dY  = NeleDerv (3) * deltaInvY                                  ! d/dy -> d/dY

        accX = - c2div2nc * dX                                          ! acceleration in x-direction
        accY = - c2div2nc * dY                                          ! acceleration in y-direction

        x = ed_time2FacesParabolicPath1D (rayX, velX, accX, xminCell, xmaxCell, ed_infiniteTime)
        y = ed_time2FacesParabolicPath1D (rayY, velY, accY, yminCell, ymaxCell, ed_infiniteTime)

        stepTime = min (x,y)                                            ! time to reach cell boundary

        if (stepTime == ed_infiniteTime .or. stepTime == 0.0) then
            call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: infinite/zero crossing time for a cell')
        end if

        stepTimeHalf = 0.5  * stepTime

        r0X = rayX + (velX + accX * stepTimeHalf) * stepTime            ! full time step x-position
        r0Y = rayY + (velY + accY * stepTimeHalf) * stepTime            ! full time step y-position

        v0X = velX + accX * stepTime                                    ! full time step x-velocity
        v0Y = velY + accY * stepTime                                    ! full time step y-velocity

        newCell = .true.

        do                                                              ! check ray path inside cell

           stepTimeFourth = 0.5 * stepTimeHalf

           r1X = rayX + (velX + accX * stepTimeFourth) * stepTimeHalf   ! 1st half time step x-position
           r1Y = rayY + (velY + accY * stepTimeFourth) * stepTimeHalf   ! 1st half time step y-position

           v1X = velX + accX * stepTimeHalf                             ! 1st half time step x-velocity
           v1Y = velY + accY * stepTimeHalf                             ! 1st half time step y-velocity

           x = (r1X - xminCell) * deltaInvX                             ! rescaled [0,1] ray x coordinate
           y = (r1Y - yminCell) * deltaInvY                             ! rescaled [0,1] ray y coordinate

           NeleDerv (1:3) = Interpolate_cubic2DFd1 (ed_cellCubicNele (1:16,i,j,1), x,y)

           a1X = - c2div2nc * NeleDerv (2) * deltaInvX                  ! 1st half time step x-acceleration
           a1Y = - c2div2nc * NeleDerv (3) * deltaInvY                  ! 1st half time step y-acceleration

           r2X = r1X + (v1X + a1X * stepTimeFourth) * stepTimeHalf      ! 2nd half time step x-position
           r2Y = r1Y + (v1Y + a1Y * stepTimeFourth) * stepTimeHalf      ! 2nd half time step y-position

           badTimeStep =     (abs (r0X - r2X) > cellStepX) &            ! decide if time step is good/bad
                        .or. (abs (r0Y - r2Y) > cellStepY)              !

           if (badTimeStep) then

               stepTime     = stepTimeHalf
               stepTimeHalf = stepTimeFourth

               r0X = r1X                                                ! new ray endpoint x-position
               r0Y = r1Y                                                ! new ray endpoint y-position
               v0X = v1X                                                ! new ray endpoint x-velocity
               v0Y = v1Y                                                ! new ray endpoint y-velocity

               newCell = .false.                                        ! we stay inside the cell
           else
               exit
           end if

        end do
!
!
!     ...At this stage we have determined the parabolic ray path such that: 1) it is accurate
!        to the demanded path accuracy, 2) it lays entirely within the current cell and 3) we
!        determined the new end position, either still within the current cell (newCell = false)
!        or on one of its walls, ready to exit the cell (newCell = true).
!
!        Next we calculate the ray's power deposition as the ray travels along that parabolic
!        path. This is done by evaluating the time integral using Gaussian quadrature. Since
!        the parabolic path is known analytically, for each Gaussian quadarature time points
!        we can calculate the ray's position, from which we can calculate the Nele and Tele
!        values at that position. The Gaussian Quadrature with weight function equal to 1 is used.
!        The associated orthogonal polynomials are the Legendre Polynomials. If the remaining ray
!        power is considered to have reached a 'zero' value, mark the ray as nonexistent and exit
!        the indefinite block loop.
!
!        Currently a 2 point Gaussian quadrature is used (accurate to third order polynomial
!        expansion of the integrand in time).
!
!
        integral = 0.0

        do n = 1,2

           gaussTime     = GaussianRoot (n) * stepTimeHalf
           gaussTimeHalf = 0.5 * gaussTime        

           x = rayX + (velX + accX * gaussTimeHalf) * gaussTime
           y = rayY + (velY + accY * gaussTimeHalf) * gaussTime

           x = (x - xminCell) * deltaInvX               ! rescaled [0,1] ray x coordinate
           y = (y - yminCell) * deltaInvY               ! rescaled [0,1] ray y coordinate

           Nele = Interpolate_cubic2DF (ed_cellCubicNele (1:16,i,j,1), x,y)
           Tele = Interpolate_cubic2DF (ed_cellCubicTele (1:16,i,j,1), x,y)

           if (Nele <= 0.0) then
               call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: Nele <= 0 for a cell')
           end if

           if (Tele <= 0.0) then
               call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: Tele <= 0 for a cell')
           end if

           lnLambda = ed_CoulombFactor (cellZbar,          &
                                        ed_electronCharge, &
                                        ed_Boltzmann,      &
                                        Tele,              &
                                        Nele               )

           nu = ed_inverseBremsstrahlungRate (cellZbar,          &
                                              ed_electronCharge, &
                                              ed_electronMass,   &
                                              ed_Boltzmann,      &
                                              Tele,              &
                                              Nele,              &
                                              rayCritDens,       &
                                              lnLambda           )
           integral = integral + nu

        end do

        integral = integral * stepTimeHalf

        powerLossFactor = exp (-integral)
        cellPower       = rayPower * (1.0 - powerLossFactor)
        cellEnergy      = cellPower * timeStep

        if (ed_depoVarIsPerMass) then
            cellEnergyDepot (i,j) = cellEnergyDepot (i,j) + cellEnergy * cellMassInv
        else
            cellEnergyDepot (i,j) = cellEnergyDepot (i,j) + cellEnergy * cellVolumeInv
        end if

        rayPower = rayPower * powerLossFactor

        if (rayPower <= ed_rayZeroPower) then
            ed_rays (RAY_BLCK,ray) = real (NONEXISTENT)
            exit
        end if
!
!
!     ...Set the new ray position and the new ray velocity. If the ray is stationary
!        (no movement), mark the ray as nonexistent and exit the block loop. In case
!        a laser IO is performed on the rays, store the current ray IO data.
!
!
        rayX = r0X
        rayY = r0Y

        velX = v0X
        velY = v0Y

        velXeq0 = (velX == 0.0)
        velYeq0 = (velY == 0.0)
 
        stationaryRay = velXeq0 .and. velYeq0

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
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = 0.0
            end if
        end if
!
!
!     ...If, at the current stage, the ray enters a new cell, we have to determine: 1) which new
!        cell (i,j) it is and 2) the appropriate nudging values on the ray's position. Due to
!        possible reflective boundary conditions on the block faces, it can happen that the ray
!        stays in the original cell. After handling the logistics inside the following 'if'
!        statement, the new cell indices i,j are either the old ones or new ones.
!
!
        if (newCell) then

            dist2minX = abs (xminCell - rayX)
            dist2maxX = abs (xmaxCell - rayX)
            dist2minY = abs (yminCell - rayY)
            dist2maxY = abs (ymaxCell - rayY)

            minDistance = min (dist2minX, dist2maxX, &
                               dist2minY, dist2maxY)

            if (minDistance > cellWallThicknessHalf) then
                call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: ray to far away from cell face')
            end if

            cellFaceMinX = (dist2minX <= cellWallThicknessHalf)
            cellFaceMaxX = (dist2maxX <= cellWallThicknessHalf)
            cellFaceMinY = (dist2minY <= cellWallThicknessHalf)
            cellFaceMaxY = (dist2maxY <= cellWallThicknessHalf)

            velXgt0 = (velX  > 0.0)
            velXlt0 = (velX  < 0.0)
            velYgt0 = (velY  > 0.0)
            velYlt0 = (velY  < 0.0)

            crossX = .false.
            crossY = .false.

            nudgeX = 0.0
            nudgeY = 0.0

            ip = i
            jp = j

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

            blockFaceMinX = (rayX == xminBlock)
            blockFaceMaxX = (rayX == xmaxBlock)
            blockFaceMinY = (rayY == yminBlock)
            blockFaceMaxY = (rayY == ymaxBlock)

            reflectX =     (blockFaceMinX .and. blockReflectMinX .and. velXlt0) &
                      .or. (blockFaceMaxX .and. blockReflectMaxX .and. velXgt0)
            reflectY =     (blockFaceMinY .and. blockReflectMinY .and. velYlt0) &
                      .or. (blockFaceMaxY .and. blockReflectMaxY .and. velYgt0)

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

            if (crossX) then
                nudgeX = (i - ip) * cellWallThicknessHalf
            end if

            if (crossY) then
                nudgeY = (j - jp) * cellWallThicknessHalf
            end if

            rayX = rayX + nudgeX
            rayY = rayY + nudgeY

            newCell = crossX .or. crossY

        end if
!
!
!     ...We are now sure about the target cell. Check, if the target cell (i,j) is still within the block.
!        If it is, we check if this is a new cell, in which case we update the cell info. If the target
!        cell is not within the block, check if the ray coordinates are still within the defined domain.
!        If not, store its latest data and mark it as nonexistent. If the ray is still within the domain
!        boundaries, exit the current block loop.
!
!
        inBlock =      (i >= iminBlock) &
                 .and. (i <= imaxBlock) &
                 .and. (j >= jminBlock) &
                 .and. (j <= jmaxBlock)

        if (inBlock) then

            if (newCell) then

                xminCell = ed_cellEdges   (i  ,1)
                xmaxCell = ed_cellEdges   (i+1,1)
                yminCell = ed_cellEdges   (j  ,2)
                ymaxCell = ed_cellEdges   (j+1,2)

                cellZbar      = ed_cellZbar    (i,j,1)
                cellDensity   = ed_cellDensity (i,j,1)
                cellVolume    = ed_cellVolume  (i,j,1)
                cellVolumeInv = 1.0 / cellVolume
                cellMass      = cellDensity * cellVolume
                cellMassInv   = 1.0 / cellMass

            end if

        else

            ed_rays (RAY_POSX,ray) = rayX
            ed_rays (RAY_POSY,ray) = rayY
            ed_rays (RAY_VELX,ray) = velX
            ed_rays (RAY_VELY,ray) = velY
            ed_rays (RAY_POWR,ray) = rayPower

            inDomain =      (rayX > ed_xminDomain) &
                      .and. (rayX < ed_xmaxDomain) &
                      .and. (rayY > ed_yminDomain) &
                      .and. (rayY < ed_ymaxDomain)

            if (.not.inDomain) then
                 ed_rays (RAY_POSX,ray) = rayX - nudgeX             ! undo the nudging, since
                 ed_rays (RAY_POSY,ray) = rayY - nudgeY             ! it is not needed anymore
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
        print *, "[ed_traceBlockRays2DRec] Ray ", ray, &
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
end subroutine ed_traceBlockRays2DRec
