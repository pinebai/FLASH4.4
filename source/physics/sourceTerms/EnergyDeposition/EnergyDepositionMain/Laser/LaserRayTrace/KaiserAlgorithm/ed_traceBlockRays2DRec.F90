!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/KaiserAlgorithm/ed_traceBlockRays2DRec
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
!!
!!***

subroutine ed_traceBlockRays2DRec (timeStep,                          &
                                   rayFirst,  rayLast,                &
                                   iminBlock, imaxBlock,              &
                                   jminBlock, jmaxBlock,              &
                                   xminBlock, xmaxBlock,              &
                                   yminBlock, ymaxBlock,              &
                                   deltaX, deltaY,                    &
                                   deltaInvX, deltaInvY,              &
                                   blockReflectMinX,                  &
                                   blockReflectMaxX,                  &
                                   blockReflectMinY,                  &
                                   blockReflectMaxY,                  &
                                                      cellEnergyDepot ) 

  use EnergyDeposition_data,  ONLY : ed_Boltzmann,                   &
                                     ed_cellCenters,                 &
                                     ed_cellDensity,                 &
                                     ed_cellEdges,                   &
                                     ed_cellGradNele,                &
                                     ed_cellGradTele,                &
                                     ed_cellNele,                    &
                                     ed_cellTele,                    &
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

  use ed_interface,           ONLY : ed_CoulombFactor,             &
                                     ed_inverseBremsstrahlungRate, &
                                     ed_time2FacesParabolicPath1D

  use ed_commInterface,       ONLY : ed_commProgressTransport, &
                                     ed_commHandleOffBlkRay,   &
                                     ed_commIncrementDeadRays

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
  real,    intent (in)    :: deltaX, deltaY
  real,    intent (in)    :: deltaInvX, deltaInvY
  logical, intent (in)    :: blockReflectMinX
  logical, intent (in)    :: blockReflectMaxX
  logical, intent (in)    :: blockReflectMinY
  logical, intent (in)    :: blockReflectMaxY
  real,    intent (inout) :: cellEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock)

  logical :: blockFaceMinX, blockFaceMaxX
  logical :: blockFaceMinY, blockFaceMaxY
  logical :: cellFaceMinX, cellFaceMaxX
  logical :: cellFaceMinY, cellFaceMaxY
  logical :: crossX, crossY
  logical :: inDomain, inBlock
  logical :: onBlockBoundaryCell
  logical :: rayCrossesBoundaryX, rayCrossesBoundaryY
  logical :: rayOutOfBlock
  logical :: reflectX, reflectY
  logical :: stationaryRay
  logical :: velXeq0, velXgt0, velXlt0
  logical :: velYeq0, velYgt0, velYlt0
  logical :: writeRay

  integer :: i,j
  integer :: ip,jp
  integer :: nRayWritePositions
  integer :: numDeadRays
  integer :: ray
  integer :: rayTag
  integer :: rayWriteIndex

  real    :: a,b,c,d
  real    :: accX,accY
  real    :: c2div1nc, c2div2nc, c2div4nc
  real    :: cellDensity
  real    :: cellEnergy
  real    :: cellMass
  real    :: cellPower
  real    :: cellVolume
  real    :: cellWallThicknessHalf
  real    :: cellZbar
  real    :: centerCellX, centerCellY
  real    :: centerNele, centerNeleOld, centerTele
  real    :: crossTime, crossTimeHalf
  real    :: dist2minX, dist2minY
  real    :: dist2maxX, dist2maxY
  real    :: distX, distY
  real    :: distXOld, distYOld
  real    :: gradNeleX, gradNeleY
  real    :: gradNeleXOld, gradNeleYOld
  real    :: gradTeleX, gradTeleY
  real    :: integral
  real    :: lnLambda
  real    :: minDistance
  real    :: Nele, NeleOld, Tele
  real    :: nu
  real    :: nudgeX, nudgeY
  real    :: powerLossFactor
  real    :: R,S,U,W
  real    :: rayCritDens, rayCritDensInv
  real    :: rayPower
  real    :: rayX, rayY
  real    :: timeX, timeY
  real    :: velX, velY
  real    :: vnewSqr
  real    :: xminCell, yminCell
  real    :: xmaxCell, ymaxCell

  real, parameter :: GaussianRoot1 = 1.577350269189626         ! is 1 + 1/sqrt(3) for integration limits [0,1]
  real, parameter :: GaussianRoot2 = 4.226497308103742e-01     ! is 1 - 1/sqrt(3) for integration limits [0,1]
!
!
!     ...Define some variables.
!
!
  numDeadRays = 0

  cellWallThicknessHalf = 0.5 * ed_cellWallThickness
!
!
!     ...Outer (threaded) loop over all rays associated with the current block.
!
!
!$omp do schedule (dynamic)
  do ray = rayFirst , rayLast
     
     call ed_commProgressTransport ()
   
     rayTag         = int (ed_rays (RAY_TAGS,ray))
     rayX           =      ed_rays (RAY_POSX,ray)
     rayY           =      ed_rays (RAY_POSY,ray)
     velX           =      ed_rays (RAY_VELX,ray)
     velY           =      ed_rays (RAY_VELY,ray)
     rayPower       =      ed_rays (RAY_POWR,ray)
     rayCritDens    =      ed_rays (RAY_DENC,ray)
     rayCritDensInv = 1.0 / rayCritDens

     c2div1nc = ed_speedOfLightSquared * rayCritDensInv
     c2div4nc = 0.25 * c2div1nc
     c2div2nc = c2div4nc + c2div4nc
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
!     ...Calculate the electron density and the electron temperature at the current ray
!        position in the (i,j) cell. These are calculated via interpolation from the
!        corresponding gradient values. If either of these values become negative, the
!        calculation must be stopped.
!
!
     centerCellX = ed_cellCenters  (i,1)
     centerCellY = ed_cellCenters  (j,2)
     centerNele  = ed_cellNele     (  i,j,1)
     centerTele  = ed_cellTele     (  i,j,1)
     gradNeleX   = ed_cellGradNele (1,i,j,1)
     gradNeleY   = ed_cellGradNele (2,i,j,1)
     gradTeleX   = ed_cellGradTele (1,i,j,1)
     gradTeleY   = ed_cellGradTele (2,i,j,1)

     distX = rayX - centerCellX
     distY = rayY - centerCellY
     Nele  = centerNele  +  gradNeleX * distX  +  gradNeleY * distY
     Tele  = centerTele  +  gradTeleX * distX  +  gradTeleY * distY

     if (Nele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: Nele <= 0 for a cell (initial)')
     end if

     if (Tele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: Tele <= 0 for a cell (initial)')
     end if
!
!
!     ...We are ready to follow the ray's path through all the cells of the current
!        block. The current cell indices (i,j) and the previous cell indices (ip,jp)
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

        xminCell    = ed_cellEdges   (i  ,1)
        xmaxCell    = ed_cellEdges   (i+1,1)
        yminCell    = ed_cellEdges   (j  ,2)
        ymaxCell    = ed_cellEdges   (j+1,2)

        cellZbar    = ed_cellZbar    (i,j,1)
        cellVolume  = ed_cellVolume  (i,j,1)
        cellDensity = ed_cellDensity (i,j,1)
        cellMass    = cellDensity * cellVolume
!
!
!     ...The ray is being refracted through the cell (i,j). We have to follow its path
!        through the cell and determine its exit cell. At that point the previous cell
!        becomes simply the current cell, that is cell (ip,jp) = cell (i,j) and the loop
!        is closed. Solve next the quadratic (parabolic) time equations to find out when
!        the ray will hit the next cell face.
!
!
        accX = - c2div2nc * gradNeleX               ! acceleration in x-direction
        accY = - c2div2nc * gradNeleY               ! acceleration in y-direction

        timeX = ed_time2FacesParabolicPath1D (rayX, velX, accX, xminCell, xmaxCell, ed_infiniteTime)
        timeY = ed_time2FacesParabolicPath1D (rayY, velY, accY, yminCell, ymaxCell, ed_infiniteTime)

        crossTime = min (timeX, timeY)

        if (crossTime == ed_infiniteTime .or. crossTime == 0.0) then
            call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: infinite/zero cell crossing time')
        end if

        crossTimeHalf = 0.5 * crossTime
!
!
!     ...We found a reasonable crossing time to a particular cell face plane. Calculate the
!        power deposition as the ray traverses the cell. This is done by evaluating an integral
!        using the initial (before crossing) velocities. The method of integral evaluation is
!        Gaussian Quadrature with weight function equal to 1. The associated orthogonal
!        polynomials are the Legendre Polynomials. If the remaining ray power is considered to
!        have reached a 'zero' value, mark the ray as nonexistent and exit the indefinite loop.
!
!
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

        U =   (     velX * gradNeleX +      velY * gradNeleY) / Nele
        W =   (     velX * gradTeleX +      velY * gradTeleY) / Tele
        R = - (gradNeleX * gradNeleX + gradNeleY * gradNeleY) * c2div4nc / Nele
        S = - (gradNeleX * gradTeleX + gradNeleY * gradTeleY) * c2div4nc / Tele

        a = GaussianRoot1 * crossTimeHalf
        b = a * a
        c = 1.0 + U * a + R * b
        d = sqrt (1.0 + W * a + S * b)

        integral = c * c / (d * d * d)

        a = GaussianRoot2 * crossTimeHalf
        b = a * a
        c = 1.0 + U * a + R * b
        d = sqrt (1.0 + W * a + S * b)

        integral = integral + c * c / (d * d * d)
        integral = integral * nu * crossTimeHalf

        powerLossFactor = exp (-integral)
        cellPower       = rayPower * (1.0 - powerLossFactor)
        cellEnergy      = cellPower * timeStep

        if (ed_depoVarIsPerMass) then
            cellEnergyDepot (i,j) = cellEnergyDepot (i,j) + cellEnergy / cellMass
        else
            cellEnergyDepot (i,j) = cellEnergyDepot (i,j) + cellEnergy / cellVolume
        end if

        rayPower = rayPower * powerLossFactor

        if (rayPower <= ed_rayZeroPower) then
            ed_rays (RAY_BLCK,ray) = real (NONEXISTENT)
            numDeadRays = numDeadRays + 1
            exit
        end if
!
!
!     ...Calculate the new ray position and velocities and check, if the ray crosses
!        a cell face. In case a laser IO is performed on the rays, store the current
!        ray IO data.
!
!
        rayX = rayX + (velX + accX * crossTimeHalf) * crossTime
        rayY = rayY + (velY + accY * crossTimeHalf) * crossTime

        dist2minX = abs (xminCell - rayX)
        dist2maxX = abs (xmaxCell - rayX)
        dist2minY = abs (yminCell - rayY)
        dist2maxY = abs (ymaxCell - rayY)

        minDistance = min (dist2minX, dist2maxX, &
                           dist2minY, dist2maxY)

        if (minDistance > cellWallThicknessHalf) then
            call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: ray too far away from cell face')
        end if

        cellFaceMinX = (dist2minX <= cellWallThicknessHalf)
        cellFaceMaxX = (dist2maxX <= cellWallThicknessHalf)
        cellFaceMinY = (dist2minY <= cellWallThicknessHalf)
        cellFaceMaxY = (dist2maxY <= cellWallThicknessHalf)

        velX = velX + accX * crossTime
        velY = velY + accY * crossTime

        velXgt0 = (velX  > 0.0)
        velXeq0 = (velX == 0.0)
        velXlt0 = (velX  < 0.0)
        velYgt0 = (velY  > 0.0)
        velYeq0 = (velY == 0.0)
        velYlt0 = (velY  < 0.0)

        stationaryRay = velXeq0 .and. velYeq0

        if (stationaryRay) then
            write (*,*) ' stationary ray detected! Removing it from list ... '
            ed_rays (RAY_BLCK,ray) = real (NONEXISTENT)
            numDeadRays = numDeadRays + 1
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
!     ...Calculate the cell indices of the new cell where the ray moves into. Again, the ray
!        is (at least) on one of the four faces and moves into the next cell from the current
!        cell, which becomes cell (ip,jp). The ray can be on:
!
!                       1 cell face  -> within the cell face
!                       2 cell faces -> in the corner between the cell faces
!
!        The goal here is now to determine into which cell (i,j) the ray will go. Cell (i,j)
!        will in most of the cases be different from cell (ip,jp), but occasionally both cells
!        can remain identical. If the ray's position is considered to be close to a face and it
!        is certain that the ray will cross that face, the ray is forced (temporarily) to be
!        exactly on that face. This is done in order to calculate the correct refraction
!        properties of the ray, if needed. The nudging value(s) and their directions to keep
!        the ray in the current cell are also set here.
!
!        The variables 'cross(X,Y)' indicate, if the ray crosses the corresponding cell boundaries
!        and Snell's law has to be applied to the X,Y components of the ray's velocity in case of
!        refraction. Note, that if at least one of the 'cross(X,Y)' is true, then the ray crosses
!        a cell boundary from cell (ip,jp) to cell (i,j).
!
!
        ip = i
        jp = j

        crossX = .false.
        crossY = .false.

        nudgeX = 0.0
        nudgeY = 0.0

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
!
!
!     ...At this stage, the ray situation is as follows:
!
!
!                                  ---------- ----------
!                                 |          |          |
!                                 |          |          |          P = previous ray position
!                                 |     *    N     +    |          N = current ray position
!                                 P          |          |        *,+ = cell centers
!                                 |   ip,jp  |    i,j   |
!                                  ---------- ----------
!
!        The new target cell (i,j) is not final yet. It depends on the outcome of Snell's law and
!        the block face boundary conditions. Check first, if a block face reflection needs to be honored.
!        After that, check if Snell's law needs to be applied to find the new velocity components.
!        For this we need to know the number of electrons in the previous old cell and the new cell
!        right at the current boundary. Again this is calculated via interpolation from the gradient
!        values. If the ray reflects, the new target cell is updated and the corresponding velocity
!        component is inverted. If the ray still crosses a cell boundary, determine the proper
!        nudging for the ray and calculate the new nudged ray positions.
!
!
        blockFaceMinX = (rayX == xminBlock)
        blockFaceMaxX = (rayX == xmaxBlock)
        blockFaceMinY = (rayY == yminBlock)
        blockFaceMaxY = (rayY == ymaxBlock)

        reflectX =      (blockFaceMinX .and. blockReflectMinX .and. velXlt0) &
                   .or. (blockFaceMaxX .and. blockReflectMaxX .and. velXgt0)
        reflectY =      (blockFaceMinY .and. blockReflectMinY .and. velYlt0) &
                   .or. (blockFaceMaxY .and. blockReflectMaxY .and. velYgt0)

        rayCrossesBoundaryX = (.not.reflectX .and. crossX)
        rayCrossesBoundaryY = (.not.reflectY .and. crossY)

        if (rayCrossesBoundaryX .or. rayCrossesBoundaryY) then

            centerNeleOld = centerNele                          ! for cell (ip,jp) at center *
            gradNeleXOld  = gradNeleX                           ! for cell (ip,jp) at center *
            gradNeleYOld  = gradNeleY                           ! for cell (ip,jp) at center *
            distXOld      = rayX - centerCellX                  ! for cell (ip,jp) at position N
            distYOld      = rayY - centerCellY                  ! for cell (ip,jp) at position N
            centerNele    = ed_cellNele       (  i,j,1)         ! for cell (i ,j ) at center +
            gradNeleX     = ed_cellGradNele   (1,i,j,1)         ! for cell (i ,j ) at center +
            gradNeleY     = ed_cellGradNele   (2,i,j,1)         ! for cell (i ,j ) at center +
            distX         = rayX - ed_cellCenters (i,1)         ! for cell (i ,j ) at position N
            distY         = rayY - ed_cellCenters (j,2)         ! for cell (i ,j ) at position N

            NeleOld = centerNeleOld + gradNeleXOld * distXOld + gradNeleYOld * distYOld
            Nele    = centerNele    + gradNeleX    * distX    + gradNeleY    * distY

            if (NeleOld <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: NeleOld < 0 for a cell (Snell)')
            end if

            if (Nele <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: Nele <= 0 for a cell (Snell)')
            end if

            if (rayCrossesBoundaryX) then
                vNewSqr  = velX * velX + (NeleOld - Nele) * c2div1nc
                reflectX = vNewSqr < 0.0
                if (.not.reflectX) then
                     velX = sign (1.,velX) * sqrt (vNewSqr)
                end if
            end if

            if (rayCrossesBoundaryY) then
                vNewSqr  = velY * velY + (NeleOld - Nele) * c2div1nc
                reflectY = vNewSqr < 0.0
                if (.not.reflectY) then
                     velY = sign (1.,velY) * sqrt (vNewSqr)
                end if
            end if

        end if

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
!
!
!     ...We are now sure about the target cell. Check, if the target cell (i,j) is still within the block.
!        If it is, we have to calculate the new electron density and the new electron temperature where
!        the ray is located in the target cell. If the target cell is not within the block, check if the
!        ray coordinates are still within the defined domain. If not, store its latest data and mark it
!        as nonexistent. If the ray is still within the domain boundaries, exit the current block loop.
!
!
        inBlock =      (i >= iminBlock) &
                 .and. (i <= imaxBlock) &
                 .and. (j >= jminBlock) &
                 .and. (j <= jmaxBlock)

        if (inBlock) then

            centerCellX = ed_cellCenters  (i,1)
            centerCellY = ed_cellCenters  (j,2)
            centerNele  = ed_cellNele     (  i,j,1)
            centerTele  = ed_cellTele     (  i,j,1)
            gradNeleX   = ed_cellGradNele (1,i,j,1)
            gradNeleY   = ed_cellGradNele (2,i,j,1)
            gradTeleX   = ed_cellGradTele (1,i,j,1)
            gradTeleY   = ed_cellGradTele (2,i,j,1)

            distX   = rayX - centerCellX
            distY   = rayY - centerCellY
            Nele    = centerNele  +  gradNeleX * distX  +  gradNeleY * distY
            Tele    = centerTele  +  gradTeleX * distX  +  gradTeleY * distY

            if (Nele <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: Nele <= 0 for a cell (target)')
            end if

            if (Tele <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays2DRec] ERROR: Tele <= 0 for a cell (target)')
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
                 numDeadRays = numDeadRays + 1
            else
                 call ed_commHandleOffBlkRay(ray)
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
     call ed_commIncrementDeadRays(numDeadRays)
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_traceBlockRays2DRec
