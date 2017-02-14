!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/PPRT/ed_traceBlockRays1DRec
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
                                     ed_xmaxDomain

  use Driver_interface,       ONLY : Driver_abortFlash

  use Interpolate_interface,  ONLY : Interpolate_cubic1DF,  &
                                     Interpolate_cubic1DFd1

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
  logical :: impossibleRay
  logical :: inDomain, inBlock
  logical :: newCell
  logical :: onBlockBoundaryCell
  logical :: rayOutOfBlock
  logical :: reflectX
  logical :: stationaryRay
  logical :: velXeq0, velXgt0, velXlt0
  logical :: writeRay

  integer :: i,n
  integer :: ip
  integer :: nRayWritePositions
  integer :: numDeadRays
  integer :: ray
  integer :: rayTag
  integer :: rayWriteIndex

  real    :: a1X
  real    :: accX
  real    :: c2div2nc
  real    :: cellDensity
  real    :: cellEnergy
  real    :: cellMass, cellMassInv
  real    :: cellPower
  real    :: cellStepX
  real    :: cellVolume, cellVolumeInv
  real    :: cellWallThicknessHalf
  real    :: cellZbar
  real    :: dist2minX, dist2maxX
  real    :: dX
  real    :: gaussTime, gaussTimeHalf
  real    :: integral
  real    :: lnLambda
  real    :: minDistance
  real    :: Nele, Tele
  real    :: nu
  real    :: nudgeX
  real    :: powerLossFactor
  real    :: r0X, r1X, r2X
  real    :: rayCritDens
  real    :: rayPower
  real    :: rayX
  real    :: stepTime, stepTimeHalf, stepTimeFourth
  real    :: v0X, v1X
  real    :: time2Face
  real    :: velX
  real    :: x
  real    :: xminCell, xmaxCell

  real    :: NeleDerv (1:2)

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
!     ...Get extra needed info about the initial cell (i).
!
!
     cellZbar      = ed_cellZbar    (i,1,1)
     cellDensity   = ed_cellDensity (i,1,1)
     cellVolume    = ed_cellVolume  (i,1,1)
     cellVolumeInv = 1.0 / cellVolume
     cellMass      = cellDensity * cellVolume
     cellMassInv   = 1.0 / cellMass
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

        NeleDerv (1:2) = Interpolate_cubic1DFd1 (ed_cellCubicNele (1:4,i,1,1), x)

        dX  = NeleDerv (2) * deltaInvX                                  ! d/dx -> d/dX

        accX = - c2div2nc * dX                                          ! acceleration in x-direction

        stepTime = ed_time2FacesParabolicPath1D (rayX, velX, accX, xminCell, xmaxCell, ed_infiniteTime)

        if (stepTime == ed_infiniteTime .or. stepTime == 0.0) then
            call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: infinite/zero crossing time for a cell')
        end if

        stepTimeHalf = 0.5  * stepTime

        r0X = rayX + (velX + accX * stepTimeHalf) * stepTime            ! full time step x-position
        v0X = velX + accX * stepTime                                    ! full time step x-velocity

        newCell = .true.

        do                                                              ! check ray path inside cell

           stepTimeFourth = 0.5 * stepTimeHalf

           r1X = rayX + (velX + accX * stepTimeFourth) * stepTimeHalf   ! 1st half time step x-position
           v1X = velX + accX * stepTimeHalf                             ! 1st half time step x-velocity

           x = (r1X - xminCell) * deltaInvX                             ! rescaled [0,1] ray x coordinate

           NeleDerv (1:2) = Interpolate_cubic1DFd1 (ed_cellCubicNele (1:4,i,1,1), x)

           a1X = - c2div2nc * NeleDerv (2) * deltaInvX                  ! 1st half time step x-acceleration
           r2X = r1X + (v1X + a1X * stepTimeFourth) * stepTimeHalf      ! 2nd half time step x-position

           badTimeStep = (abs (r0X - r2X) > cellStepX)                  ! decide if time step is good/bad

           if (badTimeStep) then

               stepTime     = stepTimeHalf
               stepTimeHalf = stepTimeFourth

               r0X = r1X                                                ! new ray endpoint x-position
               v0X = v1X                                                ! new ray endpoint x-velocity

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
           x = (x - xminCell) * deltaInvX                         ! rescaled [0,1] ray x coordinate

           Nele = Interpolate_cubic1DF (ed_cellCubicNele (1:4,i,1,1), x)
           Tele = Interpolate_cubic1DF (ed_cellCubicTele (1:4,i,1,1), x)

           if (Nele <= 0.0) then
               call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: Nele <= 0 for a cell')
           end if

           if (Tele <= 0.0) then
               call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: Tele <= 0 for a cell')
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
            cellEnergyDepot (i) = cellEnergyDepot (i) + cellEnergy * cellMassInv
        else
            cellEnergyDepot (i) = cellEnergyDepot (i) + cellEnergy * cellVolumeInv
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
        velX = v0X

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
