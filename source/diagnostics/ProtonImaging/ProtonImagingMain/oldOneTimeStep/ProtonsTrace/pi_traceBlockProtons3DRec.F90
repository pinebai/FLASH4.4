!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonsTrace/pi_traceBlockProtons3DRec
!!
!! NAME
!!
!!  pi_traceBlockProtons3DRec
!!
!! SYNOPSIS
!!
!!  call pi_traceBlockProtons3DRec (integer (in)    :: protonFirst
!!                                  integer (in)    :: protonLast,
!!                                  integer (in)    :: iminBlock,
!!                                  integer (in)    :: imaxBlock,
!!                                  integer (in)    :: jminBlock,
!!                                  integer (in)    :: jmaxBlock,
!!                                  integer (in)    :: kminBlock,
!!                                  integer (in)    :: kmaxBlock,
!!                                  real    (in)    :: xminBlock,
!!                                  real    (in)    :: xmaxBlock,
!!                                  real    (in)    :: yminBlock,
!!                                  real    (in)    :: ymaxBlock,
!!                                  real    (in)    :: zminBlock,
!!                                  real    (in)    :: zmaxBlock,
!!                                  real    (in)    :: deltaX,
!!                                  real    (in)    :: deltaY,
!!                                  real    (in)    :: deltaZ,
!!                                  real    (in)    :: deltaInvX,
!!                                  real    (in)    :: deltaInvY,
!!                                  real    (in)    :: deltaInvZ,
!!                                  logical (in)    :: blockReflectMinX,
!!                                  logical (in)    :: blockReflectMaxX,
!!                                  logical (in)    :: blockReflectMinY,
!!                                  logical (in)    :: blockReflectMaxY,
!!                                  logical (in)    :: blockReflectMinZ,
!!                                  logical (in)    :: blockReflectMaxZ)
!!
!! DESCRIPTION
!!
!!  Traces the movement of the current collection of active protons through one block
!!  for those geometries consisting formally of 3D rectangular grids (cartesian).
!!  On exit, each proton has either:
!!
!!            i)  reached a different (yet unknown) block
!!           ii)  has reached the domain boundary and exited.
!!
!! ARGUMENTS
!!
!!  protonFirst      : first proton index to be considered
!!  protonLast       : last proton index to be considered
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
!!
!! NOTES
!!        
!!  The code allows for threading to be used on the outer proton trace loop.
!!  The paths of the protons are computed using the average Lorentz force for each cell. 
!!
!!***

subroutine pi_traceBlockProtons3DRec (protonFirst,  protonLast,          &
                                      iminBlock, imaxBlock,              &
                                      jminBlock, jmaxBlock,              &
                                      kminBlock, kmaxBlock,              &
                                      xminBlock, xmaxBlock,              &
                                      yminBlock, ymaxBlock,              &
                                      zminBlock, zmaxBlock,              &
                                      deltaX, deltaY, deltaZ,            &
                                      deltaInvX, deltaInvY, deltaInvZ,   &
                                      blockReflectMinX,                  &
                                      blockReflectMaxX,                  &
                                      blockReflectMinY,                  &
                                      blockReflectMaxY,                  &
                                      blockReflectMinZ,                  &
                                      blockReflectMaxZ                   ) 

  use ProtonImaging_data,      ONLY : pi_beams,                       &
                                      pi_cellBfield,                  &
                                      pi_cellBoundary,                &
                                      pi_cellCurlBfield,              &
                                      pi_cellEdgesX,                  &
                                      pi_cellEdgesY,                  &
                                      pi_cellEdgesZ,                  &
                                      pi_cellEfield,                  &
                                      pi_cellStepTolerance,           &
                                      pi_cellWallThickness,           &
                                      pi_infiniteSpeed,               &
                                      pi_infiniteTime,                &
                                      pi_IOaddProtonsDomain2Screen,   &
                                      pi_IOmaxPointsPerBlock,         &
                                      pi_IOmaxProtonCount,            &
                                      pi_IOplotProtons,               &
                                      pi_IOprotonCount,               &
                                      pi_IOprotonPointCount,          &
                                      pi_IOprotonPoints,              &
                                      pi_IOprotonTags,                &
                                      pi_IOprotonWriteModulo,         &
                                      pi_largestPositiveReal,         &
                                      pi_opaqueBoundaries,            &
                                      pi_protonChargePerMass,         &
                                      pi_protons,                     &
                                      pi_protonsMovedIntoDomain,      &
                                      pi_RungeKuttaMethod,            &
                                      pi_screenProtonDiagnostics,     &
                                      pi_speedOfLightInv,             &
                                      pi_unitRoundoff,                &
                                      pi_useParabolicApproximation,   &
                                      pi_xminDomain,                  &
                                      pi_xmaxDomain,                  &
                                      pi_yminDomain,                  &
                                      pi_ymaxDomain,                  &
                                      pi_zminDomain,                  &
                                      pi_zmaxDomain

  use RungeKutta_interface,    ONLY : RungeKutta_stepConfined

  use Driver_interface,        ONLY : Driver_abortFlash

  use pi_traceODEfunctionData, ONLY : Bx, By, Bz,                   &
                                      BxInvC, ByInvC, BzInvC,       &
                                      CurlBx, CurlBy, CurlBz,       &
                                      Ex, Ey, Ez,                   &
                                      Qm,                           &
                                      xmaxCell, ymaxCell, zmaxCell, &
                                      xminCell, yminCell, zminCell

  use pi_interface,            ONLY : pi_maxConfinement3DRec,       &
                                      pi_minConfinement3DRec,       &
                                      pi_parabolicPathLength3D,     &
                                      pi_recordProtonOnScreen,      &
                                      pi_traceODEfunction3DRec,     &
                                      pi_time2FacesParabolicPath1D

  implicit none

#include "constants.h"
#include "Flash.h"
#include "ProtonImaging.h"

  integer, intent (in) :: protonFirst, protonLast   
  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: jminBlock, jmaxBlock
  integer, intent (in) :: kminBlock, kmaxBlock
  real,    intent (in) :: xminBlock, xmaxBlock
  real,    intent (in) :: yminBlock, ymaxBlock
  real,    intent (in) :: zminBlock, zmaxBlock
  real,    intent (in) :: deltaX, deltaY, deltaZ
  real,    intent (in) :: deltaInvX, deltaInvY, deltaInvZ
  logical, intent (in) :: blockReflectMinX
  logical, intent (in) :: blockReflectMaxX
  logical, intent (in) :: blockReflectMinY
  logical, intent (in) :: blockReflectMaxY
  logical, intent (in) :: blockReflectMinZ
  logical, intent (in) :: blockReflectMaxZ

  logical :: blockFaceMinX, blockFaceMaxX
  logical :: blockFaceMinY, blockFaceMaxY
  logical :: blockFaceMinZ, blockFaceMaxZ
  logical :: cellFaceMinX,  cellFaceMaxX
  logical :: cellFaceMinY,  cellFaceMaxY
  logical :: cellFaceMinZ,  cellFaceMaxZ
  logical :: crossX, crossY, crossZ
  logical :: inDomain, inBlock
  logical :: IOproton
  logical :: newCell
  logical :: onBlockBoundaryCell
  logical :: onScreen
  logical :: outOfCell
  logical :: parabolicPath
  logical :: protonOutOfBlock
  logical :: reflectX, reflectY, reflectZ
  logical :: velXgt0, velXlt0
  logical :: velYgt0, velYlt0
  logical :: velZgt0, velZlt0

  integer :: beam
  integer :: detector
  integer :: i,j,k
  integer :: IOprotonCount, IOpointCount
  integer :: ip,jp,kp
  integer :: proton
  integer :: protonTag
  integer :: protonBlk

  real    :: accX, accY, accZ
  real    :: addX, addY, addZ
  real    :: cellWallThicknessHalf
  real    :: cellStepError, cellStepErrorX, cellStepErrorY, cellStepErrorZ
  real    :: crossTime, crossTimeHalf
  real    :: dgJv, dgKx, dgKy, dgKz
  real    :: dist2minX, dist2minY, dist2minZ
  real    :: dist2maxX, dist2maxY, dist2maxZ
  real    :: jrkX, jrkY, jrkZ
  real    :: jerkError
  real    :: minDistance
  real    :: nudgeX, nudgeY, nudgeZ
  real    :: parabolicPathLength
  real    :: protonErrorFrac
  real    :: posX, posY, posZ
  real    :: scrX, scrY, scrZ
  real    :: stepTimeTry, stepTimeUsed, stepTimeNext
  real    :: time, timeX, timeY, timeZ
  real    :: velX, velY, velZ

  real    :: protonError     (1:10)   !
  real    :: protonErrorBase (1:10)   ! full-sized arrays including
  real    :: protonIn        (1:10)   ! possible diagnostic variables
  real    :: protonOut       (1:10)   !

  real, parameter :: sixth = 1.0/6.0
!
!
!     ...Define some variables.
!
!
  Qm = pi_protonChargePerMass

  cellStepErrorX = deltaX * pi_cellStepTolerance
  cellStepErrorY = deltaY * pi_cellStepTolerance
  cellStepErrorZ = deltaZ * pi_cellStepTolerance
  cellStepError  = min (cellStepErrorX, cellStepErrorY, cellStepErrorZ)
  cellWallThicknessHalf = 0.5 * pi_cellWallThickness

  protonErrorFrac = 1.0                          ! this means the error is controlled by the base only

  protonErrorBase (1)  = cellStepErrorX          ! error bar on posX
  protonErrorBase (2)  = cellStepErrorY          ! error bar on posY
  protonErrorBase (3)  = cellStepErrorZ          ! error bar on posZ
  protonErrorBase (4)  = pi_infiniteSpeed        ! this in effect puts no error bars on velX yet
  protonErrorBase (5)  = pi_infiniteSpeed        ! this in effect puts no error bars on velY yet
  protonErrorBase (6)  = pi_infiniteSpeed        ! this in effect puts no error bars on velZ yet
  protonErrorBase (7)  = pi_largestPositiveReal  ! this in effect puts no error bars on dgJv yet
  protonErrorBase (8)  = pi_largestPositiveReal  ! this in effect puts no error bars on dgKx yet
  protonErrorBase (9)  = pi_largestPositiveReal  ! this in effect puts no error bars on dgKy yet
  protonErrorBase (10) = pi_largestPositiveReal  ! this in effect puts no error bars on dgKz yet
!
!
!     ...Outer (threaded) loop over all protons associated with the current block.
!
!
!$omp do schedule (dynamic)
  do proton = protonFirst , protonLast

     posX = pi_protons (PROTON_POSX,proton)
     posY = pi_protons (PROTON_POSY,proton)
     posZ = pi_protons (PROTON_POSZ,proton)
     velX = pi_protons (PROTON_VELX,proton)
     velY = pi_protons (PROTON_VELY,proton)
     velZ = pi_protons (PROTON_VELZ,proton)
     time = pi_protons (PROTON_TIME,proton)
     dgJv = pi_protons (PROTON_DGJV,proton)   ! get these, even if not needed
     dgKx = pi_protons (PROTON_DGKX,proton)   ! get these, even if not needed
     dgKy = pi_protons (PROTON_DGKY,proton)   ! get these, even if not needed
     dgKz = pi_protons (PROTON_DGKZ,proton)   ! get these, even if not needed

     protonTag  = int (pi_protons (PROTON_TAGS,proton))
     protonBlk  = int (pi_protons (PROTON_BLCK,proton))
!
!
!     ...Find the indices (i,j,k) of the initial cell through which the proton will
!        enter the block. We know for sure that the proton enters the block, because
!        otherwise it would not be on the current block list. Check, on which of
!        the six possible faces the proton currently is.
!
!
     protonOutOfBlock =     (posX < xminBlock) &
                       .or. (posX > xmaxBlock) &
                       .or. (posY < yminBlock) &
                       .or. (posY > ymaxBlock) &
                       .or. (posZ < zminBlock) &
                       .or. (posZ > zmaxBlock)

     if (protonOutOfBlock) then
         call Driver_abortFlash ('[pi_traceBlockProtons3DRec] ERROR: proton found out of block')
     end if

     dist2minX = posX - xminBlock
     dist2maxX = xmaxBlock - posX
     dist2minY = posY - yminBlock
     dist2maxY = ymaxBlock - posY
     dist2minZ = posZ - zminBlock
     dist2maxZ = zmaxBlock - posZ

     minDistance = min (dist2minX, dist2maxX, &
                        dist2minY, dist2maxY, &
                        dist2minZ, dist2maxZ)

     if (minDistance > pi_cellWallThickness) then
         call Driver_abortFlash ('[pi_traceBlockProtons3DRec] ERROR: proton too far inside the block')
     end if

     i = iminBlock + int ( (posX - xminBlock) * deltaInvX )
     j = jminBlock + int ( (posY - yminBlock) * deltaInvY )
     k = kminBlock + int ( (posZ - zminBlock) * deltaInvZ )

     onBlockBoundaryCell = (     (i == iminBlock) &
                            .or. (i == imaxBlock) &
                            .or. (j == jminBlock) &
                            .or. (j == jmaxBlock) &
                            .or. (k == kminBlock) &
                            .or. (k == kmaxBlock) )

     if (.not.onBlockBoundaryCell) then
         call Driver_abortFlash ('[pi_traceBlockProtons3DRec] ERROR: proton not in a block boundary cell')
     end if
!
!
!     ...If the opaque boundary condition is on, check if the current cell belongs to
!        a boundary and stop the proton if necessary.
!
!
     if (pi_opaqueBoundaries .and. pi_cellBoundary (i,j,k) > 0.0) then
         pi_protons (PROTON_BLCK,proton) = real (NONEXISTENT)
         cycle
     end if
!
!
!     ...The proton enters the cell for sure. Proceed with the tracing.
!
!
     xminCell = pi_cellEdgesX (i  )
     xmaxCell = pi_cellEdgesX (i+1)
     yminCell = pi_cellEdgesY (j  )
     ymaxCell = pi_cellEdgesY (j+1)
     zminCell = pi_cellEdgesZ (k  )
     zmaxCell = pi_cellEdgesZ (k+1)

     dist2minX = abs (xminCell - posX)
     dist2maxX = abs (xmaxCell - posX)
     dist2minY = abs (yminCell - posY)
     dist2maxY = abs (ymaxCell - posY)
     dist2minZ = abs (zminCell - posZ)
     dist2maxZ = abs (zmaxCell - posZ)

     cellFaceMinX = (dist2minX <= cellWallThicknessHalf)
     cellFaceMaxX = (dist2maxX <= cellWallThicknessHalf)
     cellFaceMinY = (dist2minY <= cellWallThicknessHalf)
     cellFaceMaxY = (dist2maxY <= cellWallThicknessHalf)
     cellFaceMinZ = (dist2minZ <= cellWallThicknessHalf)
     cellFaceMaxZ = (dist2maxZ <= cellWallThicknessHalf)
!
!
!     ...Make sure the proton is also properly nudged into the corresponding cell.
!
!
     if (cellFaceMinX) posX = xminCell + cellWallThicknessHalf
     if (cellFaceMaxX) posX = xmaxCell - cellWallThicknessHalf
     if (cellFaceMinY) posY = yminCell + cellWallThicknessHalf
     if (cellFaceMaxY) posY = ymaxCell - cellWallThicknessHalf
     if (cellFaceMinZ) posZ = zminCell + cellWallThicknessHalf
     if (cellFaceMaxZ) posZ = zmaxCell - cellWallThicknessHalf
!
!
!     ...Get the electric and magnetic field components at the current proton
!        position in the (i,j,k) cell. If diagnostic values are wanted, get extra
!        magnetic data.
!
!
     Ex = pi_cellEfield (1,i,j,k)                  !
     Ey = pi_cellEfield (2,i,j,k)                  ! units: Gauss
     Ez = pi_cellEfield (3,i,j,k)                  !

     if (pi_screenProtonDiagnostics) then

         Bx     = pi_cellBfield (1,i,j,k)          !
         By     = pi_cellBfield (2,i,j,k)          ! units: Gauss
         Bz     = pi_cellBfield (3,i,j,k)          !

         BxInvC = Bx * pi_speedOfLightInv          !
         ByInvC = By * pi_speedOfLightInv          ! units: Gauss / c
         BzInvC = Bz * pi_speedOfLightInv          !

         CurlBx = pi_cellCurlBfield (1,i,j,k)      !
         CurlBy = pi_cellCurlBfield (2,i,j,k)      ! units: Gauss / cm
         CurlBz = pi_cellCurlBfield (3,i,j,k)      !

     else

         BxInvC = pi_cellBfield (1,i,j,k)          !
         ByInvC = pi_cellBfield (2,i,j,k)          ! units: Gauss / c
         BzInvC = pi_cellBfield (3,i,j,k)          !

     end if
!
!
!     ...Decide, if this proton is one of the IO protons. If the case, start the
!        IO proton writeout procedure. If threading is done, the IO proton count
!        must be protected from incrementation by another thread and is saved in
!        a local thread variable.
!
!
     IOproton = pi_IOplotProtons

     if (IOproton) then
         IOproton = mod (protonTag, pi_IOprotonWriteModulo) == 0
         if (IOproton) then
             !$omp critical (setIOprotonCount)
                   IOproton = pi_IOprotonCount < pi_IOmaxProtonCount
                   if (IOproton) then
                       pi_IOprotonCount = pi_IOprotonCount + 1
                       IOprotonCount = pi_IOprotonCount
                   end if
             !$omp end critical (setIOprotonCount)
         end if
     end if

     if (IOproton) then
         IOpointCount = 1
         pi_IOprotonTags       (              IOprotonCount       ) = protonTag
         pi_IOprotonPointCount (              IOprotonCount       ) = IOpointCount
         pi_IOprotonPoints     (IOpointCount, IOprotonCount, IAXIS) = posX
         pi_IOprotonPoints     (IOpointCount, IOprotonCount, JAXIS) = posY
         pi_IOprotonPoints     (IOpointCount, IOprotonCount, KAXIS) = posZ
     end if
!
!
!     ...We are ready to follow the proton's path through all the cells of the current
!        block. The current cell indices (i,j,k) and the previous cell indices (ip,jp,kp)
!        will be updated as we move through the block.
!
!
!-------------------- Loop following proton through cells in block --------------------------------------------
!
!
     stepTimeNext = pi_infiniteTime

     do                                ! indefinite loop through the block cells
                                       ! will be broken (exit) by the various conditions
                                       ! of the proton (out of block, out of domain, etc)
!
!
!     ...From the current position, velocity and accelleration of the proton, we determine
!        the parabolic crossing time to the closest cell wall. Then there are 2 possible
!        ways the proton can be traced through the cell:
!
!             1) The parabolic approximation has been found to be valid (low E/B, high speed).
!                In this case we use the parabolic crossing time directly to calculate the
!                new position and velocity of the proton on the cell wall.
!
!             2) The parabolic approximation is not valid, in which case we have to use the
!                confined RK integrator. The parabolic crossing time can be used as a first
!                approximation guess for the RK stepper. Assemble the vector of dependent
!                variables to be passed to the RK stepper. Since this will be a confined RK
!                step, we also pass the minimum/maximum cell wall coordinates. The RK stepper
!                has to stay within the cell boundaries until one of the cell walls is hit.
!                After the confined RK step has been taken, check if the errors are within the
!                error bars.
!
!
        accX = Qm * (Ex + BzInvC * velY - ByInvC * velZ)    ! parabolic acceleration in x-direction
        accY = Qm * (Ey + BxInvC * velZ - BzInvC * velX)    ! parabolic acceleration in y-direction
        accZ = Qm * (Ez + ByInvC * velX - BxInvC * velY)    ! parabolic acceleration in z-direction

        timeX = pi_time2FacesParabolicPath1D (posX, velX, accX, xminCell, xmaxCell, pi_infiniteTime)
        timeY = pi_time2FacesParabolicPath1D (posY, velY, accY, yminCell, ymaxCell, pi_infiniteTime)
        timeZ = pi_time2FacesParabolicPath1D (posZ, velZ, accZ, zminCell, zmaxCell, pi_infiniteTime)

        crossTime = min (timeX, timeY, timeZ)

        if (crossTime == pi_infiniteTime .or. crossTime == 0.0) then
            call Driver_abortFlash ('[pi_traceBlockProtons3DRec] ERROR: infinite/zero cell crossing time')
        end if

        if (pi_useParabolicApproximation) then
            jrkX = Qm * (BzInvC * accY - ByInvC * accZ)     ! jerk in x-direction
            jrkY = Qm * (BxInvC * accZ - BzInvC * accX)     ! jerk in y-direction
            jrkZ = Qm * (ByInvC * accX - BxInvC * accY)     ! jerk in z-direction
            jerkError = sixth * max (abs (jrkX), abs (jrkY), abs (jrkZ)) * crossTime * crossTime * crossTime
            parabolicPath = jerkError <= cellStepError
        else
            parabolicPath = .false.
        end if

        if (parabolicPath) then

            crossTimeHalf = 0.5 * crossTime

            addX = (velX + accX * crossTimeHalf) * crossTime
            addY = (velY + accY * crossTimeHalf) * crossTime
            addZ = (velZ + accZ * crossTimeHalf) * crossTime

            if (pi_screenProtonDiagnostics) then
                parabolicPathLength = pi_parabolicPathLength3D (velX,velY,velZ, &
                                                                accX,accY,accZ, &
                                                                crossTime       )
 
                dgJv = dgJv + addX * CurlBx + addY * CurlBy + addZ * CurlBz
                dgKx = dgKx + Bx * parabolicPathLength
                dgKy = dgKy + By * parabolicPathLength
                dgKz = dgKz + Bz * parabolicPathLength
            end if

            posX = posX + addX
            posY = posY + addY
            posZ = posZ + addZ
            velX = velX + (accX + jrkX * crossTimeHalf) * crossTime
            velY = velY + (accY + jrkY * crossTimeHalf) * crossTime
            velZ = velZ + (accZ + jrkZ * crossTimeHalf) * crossTime
            time = time + crossTime

        else

            stepTimeTry = min (crossTime, stepTimeNext)

            if (stepTimeTry == pi_infiniteTime .or. stepTimeTry == 0.0) then
                call Driver_abortFlash ('[pi_traceBlockProtons3DRec] ERROR: infinite/zero cell stepping time')
            end if

            protonIn (1)  = posX
            protonIn (2)  = posY
            protonIn (3)  = posZ
            protonIn (4)  = velX
            protonIn (5)  = velY
            protonIn (6)  = velZ

            if (pi_screenProtonDiagnostics) then

                protonIn (7)  = dgJv
                protonIn (8)  = dgKx
                protonIn (9)  = dgKy
                protonIn (10) = dgKz

                call  RungeKutta_stepConfined (pi_RungeKuttaMethod,         &
                                               pi_traceODEfunction3DRec,    &
                                               3,                           &  ! # of confined variables
                                               0.0,                         &  ! dummy time -> ODE time independent 
                                               protonIn        (1:10),      &
                                               pi_minConfinement3DRec,      &  ! lower limit confinement function
                                               pi_maxConfinement3DRec,      &  ! upper limit confinement function
                                               protonErrorFrac,             &
                                               protonErrorBase (1:10),      &
                                               stepTimeTry,                 &
                                               stepTimeUsed,                &  ! actual stepping time used
                                               stepTimeNext,                &  ! estimates of next stepping time
                                               protonOut       (1:10),      &
                                               protonError     (1:10)       )  ! can be +ve or -ve

                if (any (abs (protonError (1:10)) > protonErrorFrac * protonErrorBase (1:10))) then
                    call Driver_abortFlash ('[pi_traceBlockProtons3DRec] ERROR: error(s) in RK step out of bounds!')
                end if

                dgJv = protonOut (7)
                dgKx = protonOut (8)
                dgKy = protonOut (9)
                dgKz = protonOut (10)

            else

                call  RungeKutta_stepConfined (pi_RungeKuttaMethod,         &
                                               pi_traceODEfunction3DRec,    &
                                               3,                           &  ! # of confined variables
                                               0.0,                         &  ! dummy time -> ODE time independent 
                                               protonIn        (1:6),       &
                                               pi_minConfinement3DRec,      &  ! lower limit confinement function
                                               pi_maxConfinement3DRec,      &  ! upper limit confinement function
                                               protonErrorFrac,             &
                                               protonErrorBase (1:6),       &
                                               stepTimeTry,                 &
                                               stepTimeUsed,                &  ! actual stepping time used
                                               stepTimeNext,                &  ! estimates of next stepping time
                                               protonOut       (1:6),       &
                                               protonError     (1:6)        )  ! can be +ve or -ve

                if (any (abs (protonError (1:6)) > protonErrorFrac * protonErrorBase (1:6))) then
                    call Driver_abortFlash ('[pi_traceBlockProtons3DRec] ERROR: error(s) in RK step out of bounds!')
                end if

            end if

            posX = protonOut (1)
            posY = protonOut (2)
            posZ = protonOut (3)
            velX = protonOut (4)
            velY = protonOut (5)
            velZ = protonOut (6)
            time = time + stepTimeUsed

        end if
!
!
!     ...If this is an IO proton, write out its current position, if there is sufficient
!        buffer space. Otherwise, do nothing.
!
!
        if (IOproton) then
            IOpointCount = IOpointCount + 1
            if (IOpointCount <= pi_IOmaxPointsPerBlock) then
                pi_IOprotonPointCount (              IOprotonCount       ) = IOpointCount
                pi_IOprotonPoints     (IOpointCount, IOprotonCount, IAXIS) = posX
                pi_IOprotonPoints     (IOpointCount, IOprotonCount, JAXIS) = posY
                pi_IOprotonPoints     (IOpointCount, IOprotonCount, KAXIS) = posZ
            end if
        end if
!
!
!     ...Check where the proton is currently located.
!
!
        outOfCell =     (posX < xminCell - cellWallThicknessHalf) &    ! for debugging purposes
                   .or. (posX > xmaxCell + cellWallThicknessHalf) &    ! will be removed once the
                   .or. (posY < yminCell - cellWallThicknessHalf) &    ! code is running properly
                   .or. (posY > ymaxCell + cellWallThicknessHalf) &
                   .or. (posZ < zminCell - cellWallThicknessHalf) &
                   .or. (posZ > zmaxCell + cellWallThicknessHalf)

        if (outOfCell) then
            call Driver_abortFlash ('[pi_traceBlockProtons3DRec] ERROR: We stepped out of a cell!')
        end if

        newCell =     (posX < xminCell + cellWallThicknessHalf) &      ! checks, if the proton
                 .or. (posX > xmaxCell - cellWallThicknessHalf) &      ! is considered to be on
                 .or. (posY < yminCell + cellWallThicknessHalf) &      ! one of the cell walls
                 .or. (posY > ymaxCell - cellWallThicknessHalf) &
                 .or. (posZ < zminCell + cellWallThicknessHalf) &
                 .or. (posZ > zmaxCell - cellWallThicknessHalf)
!
!
!     ...If, at the current stage, the proton enters a new cell, we have to determine: 1) which new
!        cell (i,j,k) it is and 2) the appropriate nudging values on the proton's position. Due to
!        possible reflective boundary conditions on the block faces, it can happen that the proton
!        stays in the original cell. After handling the logistics inside the following 'if'
!        statement, the new cell indices i,j,k are either the old ones or new ones.
!
!
        if (newCell) then

            dist2minX = abs (xminCell - posX)
            dist2maxX = abs (xmaxCell - posX)
            dist2minY = abs (yminCell - posY)
            dist2maxY = abs (ymaxCell - posY)
            dist2minZ = abs (zminCell - posZ)
            dist2maxZ = abs (zmaxCell - posZ)

            minDistance = min (dist2minX, dist2maxX, &
                               dist2minY, dist2maxY, &
                               dist2minZ, dist2maxZ)

            if (minDistance > cellWallThicknessHalf) then
                call Driver_abortFlash ('[pi_traceBlockProtons3DRec] ERROR: proton to far away from cell face')
            end if

            cellFaceMinX = (dist2minX <= cellWallThicknessHalf)
            cellFaceMaxX = (dist2maxX <= cellWallThicknessHalf)
            cellFaceMinY = (dist2minY <= cellWallThicknessHalf)
            cellFaceMaxY = (dist2maxY <= cellWallThicknessHalf)
            cellFaceMinZ = (dist2minZ <= cellWallThicknessHalf)
            cellFaceMaxZ = (dist2maxZ <= cellWallThicknessHalf)

            velXgt0 = (velX > 0.0)
            velXlt0 = (velX < 0.0)
            velYgt0 = (velY > 0.0)
            velYlt0 = (velY < 0.0)
            velZgt0 = (velZ > 0.0)
            velZlt0 = (velZ < 0.0)

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

                posX   = xminCell
                nudgeX = + cellWallThicknessHalf

                if (velXlt0) then
                    i = i - 1
                    crossX = .true.
                end if

            else if (cellFaceMaxX) then

                posX   = xmaxCell
                nudgeX = - cellWallThicknessHalf

                if (velXgt0) then
                    i = i + 1
                    crossX = .true.
                end if

            end if

            if (cellFaceMinY) then

                posY   = yminCell
                nudgeY = + cellWallThicknessHalf

                if (velYlt0) then
                    j = j - 1
                    crossY = .true.
                end if

            else if (cellFaceMaxY) then

                posY   = ymaxCell
                nudgeY = - cellWallThicknessHalf

                if (velYgt0) then
                    j = j + 1
                    crossY = .true.
                end if

            end if

            if (cellFaceMinZ) then

                posZ   = zminCell
                nudgeZ = + cellWallThicknessHalf

                if (velZlt0) then
                    k = k - 1
                    crossZ = .true.
                end if

            else if (cellFaceMaxZ) then

                posZ   = zmaxCell
                nudgeZ = - cellWallThicknessHalf

                if (velZgt0) then
                    k = k + 1
                    crossZ = .true.
                end if

            end if

            blockFaceMinX = (posX == xminBlock)
            blockFaceMaxX = (posX == xmaxBlock)
            blockFaceMinY = (posY == yminBlock)
            blockFaceMaxY = (posY == ymaxBlock)
            blockFaceMinZ = (posZ == zminBlock)
            blockFaceMaxZ = (posZ == zmaxBlock)

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

            posX = posX + nudgeX
            posY = posY + nudgeY
            posZ = posZ + nudgeZ

            newCell = crossX .or. crossY .or. crossZ

        end if
!
!
!     ...We are now sure about the target cell. Check, if the target cell (i,j,k) is still within the block.
!        If it is, we have to calculate the new electron density and the new electron temperature where
!        the proton is located in the target cell. If the target cell is not within the block, check if the
!        proton coordinates are still within the defined domain. If not, record the proton on its detector
!        screen, (optionally) add its path to the IO proton plot if it hits the screen and mark it as
!        nonexistent. If the proton is still within the domain boundaries, exit the current block loop.
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

                if (pi_opaqueBoundaries .and. pi_cellBoundary (i,j,k) > 0.0) then
                    pi_protons (PROTON_BLCK,proton) = real (NONEXISTENT)
                    exit
                end if

                Ex = pi_cellEfield (1,i,j,k)
                Ey = pi_cellEfield (2,i,j,k)
                Ez = pi_cellEfield (3,i,j,k)

                if (pi_screenProtonDiagnostics) then
                    Bx     = pi_cellBfield (1,i,j,k)
                    By     = pi_cellBfield (2,i,j,k)
                    Bz     = pi_cellBfield (3,i,j,k)
                    BxInvC = Bx * pi_speedOfLightInv
                    ByInvC = By * pi_speedOfLightInv
                    BzInvC = Bz * pi_speedOfLightInv
                    CurlBx = pi_cellCurlBfield (1,i,j,k)
                    CurlBy = pi_cellCurlBfield (2,i,j,k)
                    CurlBz = pi_cellCurlBfield (3,i,j,k)
                else
                    BxInvC = pi_cellBfield (1,i,j,k)
                    ByInvC = pi_cellBfield (2,i,j,k)
                    BzInvC = pi_cellBfield (3,i,j,k)
                end if

                xminCell = pi_cellEdgesX (i  )
                xmaxCell = pi_cellEdgesX (i+1)
                yminCell = pi_cellEdgesY (j  )
                ymaxCell = pi_cellEdgesY (j+1)
                zminCell = pi_cellEdgesZ (k  )
                zmaxCell = pi_cellEdgesZ (k+1)

                stepTimeNext = pi_infiniteTime

            end if

        else

            inDomain =      (posX > pi_xminDomain) &
                      .and. (posX < pi_xmaxDomain) &
                      .and. (posY > pi_yminDomain) &
                      .and. (posY < pi_ymaxDomain) &
                      .and. (posZ > pi_zminDomain) &
                      .and. (posZ < pi_zmaxDomain)

            if (inDomain) then

                pi_protons (PROTON_POSX,proton) = posX
                pi_protons (PROTON_POSY,proton) = posY
                pi_protons (PROTON_POSZ,proton) = posZ
                pi_protons (PROTON_VELX,proton) = velX
                pi_protons (PROTON_VELY,proton) = velY
                pi_protons (PROTON_VELZ,proton) = velZ
                pi_protons (PROTON_TIME,proton) = time
                pi_protons (PROTON_DGJV,proton) = dgJv  ! store this, even if not needed
                pi_protons (PROTON_DGKX,proton) = dgKx  ! store this, even if not needed
                pi_protons (PROTON_DGKY,proton) = dgKy  ! store this, even if not needed
                pi_protons (PROTON_DGKZ,proton) = dgKz  ! store this, even if not needed

            else

                pi_protons (PROTON_BLCK,proton) = real (NONEXISTENT)  ! remove proton from domain

                posX = posX - nudgeX                                  ! undo the nudging
                posY = posY - nudgeY                                  ! since it is not
                posZ = posZ - nudgeZ                                  ! needed anymore

                detector = pi_protons (PROTON_DETC,proton)

                !$omp critical
                      call pi_recordProtonOnScreen (posX, posY, posZ,       &
                                                    velX, velY, velZ,       &
                                                    dgJv, dgKx, dgKy, dgKz, &
                                                    detector,               &
                                                                  onScreen, &
                                                          scrX, scrY, scrZ  )
                !$omp end critical

                if (IOproton .and. pi_IOaddProtonsDomain2Screen .and. onScreen) then
                    IOpointCount = IOpointCount + 1
                    if (IOpointCount <= pi_IOmaxPointsPerBlock) then
                        pi_IOprotonPointCount (              IOprotonCount       ) = IOpointCount
                        pi_IOprotonPoints     (IOpointCount, IOprotonCount, IAXIS) = scrX
                        pi_IOprotonPoints     (IOpointCount, IOprotonCount, JAXIS) = scrY
                        pi_IOprotonPoints     (IOpointCount, IOprotonCount, KAXIS) = scrZ
                    end if
                end if

            end if
     
            exit

        end if
!
!
!-------------------- End loop following proton through cells in block --------------------------------------------
!
!
     end do
!
!
!     ...Check, if we ran out of IO proton buffer space during an IO writeout.
!
!
     if (IOproton .and. (IOpointCount > pi_IOmaxPointsPerBlock)) then
         print *, '[pi_traceBlockProtons3DRec] Proton ',proton,              &
                  ' needs ',IOpointCount,' points per block, but has only ', &
                    pi_IOmaxPointsPerBlock, ' in buffer space'
     end if
!
!
!     ...Consider next proton.
!
!
  end do
!$omp end do nowait
!
!
!     ...Ready!
!
!
  return
end subroutine pi_traceBlockProtons3DRec

!     if (protonTag == 1) then
!         write (*,*) ' i,j,k = ',i,j,k
!         write (*,*) ' posX = ',posX
!         write (*,*) ' posY = ',posY
!         write (*,*) ' posZ = ',posZ
!         write (*,*) ' velX = ',velX
!         write (*,*) ' velY = ',velY
!         write (*,*) ' velZ = ',velZ
!         write (*,*) ' xminCell = ',xminCell
!         write (*,*) ' yminCell = ',yminCell
!         write (*,*) ' zminCell = ',zminCell
!         write (*,*) ' xmaxCell = ',xmaxCell
!         write (*,*) ' ymaxCell = ',ymaxCell
!         write (*,*) ' zmaxCell = ',zmaxCell
!         write (*,*) ' stepTimeTry = ',stepTimeTry
!         write (*,*) ' stepTimeUsed = ',stepTimeUsed
!     end if
