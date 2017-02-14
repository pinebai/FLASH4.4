!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/PPRT/ed_traceBlockRays2DCyl3D
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
  logical :: impossibleRay
  logical :: inDomain, inBlock
  logical :: newWedge, newWedgeIJ
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

  integer :: i,j,n
  integer :: ip,jp
  integer :: nRayWritePositions
  integer :: numDeadRays
  integer :: ray
  integer :: rayTag
  integer :: rayWriteIndex

  integer :: rayBlock

  real    :: a1X, a1Z
  real    :: accX,accZ
  real    :: c2div2nc
  real    :: dist2minX, dist2maxX
  real    :: dist2minY, dist2maxY
  real    :: dist2minZ, dist2maxZ
  real    :: dX, dZ
  real    :: gaussTime, gaussTimeHalf
  real    :: integral
  real    :: lnLambda
  real    :: minDistance
  real    :: Nele, Tele
  real    :: nu
  real    :: nudgeX, nudgeZ
  real    :: powerLossFactor
  real    :: r0X, r0Z, r1X, r1Z, r2X, r2Z
  real    :: rayCritDens
  real    :: rayPower
  real    :: rayX, rayY, rayZ
  real    :: stepTime, stepTimeHalf, stepTimeFourth
  real    :: v0X, v0Z, v1X, v1Z
  real    :: velX, velY, velZ
  real    :: w,x,y,z
  real    :: wedgeCosine, wedgeSine
  real    :: wedgeDensity
  real    :: wedgeEnergy
  real    :: wedgeMass, wedgeMassInv
  real    :: wedgePower
  real    :: wedgeSlope
  real    :: wedgeStepX, wedgeStepZ
  real    :: wedgeVolume, wedgeVolumeInv
  real    :: wedgeWallThicknessHalf
  real    :: wedgeZbar
  real    :: xminWedge, yminWedge, zminWedge
  real    :: xmaxWedge, ymaxWedge, zmaxWedge

  real    :: NeleDerv (1:3)

  real, parameter :: GaussianRoot (1:2) = (/1.577350269189626,    &  ! is 1 + 1/sqrt(3) for int limits [0,1]
                                            4.226497308103742e-01/)  ! is 1 - 1/sqrt(3) for int limits [0,1]
!
!
!     ...Define some variables.
!
!
  numDeadRays = 0

  wedgeCosine            = ed_laser3Din2DwedgeCosine
  wedgeSine              = ed_laser3Din2DwedgeSine
  wedgeSlope             = ed_laser3Din2DwedgeSlope
  wedgeWallThicknessHalf = ed_cellWallThickness * 0.5
  wedgeStepX             = ed_cellStepTolerance * deltaX
  wedgeStepZ             = ed_cellStepTolerance * deltaZ
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
!     ...Get extra needed info about the initial wedge (i,j).
!
!
     wedgeZbar      = ed_cellZbar    (i,j,1)
     wedgeDensity   = ed_cellDensity (i,j,1)
     wedgeVolume    = ed_cellVolume  (i,j,1)
     wedgeVolumeInv = 1.0 / wedgeVolume
     wedgeMass      = wedgeDensity * wedgeVolume
     wedgeMassInv   = 1.0 / wedgeMass
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
     do                                ! indefinite loop through the block wedges
                                       ! will be broken (exit) by the various conditions
                                       ! of the ray (no power, out of domain, etc)
!
!
!     ...From the current position, velocity and accelleration of the ray, we determine
!        the stepping time, such that the ray's parabolic path in x- and z-direction is
!        accurate to within the demanded position error tolerance. After this section,
!        we know the exact analytic parabolic ray path during the stepping time. We also
!        know that this path is entirely within the current wedge.
!
!
        x = (rayX - xminWedge) * deltaInvX                              ! rescaled [0,1] ray x coordinate
        z = (rayZ - zminWedge) * deltaInvZ                              ! rescaled [0,1] ray z coordinate

        NeleDerv (1:3) = Interpolate_cubic2DFd1 (ed_cellCubicNele (1:16,i,j,1), x,z)

        dX  = NeleDerv (2) * deltaInvX                                  ! d/dx -> d/dX
        dZ  = NeleDerv (3) * deltaInvZ                                  ! d/dz -> d/dZ

        accX = - c2div2nc * dX                                          ! acceleration in x-direction
        accZ = - c2div2nc * dZ                                          ! acceleration in z-direction

        ymaxWedge = wedgeSlope * rayX
        yminWedge = - ymaxWedge

        x = ed_time2FacesParabolicPath1D (rayX, velX, accX, xminWedge, xmaxWedge, ed_infiniteTime)
        z = ed_time2FacesParabolicPath1D (rayZ, velZ, accZ, zminWedge, zmaxWedge, ed_infiniteTime)

        y = ed_time2FacesParabolicPath1D (rayY,                     &
                                          velY + wedgeSlope * velX, &   ! normal vel to yminWedge
                                                 wedgeSlope * accX, &   ! normal acc to yminWedge
                                          yminWedge,                &   ! just one face here
                                          yminWedge,                &   ! (ymin/ymax not parallel)
                                          ed_infiniteTime           )

        w = ed_time2FacesParabolicPath1D (rayY,                     &
                                          velY - wedgeSlope * velX, &   ! normal vel to ymaxWedge
                                               - wedgeSlope * accX, &   ! normal acc to ymaxWedge
                                          ymaxWedge,                &   ! just one face here
                                          ymaxWedge,                &   ! (ymin/ymax not parallel)
                                          ed_infiniteTime           )

        stepTime = min (w,x,y,z)                                        ! time to reach wedge boundary

        if (stepTime == ed_infiniteTime .or. stepTime == 0.0) then
            call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: infinite/zero crossing time for a wedge')
        end if

        stepTimeHalf = 0.5  * stepTime

        r0X = rayX + (velX + accX * stepTimeHalf) * stepTime            ! full time step x-position
        r0Z = rayZ + (velZ + accZ * stepTimeHalf) * stepTime            ! full time step z-position

        v0X = velX + accX * stepTime                                    ! full time step x-velocity
        v0Z = velZ + accZ * stepTime                                    ! full time step z-velocity

        newWedge = .true.

        do                                                              ! check ray path inside wedge

           stepTimeFourth = 0.5 * stepTimeHalf

           r1X = rayX + (velX + accX * stepTimeFourth) * stepTimeHalf   ! 1st half time step x-position
           r1Z = rayZ + (velZ + accZ * stepTimeFourth) * stepTimeHalf   ! 1st half time step z-position

           v1X = velX + accX * stepTimeHalf                             ! 1st half time step x-velocity
           v1Z = velZ + accZ * stepTimeHalf                             ! 1st half time step z-velocity

           x = (r1X - xminWedge) * deltaInvX                            ! rescaled [0,1] ray x coordinate
           z = (r1Z - zminWedge) * deltaInvZ                            ! rescaled [0,1] ray z coordinate

           NeleDerv (1:3) = Interpolate_cubic2DFd1 (ed_cellCubicNele (1:16,i,j,1), x,z)

           a1X = - c2div2nc * NeleDerv (2) * deltaInvX                  ! 1st half time step x-acceleration
           a1Z = - c2div2nc * NeleDerv (3) * deltaInvZ                  ! 1st half time step z-acceleration

           r2X = r1X + (v1X + a1X * stepTimeFourth) * stepTimeHalf      ! 2nd half time step x-position
           r2Z = r1Z + (v1Z + a1Z * stepTimeFourth) * stepTimeHalf      ! 2nd half time step z-position

           badTimeStep =     (abs (r0X - r2X) > wedgeStepX) &           ! decide if time step is good/bad
                        .or. (abs (r0Z - r2Z) > wedgeStepZ)             !

           if (badTimeStep) then

               stepTime     = stepTimeHalf
               stepTimeHalf = stepTimeFourth

               r0X = r1X                                                ! new ray endpoint x-position
               r0Z = r1Z                                                ! new ray endpoint z-position
               v0X = v1X                                                ! new ray endpoint x-velocity
               v0Z = v1Z                                                ! new ray endpoint z-velocity

               newWedge = .false.                                       ! we stay inside the wedge
           else
               exit
           end if

        end do
!
!
!     ...At this stage we have determined the parabolic ray path such that: 1) it is accurate
!        to the demanded path accuracy, 2) it lays entirely within the current wedge and 3) we
!        determined the new end position, either still within the current wedge (newWedge = false)
!        or on one of its walls, ready to exit the wedge (newWedge = true).
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
           z = rayZ + (velZ + accZ * gaussTimeHalf) * gaussTime

           x = (x - xminWedge) * deltaInvX               ! rescaled [0,1] ray x coordinate
           z = (z - zminWedge) * deltaInvZ               ! rescaled [0,1] ray z coordinate

           Nele = Interpolate_cubic2DF (ed_cellCubicNele (1:16,i,j,1), x,z)
           Tele = Interpolate_cubic2DF (ed_cellCubicTele (1:16,i,j,1), x,z)

           if (Nele <= 0.0) then
               call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: Nele <= 0 for a cell')
           end if

           if (Tele <= 0.0) then
               call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: Tele <= 0 for a cell')
           end if

           lnLambda = ed_CoulombFactor (wedgeZbar,         &
                                        ed_electronCharge, &
                                        ed_Boltzmann,      &
                                        Tele,              &
                                        Nele               )

           nu = ed_inverseBremsstrahlungRate (wedgeZbar,         &
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
        wedgePower      = rayPower * (1.0 - powerLossFactor)
        wedgeEnergy     = wedgePower * timeStep

        if (ed_depoVarIsPerMass) then
            wedgeEnergyDepot (i,j) = wedgeEnergyDepot (i,j) + wedgeEnergy * wedgeMassInv
        else
            wedgeEnergyDepot (i,j) = wedgeEnergyDepot (i,j) + wedgeEnergy * wedgeVolumeInv
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
        rayY = rayY + velY * stepTime
        rayZ = r0Z 

        velX = v0X
        velZ = v0Z

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
