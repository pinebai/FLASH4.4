!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRaysCreate/ed_createRays1DRec
!!
!! NAME
!!
!!  ed_createRays1DRec
!!
!! SYNOPSIS
!!
!!  call ed_createRays1DRec (integer, intent (in) :: blockCount, 
!!                           integer, intent (in) :: blockList (:), 
!!                           real,    intent (in) :: timeStep,   
!!                           real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates one ray per beam and places it in its initial block for those geometries consisting
!!  formally of 1D rectangular grids (cartesian + spherical). On exit, all rays hitting
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
!!  Since this routine deals with 1D rectangular type of grids, the following unique features
!!  will occur:
!!
!!   1) The number of rays created must be equal to the number of active beams
!!   2) The number of distinct block ID's can be at most 2, corresponding
!!      to the two outer domain blocks.
!!
!!***

subroutine ed_createRays1DRec (blockCount, blockList, timeStep, timeSimulation)

  use Driver_interface,             ONLY : Driver_abortFlash
 
  use Grid_interface,               ONLY : Grid_getBlkBC,      &
                                           Grid_getBlkBoundBox

  use ed_interface,                 ONLY : ed_computeBeamPower

  use EnergyDeposition_data,        ONLY : ed_beams,              &
                                           ed_cellWallThickness,  &
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
                                           ed_xmaxDomain
  
  implicit none

#include "EnergyDeposition.h"
#include "constants.h"

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation

  logical :: beamActive
  logical :: blockFaceMinX, blockFaceMaxX
  logical :: blockReflectMinX, blockReflectMaxX
  logical :: ignoreBoundary
  logical :: inBlock
  logical :: rayReflects
  logical :: xminLimit, xmaxLimit

  integer :: beam
  integer :: blockID
  integer :: block
  integer :: pulseNumber

  real    :: beamCriticalDensity
  real    :: beamPower
  real    :: cellWallThicknessHalf
  real    :: criticalDensityFactor
  real    :: distToFaceMinX, distToFaceMaxX
  real    :: frequency
  real    :: initialRaySpeed
  real    :: lensX
  real    :: rayPower, raySpeed
  real    :: rayX
  real    :: velX
  real    :: xminBlock, xmaxBlock

  integer :: faces  (LOW:HIGH,1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)

  integer, allocatable :: blockBoundary (:,:)
  real,    allocatable :: blockBndBox   (:,:)
!
!
!     ...In order to avoid repeated calls to the grid get bounding box and boundary
!        condition function inside the innermost loop for each ray, we determine the
!        bounding box and boundary conditions for each block beforehand and store
!        this info away for future reference.
!
!
  allocate (blockBndBox   (1:2,1:blockCount))
  allocate (blockBoundary (1:2,1:blockCount))

  do block = 1, blockCount

     blockID = blockList (block)

     call Grid_getBlkBC       (blockID, faces)
     call Grid_getBlkBoundBox (blockID, bndBox)

     blockBndBox   (1,block) = bndBox (LOW ,IAXIS)
     blockBndBox   (2,block) = bndBox (HIGH,IAXIS)

     blockBoundary (1,block) = faces  (LOW ,IAXIS)
     blockBoundary (2,block) = faces  (HIGH,IAXIS)

  end do
!
!
!     ...Set some extra needed data.
!
!
  criticalDensityFactor =  ed_electronMass * PI / (ed_electronCharge * ed_electronCharge)
  cellWallThicknessHalf =  0.5 * ed_cellWallThickness
!
!
!     ...Outer loop over all (still active) beams. A beam is considered active, if its
!        power is above the ray zero power threshold.
!
!
  do beam = 1, ed_numberOfBeams

     pulseNumber = ed_beams (beam) % pulseNumber

     call ed_computeBeamPower (timeStep, timeSimulation, pulseNumber,   beamPower)

     beamActive = (beamPower > ed_rayZeroPower)

     if (beamActive) then

         frequency        = ed_beams (beam) % frequency
         ignoreBoundary   = ed_beams (beam) % ignoreBoundaryCondition
         initialRaySpeed  = ed_beams (beam) % initialRaySpeed
         lensX            = ed_beams (beam) % lensX

         beamCriticalDensity  = criticalDensityFactor * frequency * frequency
!
!
!     ...Create the unique ray for the current beam. Since for 1D grids the beam cross
!        sectional area is a point, there is no beam grid and all complications arising
!        due to a beam's grid are no longer present here. Only the position of the beam's
!        lens decides where the ray hits the domain.
!
!
         if (lensX < ed_xminDomain) then
             rayX = ed_xminDomain
         else if (lensX > ed_xmaxDomain) then
             rayX = ed_xmaxDomain
         else
             call Driver_abortFlash ("ed_createRays1DRec: ray does not hit domain boundary!")
         end if
!
!
!     ...loop over all blocks and see, if ray is contained in one of it. As soon as that block
!        is found, exit the block loop, since each ray can be assigned to only one block.
!        A ray belongs to a block, if its x-coordinate is such that:
!
!                        x block lower limit  <=  ray x  <  x block upper limit
!
!        except when any of the upper limits of the block coincide with the domain boundaries.
!        In that case the less '<' sign on the right must be replaced by a <= sign, otherwise
!        rays will dissapear and not accounted for, resulting in beam power loss.
!
!
         do block = 1, blockCount

            xminBlock = blockBndBox (1,block)
            xmaxBlock = blockBndBox (2,block)

            xminLimit = (rayX >= xminBlock)
            xmaxLimit = (rayX <  xmaxBlock) .or. ((rayX == xmaxBlock) .and. (xmaxBlock == ed_xmaxDomain))

            inBlock   =  xminLimit .and. xmaxLimit

            if (inBlock) then

                distToFaceMinX = abs (xminBlock - rayX)
                distToFaceMaxX = abs (xmaxBlock - rayX)

                blockFaceMinX = (distToFaceMinX < cellWallThicknessHalf)
                blockFaceMaxX = (distToFaceMaxX < cellWallThicknessHalf)

                blockReflectMinX = (blockBoundary (1,block) == REFLECTING)
                blockReflectMaxX = (blockBoundary (2,block) == REFLECTING)

                rayReflects =     (blockFaceMinX .and. blockReflectMinX) &
                             .or. (blockFaceMaxX .and. blockReflectMaxX)

                if (ignoreBoundary .or. .not.rayReflects) then

                    ed_rayCount = ed_rayCount + 1                ! this count is per processor, not block

                    if (ed_rayCount > ed_maxRayCount) then
                        call Driver_abortFlash ("ed_createRays1DRec: Not enough storage for rays array")
                    end if

                    blockID = blockList (block)

                    rayPower = beamPower
                    raySpeed = initialRaySpeed * ed_speedOfLight

                    velX = raySpeed
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

                    ed_rays (RAY_BLCK,ed_rayCount) = real (blockID)
                    ed_rays (RAY_PROC,ed_rayCount) = real (ed_globalMe)
                    ed_rays (RAY_POSX,ed_rayCount) = rayX
                    ed_rays (RAY_POSY,ed_rayCount) = 0.0
                    ed_rays (RAY_POSZ,ed_rayCount) = 0.0
                    ed_rays (RAY_VELX,ed_rayCount) = velX
                    ed_rays (RAY_VELY,ed_rayCount) = 0.0
                    ed_rays (RAY_VELZ,ed_rayCount) = 0.0
                    ed_rays (RAY_POWR,ed_rayCount) = rayPower
                    ed_rays (RAY_DENC,ed_rayCount) = beamCriticalDensity

                    ed_energyInTimestep = ed_energyInTimestep + rayPower * timeStep

                    exit

                end if      ! ray reflection condition
            end if          ! ray in block
         end do             ! block loop
     end if                 ! active beam
  end do                    ! beam loop
!
!
!     ...Deallocate block bounding box and block boundary storage array.
!
!
  deallocate (blockBndBox)
  deallocate (blockBoundary)
!
!
!     ...Ready!
!
!
  return
end subroutine ed_createRays1DRec
