!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/ed_traceRays1DRec
!!
!! NAME
!!
!!  ed_traceRays1DRec
!!
!! SYNOPSIS
!!
!!  call ed_traceRays1DRec (real, intent (in) :: timeStep)
!!
!! DESCRIPTION
!!
!!  Processes all rays on the collection of blocks on the current processor for those
!!  geometries consisting formally of 1D rectangular grids (cartesian + spherical).
!!  On exit, each ray in each block has either:
!!
!!            i)  reached a different (yet unknown) block
!!           ii)  has been absorbed -> is nonexistent
!!          iii)  has reached the domain boundary and exited -> is nonexistent.
!!
!! ARGUMENTS
!!
!!  timeStep : Current timestep value
!!
!! NOTES
!!
!!  1) This version uses monocubic expansion of the electron number density and electron
!!     temperature inside each of the block cells.
!!
!!  2) The use of threading is possible for tracing the rays through each block. 
!!
!!***

subroutine ed_traceRays1DRec (timeStep)

  use EnergyDeposition_data,  ONLY : ed_cellCubicNele,             &
                                     ed_cellCubicTele,             &
                                     ed_cellDensity,               &
                                     ed_cellEdges,                 &
                                     ed_cellVolume,                &
                                     ed_cellZbar,                  &
                                     ed_depoVar,                   &
                                     ed_gradOrder,                 &
                                     ed_rayBlockID,                &
                                     ed_rayBlockIDCount,           &
                                     ed_rayCount,                  &
                                     ed_rayNumberBlockID,          &
                                     ed_threadRayTrace

  use Driver_interface,       ONLY : Driver_abortFlash

  use Grid_interface,         ONLY : Grid_getBlkBC,          &
                                     Grid_getBlkBoundBox,    &
                                     Grid_getBlkIndexLimits, &
                                     Grid_getBlkPtr,         &
                                     Grid_getCellCoords,     &
                                     Grid_getDeltas,         &
                                     Grid_releaseBlkPtr

  use ed_interface,           ONLY : ed_blockData1DRec,      &
                                     ed_traceBlockRays1DRec

  implicit none

#include "constants.h"
#include "Flash.h"
#include "EnergyDeposition.h"

  real, intent (in) :: timeStep

  logical :: blockReflectMinX, blockReflectMaxX
  logical :: includeGuardCells
  logical :: notEnoughGuardCells

  integer :: blockID, blockNr
  integer :: iDimGuard
  integer :: iminBlock, imaxBlock
  integer :: iminData,  imaxData
  integer :: iminDerv,  imaxDerv
  integer :: iminGuard, imaxGuard
  integer :: numberOfRaysInBlock
  integer :: rayFirst, rayLast

  real    :: deltaX, deltaInvX
  real    :: xminBlock, xmaxBlock

  real    :: delta (1:MDIM)

  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)
  integer :: faces       (LOW:HIGH,1:MDIM)
  real    :: bndBox      (LOW:HIGH,1:MDIM)

  real, allocatable :: cellEnergyDepot (:)
  real, pointer     :: solnData        (:,:,:,:)
!
!
!     ...Loop over all blocks containing active rays. Get all the data associated
!        with the current block.
!
!
  rayLast  = 0

  do blockNr = 1,ed_rayBlockIDCount

     blockID = ed_rayBlockID (blockNr)

     if (blockID < 1) then
         call Driver_abortFlash ("ed_traceRays1DRec: block ID < 1 encountered! ")
     end if

     call Grid_getDeltas         (blockID, delta)
     call Grid_getBlkBC          (blockID, faces)
     call Grid_getBlkBoundBox    (blockID, bndBox)
     call Grid_getBlkPtr         (blockID, solnData, CENTER)
     call Grid_getBlkIndexLimits (blockID, blkLimits, blkLimitsGC, CENTER)

     deltaX    = delta (IAXIS)
     deltaInvX = 1.0 / deltaX

     xminBlock = bndBox (LOW ,IAXIS)
     xmaxBlock = bndBox (HIGH,IAXIS)

     iminBlock = blkLimits (LOW ,IAXIS)
     imaxBlock = blkLimits (HIGH,IAXIS)

     iminDerv  = iminBlock
     imaxDerv  = imaxBlock + 1

     iminData  = iminBlock - 2
     imaxData  = imaxBlock + 2

     iminGuard = blkLimitsGC (LOW ,IAXIS)
     imaxGuard = blkLimitsGC (HIGH,IAXIS)

     notEnoughGuardCells =     (iminData < iminGuard) &
                          .or. (imaxData > imaxGuard)

     if (notEnoughGuardCells) then
         call Driver_abortFlash ("ed_traceRays1DRec: Not enough guard cells in 1D")
     end if

     iDimGuard = imaxGuard - iminGuard + 1

     allocate (ed_cellEdges       (iDimGuard+1,1))
     allocate (ed_cellDensity     (      iminBlock:imaxBlock , 1:1 , 1:1))
     allocate (ed_cellVolume      (      iminBlock:imaxBlock , 1:1 , 1:1))
     allocate (ed_cellZbar        (      iminBlock:imaxBlock , 1:1 , 1:1))
     allocate (ed_cellCubicNele   (1:4,  iminBlock:imaxBlock , 1:1 , 1:1))
     allocate (ed_cellCubicTele   (1:4,  iminBlock:imaxBlock , 1:1 , 1:1))

     includeGuardCells = .true.

     call Grid_getCellCoords (IAXIS, blockID, FACES,  includeGuardCells, ed_cellEdges   (:,1), iDimGuard+1)

     call ed_blockData1DRec (blockID,              &
                             iminBlock, imaxBlock, &
                             iminData,  imaxData,  &
                             iminDerv,  imaxDerv,  &
                             deltaInvX,            &
                             solnData (:,:,1,1)    )

     blockReflectMinX = (faces (LOW ,IAXIS) == REFLECTING)
     blockReflectMaxX = (faces (HIGH,IAXIS) == REFLECTING)
!
!
!     ...Process all rays associated with the current block (threaded version).
!
!
     numberOfRaysInBlock = ed_rayNumberBlockID (blockNr)

     rayFirst = rayLast + 1
     rayLast  = rayLast + numberOfRaysInBlock

     !$omp parallel if (ed_threadRayTrace) &
     !$omp default (none)                  &
     !$omp private (cellEnergyDepot)       &
     !$omp shared  (timeStep,              &
     !$omp          rayFirst,  rayLast,    &
     !$omp          iminBlock, imaxBlock,  &
     !$omp          xminBlock, xmaxBlock,  &
     !$omp          deltaX,                &
     !$omp          deltaInvX,             &
     !$omp          blockReflectMinX,      &
     !$omp          blockReflectMaxX,      &
     !$omp          ed_depoVar,            &
     !$omp          solnData               )
              
          allocate (cellEnergyDepot (iminBlock:imaxBlock))

          cellEnergyDepot (iminBlock:imaxBlock) = 0.0

          call ed_traceBlockRays1DRec (timeStep,                          &
                                       rayFirst,  rayLast,                &
                                       iminBlock, imaxBlock,              &
                                       xminBlock, xmaxBlock,              &
                                       deltaX,                            &
                                       deltaInvX,                         &
                                       blockReflectMinX,                  &
                                       blockReflectMaxX,                  &
                                                          cellEnergyDepot )
!
!
!     ...The next section is a critical floating point accumulation from the various threads
!        and so the answer will depend on the number of threads and the thread schedule.
!
!
          !$omp critical (EnergyDepotTransfer)
                solnData  (ed_depoVar,iminBlock:imaxBlock, 1, 1) = &
                solnData  (ed_depoVar,iminBlock:imaxBlock, 1, 1) + &
                cellEnergyDepot (     iminBlock:imaxBlock)
          !$omp end critical (EnergyDepotTransfer)

          deallocate (cellEnergyDepot)

     !$omp end parallel
!
!
!     ...Consider next block.
!
!
     deallocate (ed_cellCubicNele  )
     deallocate (ed_cellCubicTele  )
     deallocate (ed_cellEdges      )
     deallocate (ed_cellDensity    )
     deallocate (ed_cellVolume     )
     deallocate (ed_cellZbar       )

     call Grid_releaseBlkPtr (blockID, solnData, CENTER)

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_traceRays1DRec
