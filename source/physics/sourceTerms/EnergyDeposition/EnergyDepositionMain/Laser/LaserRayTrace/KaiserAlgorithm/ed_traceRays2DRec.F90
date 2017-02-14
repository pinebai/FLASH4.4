!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/KaiserAlgorithm/ed_traceRays2DRec
!!
!! NAME
!!
!!  ed_traceRays2DRec
!!
!! SYNOPSIS
!!
!!  call ed_traceRays2DRec (real, intent (in) :: timeStep)
!!
!! DESCRIPTION
!!
!!  Processes all rays on the collection of blocks on the current processor for those
!!  geometries consisting formally of 2D rectangular grids (cartesian + cylindrical).
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
!!  The use of threading is possible for tracing the rays through each block. 
!!
!!***

subroutine ed_traceRays2DRec (timeStep)

  use EnergyDeposition_data,  ONLY : ed_cellCenters,      &
                                     ed_cellDensity,      &
                                     ed_cellEdges,        &
                                     ed_cellGradNele,     &
                                     ed_cellGradTele,     &
                                     ed_cellNele,         &
                                     ed_cellTele,         &
                                     ed_cellVolume,       &
                                     ed_cellZbar,         &
                                     ed_depoVar,          &
                                     ed_gradOrder,        &
                                     ed_rayBlockID,       &
                                     ed_rayBlockIDCount,  &
                                     ed_rayCount,         &
                                     ed_rayNumberBlockID, &
                                     ed_threadRayTrace

  use Driver_interface,       ONLY : Driver_abortFlash

  use Grid_interface,         ONLY : Grid_getBlkBC,          &
                                     Grid_getBlkBoundBox,    &
                                     Grid_getBlkIndexLimits, &
                                     Grid_getBlkPtr,         &
                                     Grid_getCellCoords,     &
                                     Grid_getDeltas,         &
                                     Grid_releaseBlkPtr

  use ed_interface,           ONLY : ed_blockData2DRec,      &
                                     ed_traceBlockRays2DRec

  implicit none

#include "constants.h"
#include "Flash.h"
#include "EnergyDeposition.h"

  real, intent (in) :: timeStep

  logical :: blockReflectMinX, blockReflectMaxX
  logical :: blockReflectMinY, blockReflectMaxY
  logical :: includeGuardCells
  logical :: notEnoughGuardCells

  integer :: blockID, blockNr
  integer :: iDimGuard, jDimGuard
  integer :: imaxBlock, jmaxBlock
  integer :: imaxData,  jmaxData
  integer :: imaxDerv,  jmaxDerv
  integer :: imaxGuard, jmaxGuard
  integer :: iminBlock, jminBlock
  integer :: iminData,  jminData
  integer :: iminDerv,  jminDerv
  integer :: iminGuard, jminGuard
  integer :: maxDimGuard
  integer :: numberOfRaysInBlock
  integer :: rayFirst, rayLast

  real    :: deltaInvX, deltaInvY
  real    :: deltaX, deltaY
  real    :: xminBlock, yminBlock
  real    :: xmaxBlock, ymaxBlock

  real    :: delta (1:MDIM)

  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)
  integer :: faces       (LOW:HIGH,1:MDIM)
  real    :: bndBox      (LOW:HIGH,1:MDIM)

  real, allocatable :: cellEnergyDepot (:,:)
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
         call Driver_abortFlash ("ed_traceRays2DRec: block ID < 1 encountered! ")
     end if

     call Grid_getDeltas         (blockID, delta)
     call Grid_getBlkBC          (blockID, faces)
     call Grid_getBlkBoundBox    (blockID, bndBox)
     call Grid_getBlkPtr         (blockID, solnData, CENTER)
     call Grid_getBlkIndexLimits (blockID, blkLimits, blkLimitsGC, CENTER)

     deltaX    = delta (IAXIS)
     deltaY    = delta (JAXIS)
     deltaInvX = 1.0 / deltaX
     deltaInvY = 1.0 / deltaY

     xminBlock = bndBox (LOW ,IAXIS)
     xmaxBlock = bndBox (HIGH,IAXIS)
     yminBlock = bndBox (LOW ,JAXIS)
     ymaxBlock = bndBox (HIGH,JAXIS)

     iminBlock = blkLimits (LOW ,IAXIS)
     imaxBlock = blkLimits (HIGH,IAXIS)
     jminBlock = blkLimits (LOW ,JAXIS)
     jmaxBlock = blkLimits (HIGH,JAXIS)

     iminDerv  = iminBlock - 1
     imaxDerv  = imaxBlock + 1
     jminDerv  = jminBlock - 1
     jmaxDerv  = jmaxBlock + 1

     iminData  = iminBlock - ed_gradOrder
     imaxData  = imaxBlock + ed_gradOrder
     jminData  = jminBlock - ed_gradOrder
     jmaxData  = jmaxBlock + ed_gradOrder

     iminGuard = blkLimitsGC (LOW ,IAXIS)
     imaxGuard = blkLimitsGC (HIGH,IAXIS)
     jminGuard = blkLimitsGC (LOW ,JAXIS)
     jmaxGuard = blkLimitsGC (HIGH,JAXIS)

     notEnoughGuardCells =     (iminDerv < iminGuard) &
                          .or. (imaxDerv > imaxGuard) &
                          .or. (jminDerv < jminGuard) &
                          .or. (jmaxDerv > jmaxGuard) &
                          .or. (iminData < iminGuard) &
                          .or. (imaxData > imaxGuard) &
                          .or. (jminData < jminGuard) &
                          .or. (jmaxData > jmaxGuard)

     if (notEnoughGuardCells) then
         call Driver_abortFlash ("ed_traceRays2DRec: Not enough guard cells in 2D")
     end if

     iDimGuard   = imaxGuard - iminGuard + 1
     jDimGuard   = jmaxGuard - jminGuard + 1

     maxDimGuard = max (iDimGuard,jDimGuard)

     allocate (ed_cellCenters     (maxDimGuard  ,2))
     allocate (ed_cellEdges       (maxDimGuard+1,2))
     allocate (ed_cellDensity     (      iminBlock:imaxBlock , jminBlock:jmaxBlock , 1:1))
     allocate (ed_cellVolume      (      iminBlock:imaxBlock , jminBlock:jmaxBlock , 1:1))
     allocate (ed_cellZbar        (      iminBlock:imaxBlock , jminBlock:jmaxBlock , 1:1))
     allocate (ed_cellNele        (      iminData :imaxData  , jminData :jmaxData  , 1:1))
     allocate (ed_cellTele        (      iminData :imaxData  , jminData :jmaxData  , 1:1))
     allocate (ed_cellGradNele    (1:2 , iminDerv :imaxDerv  , jminDerv :jmaxDerv  , 1:1))
     allocate (ed_cellGradTele    (1:2 , iminDerv :imaxDerv  , jminDerv :jmaxDerv  , 1:1))

     includeGuardCells = .true.

     call Grid_getCellCoords (IAXIS, blockID, FACES,  includeGuardCells, ed_cellEdges   (:,1), iDimGuard+1)
     call Grid_getCellCoords (JAXIS, blockID, FACES,  includeGuardCells, ed_cellEdges   (:,2), jDimGuard+1)
     call Grid_getCellCoords (IAXIS, blockID, CENTER, includeGuardCells, ed_cellCenters (:,1), iDimGuard  )
     call Grid_getCellCoords (JAXIS, blockID, CENTER, includeGuardCells, ed_cellCenters (:,2), jDimGuard  )

     call ed_blockData2DRec (blockID,              &
                             iminBlock, imaxBlock, &
                             jminBlock, jmaxBlock, &
                             iminData,  imaxData,  &
                             jminData,  jmaxData,  &
                             iminDerv,  imaxDerv,  &
                             jminDerv,  jmaxDerv,  &
                             deltaX,    deltaY,    &
                             deltaInvX, deltaInvY, &
                             solnData (:,:,:,1)    )

     blockReflectMinX = (faces (LOW ,IAXIS) == REFLECTING)
     blockReflectMaxX = (faces (HIGH,IAXIS) == REFLECTING)
     blockReflectMinY = (faces (LOW ,JAXIS) == REFLECTING)
     blockReflectMaxY = (faces (HIGH,JAXIS) == REFLECTING)
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
     !$omp          jminBlock, jmaxBlock,  &
     !$omp          xminBlock, xmaxBlock,  &
     !$omp          yminBlock, ymaxBlock,  &
     !$omp          deltaX, deltaY,        &
     !$omp          deltaInvX, deltaInvY,  &
     !$omp          blockReflectMinX,      &
     !$omp          blockReflectMaxX,      &
     !$omp          blockReflectMinY,      &
     !$omp          blockReflectMaxY,      &
     !$omp          ed_depoVar,            &
     !$omp          solnData               )
              
          allocate (cellEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock))

          cellEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock) = 0.0

          call ed_traceBlockRays2DRec (timeStep,                          &
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
!
!
!     ...The next section is a critical floating point accumulation from the various threads
!        and so the answer will depend on the number of threads and the thread schedule.
!
!
          !$omp critical (EnergyDepotTransfer)
                solnData  (ed_depoVar,iminBlock:imaxBlock, jminBlock:jmaxBlock, 1) = &
                solnData  (ed_depoVar,iminBlock:imaxBlock, jminBlock:jmaxBlock, 1) + &
                cellEnergyDepot (     iminBlock:imaxBlock, jminBlock:jmaxBlock)
          !$omp end critical (EnergyDepotTransfer)

          deallocate (cellEnergyDepot)

     !$omp end parallel
!
!
!     ...Consider next block.
!
!
     deallocate (ed_cellCenters    )
     deallocate (ed_cellEdges      )
     deallocate (ed_cellDensity    )
     deallocate (ed_cellVolume     )
     deallocate (ed_cellZbar       )
     deallocate (ed_cellNele       )
     deallocate (ed_cellTele       )
     deallocate (ed_cellGradNele   )
     deallocate (ed_cellGradTele   )

     call Grid_releaseBlkPtr (blockID, solnData, CENTER)

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_traceRays2DRec
