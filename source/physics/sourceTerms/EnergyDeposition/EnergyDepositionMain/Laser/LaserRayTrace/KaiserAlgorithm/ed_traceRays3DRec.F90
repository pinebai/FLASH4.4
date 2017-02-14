!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/KaiserAlgorithm/ed_traceRays3DRec
!!
!! NAME
!!
!!  ed_traceRays3DRec
!!
!! SYNOPSIS
!!
!!  call ed_traceRays3DRec (real, intent (in) :: timeStep)
!!
!! DESCRIPTION
!!
!!  Processes all rays on the collection of blocks on the current processor for
!!  those geometries consisting formally of 3D rectangular grids (cartesian).
!!  On exit, each ray has either:
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

subroutine ed_traceRays3DRec (timeStep)

  use EnergyDeposition_data,  ONLY : ed_cellCenters,      &
                                     ed_cellDensity,      &
                                     ed_cellEdges,        &
                                     ed_cellGradNele,     &
                                     ed_cellGradTele,     &
                                     ed_cellNele,         &
                                     ed_cellTele,         &
                                     ed_cellZbar,         &
                                     ed_depoVar,          &
                                     ed_gradOrder,        &
                                     ed_irradVar,         &
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

  use ed_interface,           ONLY : ed_blockData3DRec,      &
                                     ed_traceBlockRays3DRec

  implicit none

#include "constants.h"
#include "Flash.h"
#include "EnergyDeposition.h"

  real, intent (in) :: timeStep

  logical :: blockReflectMinX, blockReflectMaxX
  logical :: blockReflectMinY, blockReflectMaxY
  logical :: blockReflectMinZ, blockReflectMaxZ
  logical :: includeGuardCells
  logical :: notEnoughGuardCells

  integer :: blockID, blockNr
  integer :: iDimGuard, jDimGuard, kDimGuard
  integer :: imaxBlock, jmaxBlock, kmaxBlock
  integer :: imaxData,  jmaxData,  kmaxData
  integer :: imaxDerv,  jmaxDerv,  kmaxDerv
  integer :: imaxGuard, jmaxGuard, kmaxGuard
  integer :: iminBlock, jminBlock, kminBlock
  integer :: iminData,  jminData,  kminData
  integer :: iminDerv,  jminDerv,  kminDerv
  integer :: iminGuard, jminGuard, kminGuard
  integer :: maxDimGuard
  integer :: numberOfRaysInBlock
  integer :: rayFirst, rayLast

  real    :: deltaX,    deltaY,    deltaZ
  real    :: deltaInvX, deltaInvY, deltaInvZ
  real    :: xminBlock, yminBlock, zminBlock
  real    :: xmaxBlock, ymaxBlock, zmaxBlock

  real    :: delta (1:MDIM)

  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)
  integer :: faces       (LOW:HIGH,1:MDIM)
  real    :: bndBox      (LOW:HIGH,1:MDIM)

  real, allocatable :: cellEnergyDepot (:,:,:)
  real, allocatable :: cellIntensityDepot (:,:,:)
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
         call Driver_abortFlash ("ed_traceRays3DRec: block ID < 1 encountered! ")
     end if

     call Grid_getDeltas         (blockID, delta)
     call Grid_getBlkBC          (blockID, faces)
     call Grid_getBlkBoundBox    (blockID, bndBox)
     call Grid_getBlkPtr         (blockID, solnData, CENTER)
     call Grid_getBlkIndexLimits (blockID, blkLimits, blkLimitsGC, CENTER)

     deltaX    = delta (IAXIS)
     deltaY    = delta (JAXIS)
     deltaZ    = delta (KAXIS)
     deltaInvX = 1.0 / deltaX
     deltaInvY = 1.0 / deltaY
     deltaInvZ = 1.0 / deltaZ

     xminBlock = bndBox (LOW ,IAXIS)
     xmaxBlock = bndBox (HIGH,IAXIS)
     yminBlock = bndBox (LOW ,JAXIS)
     ymaxBlock = bndBox (HIGH,JAXIS)
     zminBlock = bndBox (LOW ,KAXIS)
     zmaxBlock = bndBox (HIGH,KAXIS)

     iminBlock = blkLimits (LOW ,IAXIS)
     imaxBlock = blkLimits (HIGH,IAXIS)
     jminBlock = blkLimits (LOW ,JAXIS)
     jmaxBlock = blkLimits (HIGH,JAXIS)
     kminBlock = blkLimits (LOW ,KAXIS)
     kmaxBlock = blkLimits (HIGH,KAXIS)

     iminDerv  = iminBlock - 1
     imaxDerv  = imaxBlock + 1
     jminDerv  = jminBlock - 1
     jmaxDerv  = jmaxBlock + 1
     kminDerv  = kminBlock - 1
     kmaxDerv  = kmaxBlock + 1

     iminData  = iminBlock - ed_gradOrder
     imaxData  = imaxBlock + ed_gradOrder
     jminData  = jminBlock - ed_gradOrder
     jmaxData  = jmaxBlock + ed_gradOrder
     kminData  = kminBlock - ed_gradOrder
     kmaxData  = kmaxBlock + ed_gradOrder

     iminGuard = blkLimitsGC (LOW ,IAXIS)
     imaxGuard = blkLimitsGC (HIGH,IAXIS)
     jminGuard = blkLimitsGC (LOW ,JAXIS)
     jmaxGuard = blkLimitsGC (HIGH,JAXIS)
     kminGuard = blkLimitsGC (LOW ,KAXIS)
     kmaxGuard = blkLimitsGC (HIGH,KAXIS)

     notEnoughGuardCells =     (iminDerv < iminGuard) &
                          .or. (imaxDerv > imaxGuard) &
                          .or. (jminDerv < jminGuard) &
                          .or. (jmaxDerv > jmaxGuard) &
                          .or. (kminDerv < kminGuard) &
                          .or. (kmaxDerv > kmaxGuard) &
                          .or. (iminData < iminGuard) &
                          .or. (imaxData > imaxGuard) &
                          .or. (jminData < jminGuard) &
                          .or. (jmaxData > jmaxGuard) &
                          .or. (kminData < kminGuard) &
                          .or. (kmaxData > kmaxGuard)

     if (notEnoughGuardCells) then
         call Driver_abortFlash ("ed_traceRays3DRec: Not enough guard cells in 3D")
     end if

     iDimGuard   = imaxGuard - iminGuard + 1
     jDimGuard   = jmaxGuard - jminGuard + 1
     kDimGuard   = kmaxGuard - kminGuard + 1

     maxDimGuard = max (iDimGuard, jDimGuard, kDimGuard)

     allocate (ed_cellCenters     (maxDimGuard  ,3))
     allocate (ed_cellEdges       (maxDimGuard+1,3))
     allocate (ed_cellDensity     (      iminBlock:imaxBlock , jminBlock:jmaxBlock , kminBlock:kmaxBlock))
     allocate (ed_cellZbar        (      iminBlock:imaxBlock , jminBlock:jmaxBlock , kminBlock:kmaxBlock))
     allocate (ed_cellNele        (      iminData :imaxData  , jminData :jmaxData  , kminData :kmaxData ))
     allocate (ed_cellTele        (      iminData :imaxData  , jminData :jmaxData  , kminData :kmaxData ))
     allocate (ed_cellGradNele    (1:3 , iminDerv :imaxDerv  , jminDerv :jmaxDerv  , kminDerv :kmaxDerv ))
     allocate (ed_cellGradTele    (1:3 , iminDerv :imaxDerv  , jminDerv :jmaxDerv  , kminDerv :kmaxDerv ))

     includeGuardCells = .true.

     call Grid_getCellCoords (IAXIS, blockID, FACES,  includeGuardCells, ed_cellEdges   (:,1), iDimGuard+1)
     call Grid_getCellCoords (JAXIS, blockID, FACES,  includeGuardCells, ed_cellEdges   (:,2), jDimGuard+1)
     call Grid_getCellCoords (KAXIS, blockID, FACES,  includeGuardCells, ed_cellEdges   (:,3), kDimGuard+1)
     call Grid_getCellCoords (IAXIS, blockID, CENTER, includeGuardCells, ed_cellCenters (:,1), iDimGuard  )
     call Grid_getCellCoords (JAXIS, blockID, CENTER, includeGuardCells, ed_cellCenters (:,2), jDimGuard  )
     call Grid_getCellCoords (KAXIS, blockID, CENTER, includeGuardCells, ed_cellCenters (:,3), kDimGuard  )

     call ed_blockData3DRec (iminBlock, imaxBlock,            &
                             jminBlock, jmaxBlock,            &
                             kminBlock, kmaxBlock,            &
                             iminData,  imaxData,             &
                             jminData,  jmaxData,             &
                             kminData,  kmaxData,             &
                             iminDerv,  imaxDerv,             &
                             jminDerv,  jmaxDerv,             &
                             kminDerv,  kmaxDerv,             &
                             deltaX,    deltaY,    deltaZ,    &
                             deltaInvX, deltaInvY, deltaInvZ, &
                             solnData (:,:,:,:)    )

     blockReflectMinX = (faces (LOW ,IAXIS) == REFLECTING)
     blockReflectMaxX = (faces (HIGH,IAXIS) == REFLECTING)
     blockReflectMinY = (faces (LOW ,JAXIS) == REFLECTING)
     blockReflectMaxY = (faces (HIGH,JAXIS) == REFLECTING)
     blockReflectMinZ = (faces (LOW ,KAXIS) == REFLECTING)
     blockReflectMaxZ = (faces (HIGH,KAXIS) == REFLECTING)
!
!
!     ...Process all rays associated with the current block (threaded version).
!
!
     numberOfRaysInBlock = ed_rayNumberBlockID (blockNr)

     rayFirst = rayLast + 1
     rayLast  = rayLast + numberOfRaysInBlock

     !$omp parallel if (ed_threadRayTrace)           &
     !$omp default (none)                            &
     !$omp private (cellEnergyDepot,                 &
     !$omp          cellIntensityDepot)              &
     !$omp shared  (timeStep,                        &
     !$omp          rayFirst,  rayLast,              &
     !$omp          iminBlock, imaxBlock,            &
     !$omp          jminBlock, jmaxBlock,            &
     !$omp          kminBlock, kmaxBlock,            &
     !$omp          xminBlock, xmaxBlock,            &
     !$omp          yminBlock, ymaxBlock,            &
     !$omp          zminBlock, zmaxBlock,            &
     !$omp          deltaX, deltaY, deltaZ,          &
     !$omp          deltaInvX, deltaInvY, deltaInvZ, &
     !$omp          blockReflectMinX,                &
     !$omp          blockReflectMaxX,                &
     !$omp          blockReflectMinY,                &
     !$omp          blockReflectMaxY,                &
     !$omp          blockReflectMinZ,                &
     !$omp          blockReflectMaxZ,                &
     !$omp          ed_depoVar,                      &
     !$omp          ed_irradVar,                     &
     !$omp          solnData                         )
              
          allocate (cellEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock,kminBlock:kmaxBlock))

          cellEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock,kminBlock:kmaxBlock) = 0.0

          if (ed_irradVar > 0) then
             allocate (cellIntensityDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock,kminBlock:kmaxBlock))

             cellIntensityDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock,kminBlock:kmaxBlock) = 0.0

             call ed_traceBlockRays3DRec (timeStep,                          &
                                       rayFirst,  rayLast,                &
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
                                       blockReflectMaxZ,                  &
                                                         cellEnergyDepot, &
                                                       cellIntensityDepot )

          else

             call ed_traceBlockRays3DRec (timeStep,                          &
                                       rayFirst,  rayLast,                &
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
                                       blockReflectMaxZ,                  &
                                                          cellEnergyDepot )

          end if
!
!
!     ...The next section is a critical floating point accumulation from the various threads
!        and so the answer will depend on the number of threads and the thread schedule.
!
!
          !$omp critical (EnergyDepotTransfer)

                solnData  (ed_depoVar,iminBlock:imaxBlock, jminBlock:jmaxBlock, kminBlock:kmaxBlock) = &
                solnData  (ed_depoVar,iminBlock:imaxBlock, jminBlock:jmaxBlock, kminBlock:kmaxBlock) + &
                cellEnergyDepot      (iminBlock:imaxBlock, jminBlock:jmaxBlock, kminBlock:kmaxBlock)

          !$omp end critical (EnergyDepotTransfer)

          deallocate (cellEnergyDepot)

          if (ed_irradVar > 0) then

             !$omp critical (IntensityDepotTransfer)
                solnData  (ed_irradVar,iminBlock:imaxBlock, jminBlock:jmaxBlock, kminBlock:kmaxBlock) = &
                solnData  (ed_irradVar,iminBlock:imaxBlock, jminBlock:jmaxBlock, kminBlock:kmaxBlock) + &
                cellIntensityDepot    (iminBlock:imaxBlock, jminBlock:jmaxBlock, kminBlock:kmaxBlock)
             !$omp end critical (IntensityDepotTransfer)

             deallocate (cellIntensityDepot)
          end if


     !$omp end parallel
!
!
!     ...Consider next block.
!
!
     deallocate (ed_cellCenters    )
     deallocate (ed_cellEdges      )
     deallocate (ed_cellDensity    )
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
end subroutine ed_traceRays3DRec
