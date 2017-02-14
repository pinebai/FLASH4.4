!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonsTrace/pi_traceProtons3DRec
!!
!! NAME
!!
!!  pi_traceProtons3DRec
!!
!! SYNOPSIS
!!
!!  call pi_traceProtons3DRec ()
!!
!! DESCRIPTION
!!
!!  Processes all protons on the collection of blocks on the current processor for
!!  those geometries consisting formally of 3D rectangular grids (cartesian).
!!  On exit, each proton has either:
!!
!!            i)  reached a different (yet unknown) block
!!           ii)  spent exactly the time step time in the domain -> write it to disk
!!          iii)  has reached the domain boundary and exited.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!  The use of threading is possible for tracing the protons through each block. 
!!
!!***

subroutine pi_traceProtons3DRec ()

  use ProtonImaging_data,  ONLY : pi_cellBfield,              &
                                  pi_cellBoundary,            &
                                  pi_cellCurlBfield,          &
                                  pi_cellEdgesX,              &
                                  pi_cellEdgesY,              &
                                  pi_cellEdgesZ,              &
                                  pi_cellEfield,              &
                                  pi_protonBlockID,           &
                                  pi_protonBlockIDCount,      &
                                  pi_protonCount,             &
                                  pi_protonNumberBlockID,     &
                                  pi_recalculateCellData,     &
                                  pi_screenProtonDiagnostics, &
                                  pi_threadProtonTrace

  use Driver_interface,    ONLY : Driver_abortFlash

  use Grid_interface,      ONLY : Grid_getBlkBC,          &
                                  Grid_getBlkBoundBox,    &
                                  Grid_getBlkIndexLimits, &
                                  Grid_getBlkPtr,         &
                                  Grid_getCellCoords,     &
                                  Grid_getDeltas,         &
                                  Grid_releaseBlkPtr

  use pi_interface,        ONLY : pi_blockData3DRec,        &
                                  pi_calculateCellData,     &
                                  pi_traceBlockProtons3DRec

  implicit none

#include "constants.h"
#include "Flash.h"
#include "ProtonImaging.h"

  logical :: blockReflectMinX, blockReflectMaxX
  logical :: blockReflectMinY, blockReflectMaxY
  logical :: blockReflectMinZ, blockReflectMaxZ
  logical :: includeGuardCells

  integer :: blockID, blockNr
  integer :: iDimBlock, jDimBlock, kDimBlock
  integer :: imaxBlock, jmaxBlock, kmaxBlock
  integer :: iminBlock, jminBlock, kminBlock
  integer :: numberOfProtonsInBlock
  integer :: protonFirst, protonLast

  real    :: deltaX,    deltaY,    deltaZ
  real    :: deltaInvX, deltaInvY, deltaInvZ
  real    :: xminBlock, yminBlock, zminBlock
  real    :: xmaxBlock, ymaxBlock, zmaxBlock

  real    :: delta (1:MDIM)

  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)
  integer :: faces       (LOW:HIGH,1:MDIM)
  real    :: bndBox      (LOW:HIGH,1:MDIM)

  real, pointer :: solnData (:,:,:,:)
!
!
!     ...Loop over all blocks containing active protons. Get all the data associated
!        with the current block.
!
!
  protonLast  = 0

  do blockNr = 1,pi_protonBlockIDCount

     blockID = pi_protonBlockID (blockNr)

     if (blockID < 1) then
         call Driver_abortFlash ("pi_traceProtons3DRec: block ID < 1 encountered! ")
     end if

     call Grid_getDeltas         (blockID, delta)
     call Grid_getBlkBC          (blockID, faces)
     call Grid_getBlkBoundBox    (blockID, bndBox)
     call Grid_getBlkPtr         (blockID, solnData, CENTER)
     call Grid_getBlkIndexLimits (blockID, blkLimits, blkLimitsGC, CENTER)

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

     iDimBlock = imaxBlock - iminBlock + 1
     jDimBlock = jmaxBlock - jminBlock + 1
     kDimBlock = kmaxBlock - kminBlock + 1

     allocate (pi_cellEdgesX     (iminBlock:imaxBlock+1))
     allocate (pi_cellEdgesY     (jminBlock:jmaxBlock+1))
     allocate (pi_cellEdgesZ     (kminBlock:kmaxBlock+1))
     allocate (pi_cellBoundary   (      iminBlock:imaxBlock , jminBlock:jmaxBlock , kminBlock:kmaxBlock))
     allocate (pi_cellBfield     (1:3 , iminBlock:imaxBlock , jminBlock:jmaxBlock , kminBlock:kmaxBlock))
     allocate (pi_cellEfield     (1:3 , iminBlock:imaxBlock , jminBlock:jmaxBlock , kminBlock:kmaxBlock))

     if (pi_screenProtonDiagnostics) then
         allocate (pi_cellCurlBfield (1:3 , iminBlock:imaxBlock , jminBlock:jmaxBlock , kminBlock:kmaxBlock))
     end if

     if (pi_recalculateCellData) then
         call pi_calculateCellData (xminBlock, xmaxBlock, iminBlock, imaxBlock, deltaX, pi_cellEdgesX (:))
         call pi_calculateCellData (yminBlock, ymaxBlock, jminBlock, jmaxBlock, deltaY, pi_cellEdgesY (:))
         call pi_calculateCellData (zminBlock, zmaxBlock, kminBlock, kmaxBlock, deltaZ, pi_cellEdgesZ (:))
     else
         deltaX = delta (IAXIS)
         deltaY = delta (JAXIS)
         deltaZ = delta (KAXIS)
         includeGuardCells = .false.
         call Grid_getCellCoords (IAXIS, blockID, FACES,  includeGuardCells, pi_cellEdgesX (:), iDimBlock+1)
         call Grid_getCellCoords (JAXIS, blockID, FACES,  includeGuardCells, pi_cellEdgesY (:), jDimBlock+1)
         call Grid_getCellCoords (KAXIS, blockID, FACES,  includeGuardCells, pi_cellEdgesZ (:), kDimBlock+1)
     end if

     deltaInvX = 1.0 / deltaX
     deltaInvY = 1.0 / deltaY
     deltaInvZ = 1.0 / deltaZ

     call pi_blockData3DRec (iminBlock, imaxBlock,            &
                             jminBlock, jmaxBlock,            &
                             kminBlock, kmaxBlock,            &
                             deltaInvX, deltaInvY, deltaInvZ, &
                             solnData (:,:,:,:)               )

     blockReflectMinX = (faces (LOW ,IAXIS) == REFLECTING)
     blockReflectMaxX = (faces (HIGH,IAXIS) == REFLECTING)
     blockReflectMinY = (faces (LOW ,JAXIS) == REFLECTING)
     blockReflectMaxY = (faces (HIGH,JAXIS) == REFLECTING)
     blockReflectMinZ = (faces (LOW ,KAXIS) == REFLECTING)
     blockReflectMaxZ = (faces (HIGH,KAXIS) == REFLECTING)
!
!
!     ...Process all protons associated with the current block (threaded version).
!
!
     numberOfProtonsInBlock = pi_protonNumberBlockID (blockNr)

     protonFirst = protonLast + 1
     protonLast  = protonLast + numberOfProtonsInBlock

     !$omp parallel if (pi_threadProtonTrace)        &
     !$omp default (none)                            &
     !$omp shared  (protonFirst, protonLast,         &
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
     !$omp          blockReflectMaxZ                 )

          call pi_traceBlockProtons3DRec (protonFirst, protonLast,           &
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

     !$omp end parallel
!
!
!     ...Consider next block.
!
!
     deallocate (pi_cellEdgesX)
     deallocate (pi_cellEdgesY)
     deallocate (pi_cellEdgesZ)
     deallocate (pi_cellEfield)
     deallocate (pi_cellBfield)
     deallocate (pi_cellBoundary)

     if (allocated (pi_cellCurlBfield)) deallocate (pi_cellCurlBfield)

     call Grid_releaseBlkPtr (blockID, solnData, CENTER)

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine pi_traceProtons3DRec
