!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDisk/pi_readDiskProtons3DRec
!!
!! NAME
!!
!!  pi_readDiskProtons3DRec
!!
!! SYNOPSIS
!!
!!  call pi_readDiskProtons3DRec (integer, intent (in)    :: blockCount,
!!                                integer, intent (in)    :: blockList (:),
!!                                logical, intent (inout) :: moreOnDisk)
!!
!! DESCRIPTION
!!
!!  Reads in a batch of disk protons from the old disk proton file for 3D rectangular
!!  (cartesian) geometries and determines the list of protons that are present in one
!!  of the blocks on the current processor. Their block ID's are not ordered.
!!
!! ARGUMENTS
!!
!!  blockCount : Number of blocks on current processor
!!  blockList  : All block ID numbers
!!  moreOnDisk : if true, there are more disk proton batches on the old disk proton file
!!
!! NOTES
!!
!!***

subroutine pi_readDiskProtons3DRec (blockCount, blockList, moreOnDisk)

  use Driver_interface,    ONLY : Driver_abortFlash
 
  use Grid_interface,      ONLY : Grid_getBlkBC,      &
                                  Grid_getBlkBoundBox
 
  use ProtonImaging_data,  ONLY : pi_cellWallThickness,         &
                                  pi_diskProtonOldFileID,       &
                                  pi_diskProtonOldFileNbatches, &
                                  pi_diskProtons,               &
                                  pi_globalMe,                  &
                                  pi_maxProtonCount,            &
                                  pi_protonCount,               &
                                  pi_protons,                   &
                                  pi_xmaxDomain,                &
                                  pi_ymaxDomain,                &
                                  pi_zmaxDomain
  
  implicit none

#include "ProtonImaging.h"
#include "constants.h"

  integer, intent (in)    :: blockCount
  integer, intent (in)    :: blockList (1:blockCount)
  logical, intent (inout) :: moreOnDisk

  logical :: blockFaceMinX, blockFaceMaxX
  logical :: blockFaceMinY, blockFaceMaxY
  logical :: blockFaceMinZ, blockFaceMaxZ
  logical :: inBlock
  logical :: xmaxLimit, ymaxLimit, zmaxLimit
  logical :: xminLimit, yminLimit, zminLimit

  integer :: block, blockID
  integer :: proton
  integer :: numberOfProtons

  real    :: cellWallThicknessHalf
  real    :: distToFaceMinX, distToFaceMinY, distToFaceMinZ
  real    :: distToFaceMaxX, distToFaceMaxY, distToFaceMaxZ
  real    :: posX, posY, posZ
  real    :: xminBlock, yminBlock, zminBlock
  real    :: xmaxBlock, ymaxBlock, zmaxBlock

  integer :: faces  (LOW:HIGH,1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)

  integer, allocatable :: blockBoundary (:,:)
  real,    allocatable :: blockBndBox   (:,:)
!
!
!     ...Read in the current disk proton batch.
!
!
  if (pi_diskProtonOldFileNbatches <= 0) then
      call Driver_abortFlash ("pi_readDiskProtons: # of disk proton batches on old disk proton file <= 0 !")
  end if

  read (pi_diskProtonOldFileID) numberOfProtons

  if (numberOfProtons > 0) then
      read (pi_diskProtonOldFileID) pi_diskProtons (1:PROTON_ATTRCOUNT,1:numberOfProtons)
      pi_diskProtonOldFileNbatches = pi_diskProtonOldFileNbatches - 1
  end if

  moreOnDisk = (pi_diskProtonOldFileNbatches > 0)

!  write (*,*) ' Proc, # of protons read = ',pi_globalMe,numberOfProtons
!
!
!     ...In order to avoid repeated calls to the grid get bounding box and boundary
!        condition function inside the innermost loop for each proton, we determine the
!        bounding box and boundary consitions for each block beforehand and store
!        this info away for future reference.
!
!
  cellWallThicknessHalf = 0.5 * pi_cellWallThickness

  allocate (blockBndBox   (1:6,1:blockCount))
  allocate (blockBoundary (1:6,1:blockCount))

  do block = 1, blockCount

     blockID = blockList (block)

     call Grid_getBlkBC       (blockID, faces)
     call Grid_getBlkBoundBox (blockID, bndBox)

     blockBndBox   (1,block) = bndBox (LOW ,IAXIS)
     blockBndBox   (2,block) = bndBox (HIGH,IAXIS)
     blockBndBox   (3,block) = bndBox (LOW ,JAXIS)
     blockBndBox   (4,block) = bndBox (HIGH,JAXIS)
     blockBndBox   (5,block) = bndBox (LOW ,KAXIS)
     blockBndBox   (6,block) = bndBox (HIGH,KAXIS)

     blockBoundary (1,block) = faces  (LOW ,IAXIS)
     blockBoundary (2,block) = faces  (HIGH,IAXIS)
     blockBoundary (3,block) = faces  (LOW ,JAXIS)
     blockBoundary (4,block) = faces  (HIGH,JAXIS)
     blockBoundary (5,block) = faces  (LOW ,KAXIS)
     blockBoundary (6,block) = faces  (HIGH,KAXIS)

  end do
!
!
!     ...Outer loop over all disk protons in current batch.
!
!
  do proton = 1, numberOfProtons

     posX = pi_diskProtons (PROTON_POSX, proton)
     posY = pi_diskProtons (PROTON_POSY, proton)
     posZ = pi_diskProtons (PROTON_POSZ, proton)
!
!
!     ...loop over all blocks and see, if the proton is contained in one of it. As soon as
!        that block is found, exit the block loop, since each proton can be assigned to only
!        one block. A proton belongs to a block, if its x,y,z coordinate is such that:
!
!                        x,y,z block lower limit  <=  proton x,y,z  <  x,y,z block upper limit
!
!        except when any of the upper limits of the block coincide with the domain boundaries.
!        In that case the less '<' sign on the right must be replaced by a <= sign, otherwise
!        protons will dissapear and not accounted for, resulting in proton loss.
!
!
     do block = 1, blockCount

        xminBlock = blockBndBox (1,block)
        xmaxBlock = blockBndBox (2,block)
        yminBlock = blockBndBox (3,block)
        ymaxBlock = blockBndBox (4,block)
        zminBlock = blockBndBox (5,block)
        zmaxBlock = blockBndBox (6,block)

        xminLimit = (posX >= xminBlock)
        yminLimit = (posY >= yminBlock)
        zminLimit = (posZ >= zminBlock)
        xmaxLimit = (posX <  xmaxBlock) .or. ((posX == xmaxBlock) .and. (xmaxBlock == pi_xmaxDomain))
        ymaxLimit = (posY <  ymaxBlock) .or. ((posY == ymaxBlock) .and. (ymaxBlock == pi_ymaxDomain))
        zmaxLimit = (posZ <  zmaxBlock) .or. ((posZ == zmaxBlock) .and. (zmaxBlock == pi_zmaxDomain))

        inBlock  =      xminLimit &
                  .and. xmaxLimit &
                  .and. yminLimit &
                  .and. ymaxLimit &
                  .and. zminLimit &
                  .and. zmaxLimit

        if (inBlock) then

            pi_protonCount = pi_protonCount + 1        ! this count is per processor, not block

            if (pi_protonCount > pi_maxProtonCount) then
                call Driver_abortFlash ("pi_readDiskProtons3DRec: No storage left for proton array")
            end if

            blockID = blockList (block)
!
!
!     ...We now need to make sure, that the proton is placed inside the block and not within the block
!        wall. If this is not done, then initial cell index determination can suffer from roundoff
!        errors with placement of the proton in a wrong cell ouside the current block. Note, that even
!        after pushing the proton inside the block is done, there is still the possibility, that the proton
!        sits exactly on a cell face boundary, but this will happen inside the block, so there is no
!        danger of placing the proton initially in the wrong block. The appropriate (initial) nudging for
!        the cell will be done in the corresponding proton tracing routine, as here we do not have the
!        individual cell information.
!
!
            distToFaceMinX = abs (xminBlock - posX)
            distToFaceMaxX = abs (xmaxBlock - posX)
            distToFaceMinY = abs (yminBlock - posY)
            distToFaceMaxY = abs (ymaxBlock - posY)
            distToFaceMinZ = abs (zminBlock - posZ)
            distToFaceMaxZ = abs (zmaxBlock - posZ)

            blockFaceMinX = (distToFaceMinX < cellWallThicknessHalf)
            blockFaceMaxX = (distToFaceMaxX < cellWallThicknessHalf)
            blockFaceMinY = (distToFaceMinY < cellWallThicknessHalf)
            blockFaceMaxY = (distToFaceMaxY < cellWallThicknessHalf)
            blockFaceMinZ = (distToFaceMinZ < cellWallThicknessHalf)
            blockFaceMaxZ = (distToFaceMaxZ < cellWallThicknessHalf)

            if (blockFaceMinX) then
                posX = xminBlock + cellWallThicknessHalf
            else if (blockFaceMaxX) then
                posX = xmaxBlock - cellWallThicknessHalf
            end if

            if (blockFaceMinY) then
                posY = yminBlock + cellWallThicknessHalf
            else if (blockFaceMaxY) then
                posY = ymaxBlock - cellWallThicknessHalf
            end if

            if (blockFaceMinZ) then
                posZ = zminBlock + cellWallThicknessHalf
            else if (blockFaceMaxZ) then
                posZ = zmaxBlock - cellWallThicknessHalf
            end if

            pi_protons (PROTON_POSX,pi_protonCount) = posX
            pi_protons (PROTON_POSY,pi_protonCount) = posY
            pi_protons (PROTON_POSZ,pi_protonCount) = posZ
            pi_protons (PROTON_VELX,pi_protonCount) = pi_diskProtons (PROTON_VELX, proton)
            pi_protons (PROTON_VELY,pi_protonCount) = pi_diskProtons (PROTON_VELY, proton)
            pi_protons (PROTON_VELZ,pi_protonCount) = pi_diskProtons (PROTON_VELZ, proton)
            pi_protons (PROTON_TIME,pi_protonCount) = 0.0
            pi_protons (PROTON_BLCK,pi_protonCount) = real (blockID)
            pi_protons (PROTON_PROC,pi_protonCount) = real (pi_globalMe)
            pi_protons (PROTON_TAGS,pi_protonCount) = pi_diskProtons (PROTON_TAGS, proton)
            pi_protons (PROTON_BEAM,pi_protonCount) = pi_diskProtons (PROTON_BEAM, proton)
            pi_protons (PROTON_DETC,pi_protonCount) = pi_diskProtons (PROTON_DETC, proton)
            pi_protons (PROTON_DGJV,pi_protonCount) = pi_diskProtons (PROTON_DGJV, proton)
            pi_protons (PROTON_DGKX,pi_protonCount) = pi_diskProtons (PROTON_DGKX, proton)
            pi_protons (PROTON_DGKY,pi_protonCount) = pi_diskProtons (PROTON_DGKY, proton)
            pi_protons (PROTON_DGKZ,pi_protonCount) = pi_diskProtons (PROTON_DGKZ, proton)

            exit                                 ! exits the block loop

        end if                                   ! proton in block
     end do                                      ! block loop
  end do                                         ! batch disk protons loop
!
!
!     ...Deallocate all allocated arrays.
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
end subroutine pi_readDiskProtons3DRec
