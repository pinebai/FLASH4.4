!!****h* source/physics/TreeRay/TreeRay_interface
!!
!! This is the header file for the treeray module that defines its
!! public interfaces.
!!***

Module TreeRay_interface

  interface TreeRay_finalize
     subroutine TreeRay_finalize()
     end subroutine TreeRay_finalize
  end interface

  interface TreeRay_init
     subroutine TreeRay_init()
       
     end subroutine TreeRay_init
  end interface

  interface TreeRay_potentialListOfBlocks
     subroutine TreeRay_potentialListOfBlocks(blockCount,blockList)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
     end subroutine TreeRay_potentialListOfBlocks
  end interface

  interface TreeRay_unitTest
     subroutine TreeRay_unitTest( fileUnit, perfect)
       implicit none
       integer, intent(in) :: fileUnit
       logical, intent(out) :: perfect
     end subroutine TreeRay_unitTest
  end interface

! The following subroutines and functions are needed by the tree solver
#include "constants.h"
#include "FortranLangFeatures.fh"
  interface TreeRay_bhAccBotNode
    subroutine TreeRay_bhAccBotNode(blockno, point, blkLimits, locCoords, solnData, botnode, accnode)
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, dimension(:,:,:) :: locCoords
      real, DIMENSION(:,:,:,:), POINTER :: solnData
      real, dimension(:), intent(IN) :: botnode
      real, dimension(:), intent(INOUT) :: accnode
    end subroutine TreeRay_bhAccBotNode
  end interface

  interface TreeRay_bhAccNode
    subroutine TreeRay_bhAccNode(subnode, accnode)
      real, dimension(:), intent(IN)  :: subnode
      real, dimension(:), intent(INOUT) :: accnode
    end subroutine TreeRay_bhAccNode
  end interface

!  interface TreeRay_bhBotNodeContrib
!    subroutine TreeRay_bhBotNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)
!      real, dimension(:), intent(IN) :: node
!      real, intent(IN) :: ndSize
!      real, dimension(MDIM+2), intent(IN) :: dr
!      integer, intent(IN) :: blockno
!      integer, dimension(MDIM), intent(IN) :: point
!      integer, dimension(2,MDIM)   :: blkLimits
!      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
!    end subroutine TreeRay_bhBotNodeContrib
!  end interface

  interface TreeRay_bhBTInit
    subroutine TreeRay_bhBTInit()
    end subroutine TreeRay_bhBTInit
  end interface

  interface
    logical function TreeRay_bhCheckAccuracy(block, point, blkLimits, solnData)
    implicit none
    integer, intent(IN) :: block
    integer, dimension(MDIM), intent(IN) :: point
    integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end function TreeRay_bhCheckAccuracy
  end interface

  interface TreeRay_bhFillBotNode
    subroutine TreeRay_bhFillBotNode(blockno, point, blkLimits, solnData, botnode)
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER :: solnData
      real, dimension(:), intent(INOUT) :: botnode
    end subroutine TreeRay_bhFillBotNode
  end interface

  interface TreeRay_bhFinalizeIter
    subroutine TreeRay_bhFinalizeIter()
    end subroutine TreeRay_bhFinalizeIter
  end interface

  interface TreeRay_bhNodeContrib
    subroutine TreeRay_bhNodeContrib(node, trLevel, reflevel, dr, blockno, point, blkLimits, solnData)
      real, dimension(:), intent(IN) :: node
      integer, intent(IN) :: trLevel, refLevel, blockno
      integer, dimension(MDIM), intent(IN) :: point
      real, dimension(MDIM+2), intent(IN) :: dr
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine TreeRay_bhNodeContrib
  end interface

  interface TreeRay_bhNormalizeNode
    subroutine TreeRay_bhNormalizeNode(smr, node)
      real, dimension(MDIM), intent(IN) :: smr
      real, dimension(:), intent(INOUT) :: node
    end subroutine TreeRay_bhNormalizeNode
  end interface

  interface TreeRay_bhPostprocNode
    subroutine TreeRay_bhPostprocNode(ndSize, node)
      real, intent(IN) :: ndSize
      real, dimension(:), intent(INOUT) :: node
    end subroutine TreeRay_bhPostprocNode
  end interface

  interface TreeRay_bhStartBlock
    subroutine TreeRay_bhStartBlock(blockno, blkLimits, solnData)
      integer, intent(IN) :: blockno
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine TreeRay_bhStartBlock
  end interface

  interface TreeRay_bhTreeWalkEnd
    subroutine TreeRay_bhTreeWalkEnd(iterate)
      logical,intent(OUT) :: iterate
    end subroutine TreeRay_bhTreeWalkEnd
  end interface

  interface TreeRay_bhFinalizeBlock
    subroutine TreeRay_bhFinalizeBlock(blockno, blkLimits, solnData)
      integer, intent(IN) :: blockno
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER :: solnData
    end subroutine TreeRay_bhFinalizeBlock
  end interface


  interface TreeRay_bhInitFieldVar
    subroutine TreeRay_bhInitFieldVar(gpotVar)
      implicit none
      integer, intent(IN) :: gpotVar
    end subroutine TreeRay_bhInitFieldVar
  end interface

  interface TreeRay_bhMAC
    logical function TreeRay_bhMAC(node, ndSize2, dr, blockno, point, blkLimits, solnData)
    implicit none
    real, dimension(:), intent(IN) :: node
    real, intent(IN) :: ndSize2
    integer, dimension(MDIM), intent(IN) :: point
    real, dimension(MDIM+2), intent(IN) :: dr
    integer, intent(IN) :: blockno
    integer, dimension(2,MDIM)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER :: solnData
    end function TreeRay_bhMAC
  end interface
  
  interface TreeRay_bhPartErr
    subroutine TreeRay_bhPartErr(node, ndSize, dr, perr)
    implicit none
    real, dimension(:), intent(IN) :: node
    real, intent(IN) :: ndSize
    real, dimension(MDIM+2), intent(IN) :: dr
    real :: perr
    end subroutine TreeRay_bhPartErr
  end interface
  
  interface TreeRay_bhSelfContrib
    subroutine TreeRay_bhSelfContrib(node, cellsize, blockno, point, blkLimits, solnData)
      implicit none
      real, dimension(:), intent(IN) :: node
      real, dimension(MDIM), intent(IN) :: cellsize
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine TreeRay_bhSelfContrib
  end interface

  interface TreeRay_bhGetNodeStruct
    subroutine TreeRay_bhGetNodeStruct(im, ix, iy, iz, nsize, bnsize)
      implicit none
      integer, intent(IN) :: im, ix, iy, iz
      integer, intent(INOUT) :: bnsize, nsize
    end subroutine TreeRay_bhGetNodeStruct
  end interface




! end of tree solver subroutines/functions



end Module TreeRay_interface
