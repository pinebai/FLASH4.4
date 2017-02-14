!!****ih* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhLocalInterface
!!
!! This is the header file for the Barnes-Hut tree solver that defines
!! additional interfaces private to the BHTree implementation.
!!***
Module gr_bhLocalInterface

#include "constants.h"
#include "FortranLangFeatures.fh"

  interface
    subroutine gr_bhBuildTreeBlock(block)
    implicit none
    integer,intent(in) :: block
    end subroutine gr_bhBuildTreeBlock
  end interface


  interface
    subroutine gr_bhCalcBlockTreePos()
    implicit none
    end subroutine gr_bhCalcBlockTreePos
  end interface


  interface
    subroutine gr_bhComParentTree(level)
    implicit none
    integer, intent(IN) :: level
    end subroutine gr_bhComParentTree
  end interface


  interface
    subroutine gr_bhDestroyTree()
    implicit none
    end subroutine gr_bhDestroyTree
  end interface

  interface
    subroutine gr_bhPeriodicDr(dr)
    implicit none
    real, dimension(MDIM+2), intent(INOUT) :: dr
    end subroutine gr_bhPeriodicDr
  end interface


  interface
    subroutine gr_bhFinalize()
    implicit none
    end subroutine gr_bhFinalize
  end interface

  interface
    subroutine gr_bhFinalizeIter()
    implicit none
    end subroutine gr_bhFinalizeIter
  end interface


  interface
    integer function gr_bhGetTreePos(level, mi)
      use gr_bhData, ONLY: gr_bhTreeLevels
      implicit none
      integer,intent(in) :: level, mi(1:gr_bhTreeLevels)
    end function gr_bhGetTreePos
  end interface


  interface
    integer function gr_bhGetTreeSize(level)
    implicit none
    integer,intent(in) :: level
    end function gr_bhGetTreeSize
  end interface

  interface
    subroutine gr_bhLeafContrib(block, tr, cpu, solnData, blkLimits, cellcnt, nodecnt)
      implicit none
      integer,intent(in) :: block, tr, cpu
      real,intent(out)   :: cellcnt, nodecnt
      real, dimension(:,:,:,:), POINTER_INTENT_IN :: solnData
      integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    end subroutine gr_bhLeafContrib
  end interface

  interface
    subroutine gr_bhParentContrib(block, tr, cpu, solnData, blkLimits)
      implicit none
      integer,intent(in) :: block, tr, cpu
      real, dimension(:,:,:,:), POINTER_INTENT_IN :: solnData
      integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    end subroutine gr_bhParentContrib
  end interface


  interface
    subroutine gr_bhTreeWalkBlock(block)
    implicit none
    integer, intent(IN) :: block
    end subroutine gr_bhTreeWalkBlock
  end interface

  interface
    subroutine gr_bhTreeWalkBlockUnified(block)
    implicit none
    integer, intent(IN) :: block
    end subroutine gr_bhTreeWalkBlockUnified
  end interface

  interface
    subroutine gr_bhTreeWalkBlockPQ(block)
    implicit none
    integer, intent(IN) :: block
    end subroutine gr_bhTreeWalkBlockPQ
  end interface

  interface
    subroutine gr_bhTreeWalkPoint(x, y, z, block, point, blkLimits, solnData, dcount)
    implicit none
    real, intent(IN) :: x, y, z
    integer, intent(IN)    :: block
    integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    integer, dimension(MDIM), intent(IN) :: point
    real, POINTER, DIMENSION(:,:,:,:) :: solnData
    real, intent(INOUT) :: dcount(1:3)
    end subroutine gr_bhTreeWalkPoint
  end interface

!! Subroutines calling equivalent subroutines of physical modules that solve a specific problem using
!! the tree code

  interface
    subroutine gr_bhAccBotNode(blockno, point, blkLimits, locCoords, solnData, botnode, accnode)
    implicit none
    integer, intent(IN) :: blockno
    integer, dimension(MDIM), intent(IN) :: point
    integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    real, dimension(:,:,:), intent(IN) :: locCoords
    real, DIMENSION(:,:,:,:), POINTER :: solnData
    real, dimension(:), intent(IN) :: botnode
    real, dimension(:), intent(OUT) :: accnode
    end subroutine gr_bhAccBotNode
  end interface

  interface
    subroutine gr_bhAccNode(subnode, accnode)
    implicit none
    real, dimension(:), intent(IN)  :: subnode
    real, dimension(:), intent(OUT) :: accnode
    end subroutine gr_bhAccNode
  end interface

  interface
    subroutine gr_bhBotNodeContrib(node, trLevel, refLevel, dr, blockno, point, blkLimits, solnData)
    implicit none
    real, dimension(:), intent(IN) :: node
    integer, intent(IN) :: trLevel, refLevel
    real, dimension(MDIM+2), intent(IN) :: dr
    integer, intent(IN) :: blockno
    integer, dimension(MDIM), intent(IN) :: point
    integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine gr_bhBotNodeContrib
  end interface

  interface
    subroutine gr_bhBTInit()
    implicit none
    end subroutine gr_bhBTInit
  end interface

  interface
    logical function gr_bhCheckAccuracy(block, point, blkLimits, solnData)
    implicit none
    integer, intent(IN) :: block
    integer, dimension(MDIM), intent(IN) :: point
    integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end function gr_bhCheckAccuracy
  end interface

  interface
    subroutine gr_bhFillBotNode(blockno, point, blkLimits, solnData, botnode)
    implicit none
    integer, intent(IN) :: blockno
    integer, dimension(MDIM), intent(IN) :: point
    integer, dimension(2,MDIM)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    real, dimension(:), intent(OUT) :: botnode
    end subroutine gr_bhFillBotNode
  end interface

  interface
    subroutine gr_bhFinalizeBlock(blockno, blkLimits, solnData)
    implicit none
    integer, intent(IN) :: blockno
    integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER :: solnData
    end subroutine gr_bhFinalizeBlock
  end interface

  interface
    logical function gr_bhMAC(physMAC, node, ndSize2, drGC, dr, &
    & blockno, point, blkLimits, solnData)
    implicit none
    logical, intent(IN) :: physMAC
    real, dimension(:), intent(IN) :: node
    real, intent(IN) :: ndSize2
    integer, dimension(MDIM), intent(IN) :: point
    real, dimension(MDIM), intent(IN) :: drGC
    real, dimension(MDIM+2), intent(IN) :: dr
    integer, intent(IN) :: blockno
    integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end function gr_bhMAC
  end interface

  interface
    subroutine gr_bhNodeContrib(node, trLevel, refLevel, dr, blockno, point, blkLimits, solnData)
    implicit none
    real, dimension(:), intent(IN) :: node
    integer, intent(IN) :: trLevel, refLevel
    integer, intent(IN) :: blockno
    integer, dimension(MDIM), intent(IN) :: point
    integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    real, dimension(MDIM+2), intent(IN) :: dr
    real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine gr_bhNodeContrib
  end interface

  interface
    subroutine gr_bhNormalizeNode(node)
    implicit none
    real, dimension(:), intent(OUT) :: node
    end subroutine gr_bhNormalizeNode
  end interface

  interface
    subroutine gr_bhPartErr(node, ndSize, dr, perr)
    implicit none
    real, dimension(:), intent(IN) :: node
    real, intent(IN) :: ndSize
    real, dimension(MDIM+2), intent(IN) :: dr
    real, dimension(2), intent(OUT) :: perr
    end subroutine gr_bhPartErr
  end interface

  interface
    subroutine gr_bhPostprocNode(ndSize, node)
    implicit none
    real, intent(IN) :: ndSize
    real, dimension(:), intent(INOUT) :: node
    end subroutine gr_bhPostprocNode
  end interface

  interface
    subroutine gr_bhPQInsert(block, tr, cpu, btp, point_mbLp1, phys_point, node)
    implicit none
    integer, intent(IN)               :: block, tr, cpu, btp
    integer, dimension(MDIM), intent(IN) :: point_mbLp1
    real, dimension(MDIM), intent(IN) :: phys_point
    real, dimension(:), intent(INOUT)  :: node
    end subroutine gr_bhPQInsert
  end interface

  interface
    subroutine gr_bhPQPull(tr, cpu, btp, int_type, dr, perr)
    implicit none
    integer, intent(OUT)                :: tr, cpu, btp, int_type
    real, dimension(MDIM+2),intent(OUT) :: dr
    real, intent(OUT)                   :: perr
    end subroutine gr_bhPQPull
  end interface

  interface
    subroutine gr_bhSelfContrib(node, cellsize, blockno, point, blkLimits, solnData)
    implicit none
    real, dimension(:), intent(IN) :: node
    real, dimension(MDIM), intent(IN) :: cellsize
    integer, intent(IN) :: blockno
    integer, dimension(MDIM), intent(IN) :: point
    integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine gr_bhSelfContrib
  end interface

  interface
    subroutine gr_bhStartBlock(blockno, blkLimits, solnData)
    implicit none
    integer, intent(IN) :: blockno
    integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine gr_bhStartBlock
  end interface
end Module gr_bhLocalInterface


