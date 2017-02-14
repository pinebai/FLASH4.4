!!****ih* source/physics/TreeRay/TreeRayMain/tr_odInterface
!!
!!***

Module tr_odInterface

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "FortranLangFeatures.fh"

  interface tr_odAccBotNode
     subroutine tr_odAccBotNode(blockno, point, blkLimits, locCoords, solnData, &
          & botnode, accnode)
       implicit none
       integer, intent(IN) :: blockno
       integer, dimension(MDIM), intent(IN) :: point
       integer, dimension(2,MDIM)   :: blkLimits
       real, dimension(:,:,:), intent(IN) :: locCoords
       real, DIMENSION(:,:,:,:), POINTER, intent(IN) :: solnData
       real, dimension(:), intent(IN) :: botnode
       real, dimension(:), intent(INOUT) :: accnode
     end subroutine tr_odAccBotNode
  end interface

  interface tr_odAccNode
     subroutine tr_odAccNode(subnode, accnode)
       implicit none
        real, dimension(:), intent(IN)  :: subnode
        real, dimension(:), intent(OUT) :: accnode
      end subroutine tr_odAccNode
  end interface

  interface tr_odBTInit
     subroutine tr_odBTInit()
     implicit none
     end subroutine tr_odBTInit
  end interface

  interface tr_odFillBotNode
    subroutine tr_odFillBotNode(blockno, point, blkLimits, solnData, botnode)
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER :: solnData
      real, dimension(:), intent(INOUT) :: botnode
    end subroutine tr_odFillBotNode
  end interface

  interface tr_odFinalize
    subroutine tr_odFinalize()
      implicit none
    end subroutine tr_odFinalize
  end interface

  interface tr_odFinalizeCell
    subroutine tr_odFinalizeCell(solnpoint, dl_poc, cdMaps)
      use TreeRay_data, ONLY : tr_ncd, tr_nPix
      implicit none
      real, DIMENSION(:), POINTER, intent(INOUT) :: solnPoint
      real, DIMENSION(MDIM), intent(IN) :: dl_poc
      real, DIMENSION(tr_nCd, 0:tr_nPix-1), intent(IN) :: cdMaps
    end subroutine tr_odFinalizeCell
  end interface

  interface tr_odGetNodeStruct
    subroutine tr_odGetNodeStruct(nsize, bnsize)
      implicit none
      integer, intent(INOUT) :: bnsize, nsize
    end subroutine tr_odGetNodeStruct
  end interface

  interface tr_odInit
    subroutine tr_odInit()
      implicit none
    end subroutine tr_odInit
  end interface

  interface tr_odInitFieldVar
    subroutine tr_odInitFieldVar()
       implicit none
    end subroutine tr_odInitFieldVar
  end interface

  interface tr_odIntegrateRay
    subroutine tr_odIntegrateRay(srcf_ray, rho_ray, erad_ray, eflux)
      use TreeRay_data, ONLY : tr_bhNR
      implicit none
       real,dimension(0:tr_bhNR),intent(INOUT) :: rho_ray, srcf_ray, erad_ray
       real,intent(OUT) :: eflux
     end subroutine tr_odIntegrateRay
  end interface
  
  interface tr_odNodeContrib
    subroutine tr_odNodeContrib(node, trLevel, refLevel, ii, jj, kk, ins, iph, ith, dr)
       implicit none
       real, dimension(:), intent(IN) :: node
       integer, intent(IN) :: trLevel, refLevel, ins, iph, ith
       integer, intent(IN) :: ii, jj, kk
      real, dimension(MDIM+2), intent(IN) :: dr
    end subroutine tr_odNodeContrib
  end interface

  interface tr_odRadToGas
    subroutine tr_odRadToGas(solnPoint, vol_poc, area_poc)
      implicit none
      real, DIMENSION(:), POINTER, intent(INOUT) :: solnPoint
      real, intent(IN) :: vol_poc, area_poc
    end subroutine tr_odRadToGas
  end interface

  interface tr_odSelfContrib
    subroutine tr_odSelfContrib(node, ii,jj,kk)
      implicit none
       real, dimension(:), intent(IN) :: node
       integer, intent(IN) :: ii, jj, kk
    end subroutine tr_odSelfContrib
  end interface

  interface tr_odStartBlock
    subroutine tr_odStartBlock(blockno, blkLimits, solnData)
      integer, intent(IN) :: blockno
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine tr_odStartBlock
  end interface

end Module tr_odInterface

