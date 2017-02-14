!!****ih* source/Grid/localAPI/gr_interface
!!
!! NAME
!!
!!  gr_interface
!!
!! SYNOPSIS
!!
!!  use gr_interface
!!
!! DESCRIPTION
!!
!!  Interfaces for some subprograms private to the GridMain subunit.
!!
!!  Currently, these are mostly subprograms used in the PARAMESH3/4 Grid implementations.
!!
!!***

! Modification history:
!     Created gr_flashWithPM3_interfaces          February 2007  KW
!     Renamed gr_pm3Interface                     January  2008  AD
!     Renamed gr_interface and moved to localAPI  June     2008  KW
!     Added gr_findMean, alphabetized functions   July     2008  LBR
!     Added gr_findWhichChildren,gr_findAllNeghID, 
!           gr_checkGridState                     Nov      2008  CD
!     Added gr_updateRefinement                   December 2008  KW
!     Modified gr_setGcFillNLayers                June     2009  KW
!     Added gr_createBlock                        October  2012  KW

module gr_interface
#include "constants.h"
#include "Flash.h"
  implicit none

  interface
     subroutine gr_createBlock(blockImin, blockImax, &
          blockJmin, blockJmax, blockKmin, blockKmax,blockID)
       implicit none
       real,intent(IN) :: blockImin, blockImax, blockJmin, blockJmax, blockKmin, blockKmax
       integer,intent(IN) :: blockID
     end subroutine gr_createBlock
  end interface

  interface
     integer function gr_extractBCForDirection(packedBCs,axis,leftOrRight)
     ! implementation in GridMain/paramesh/paramesh4/gr_packBCs.F90
       implicit none
       integer,intent(IN) :: packedBCs,axis,leftOrRight
     end function gr_extractBCForDirection
  end interface

  interface
     subroutine gr_findBlock(blkList,blkCount,pos,blockID)
       implicit none
       integer,intent(IN) :: blkCount
       integer,dimension(blkCount),intent(IN) :: blkList
       real,dimension(MDIM),intent(IN) :: pos
       integer,intent(INOUT) :: blockID
     end subroutine gr_findBlock
  end interface

  interface
     subroutine gr_findMean(iSrc, iType, bGuardcell, mean)
       implicit none
       integer, intent(in) :: iSrc, iType
       logical, intent(in) :: bGuardcell
       real, intent(out) :: mean
     end subroutine gr_findMean
  end interface

  interface
     subroutine gr_findWhichChild(pos,bndBox,negh,whichChild)
       implicit none
       real,dimension(MDIM), intent(IN) :: pos
       real,dimension(LOW:HIGH,MDIM),intent(IN) :: bndBox
       integer, dimension(MDIM),intent(IN) :: negh
       integer, intent(OUT) :: whichChild
     end subroutine gr_findWhichChild
  end interface

  interface
     subroutine gr_findNeghID(blockID,pos,negh,neghID)
       implicit none
       integer,intent(IN) :: blockID
       real,dimension(MDIM), intent(IN) :: pos
       integer,dimension(MDIM),intent(IN) :: negh
       integer,dimension(BLKNO:PROCNO),intent(OUT) :: neghID
     end subroutine gr_findNeghID
  end interface

  interface
     subroutine gr_getBlkHandle(remoteBlockID, pe, blockHandle)
     ! implementation in GridMain/paramesh
       implicit none
       integer, intent(in) :: remoteBlockID, pe
       integer, intent(INOUT) :: blockHandle
     end subroutine gr_getBlkHandle
  end interface

  interface
     integer function gr_packBCs(bcILeft, bcIRight, bcJLeft, bcJRight, bcKLeft, bcKRight)
     ! implementation in GridMain/paramesh/paramesh4
       implicit none
       integer,intent(IN) :: bcILeft, bcIRight, bcJLeft, bcJRight, bcKLeft, bcKRight
     end function gr_packBCs
  end interface
  
  interface
     subroutine gr_setGcFillNLayers(layers, idir, guard, minLayers, returnLayers)
     ! implementation in GridMain/paramesh
       implicit none
       integer,dimension(MDIM), intent(OUT) :: layers
       integer, intent(IN)  :: idir, guard
       integer, intent(IN),OPTIONAL  :: minLayers
       integer,intent(OUT),OPTIONAL  :: returnLayers(MDIM)
     end subroutine gr_setGcFillNLayers
  end interface

  interface
     subroutine gr_neghAtSameLevel(blockID,atSameLevel)
       integer,intent(IN) :: blockID
       logical,dimension(LEFT_EDGE:RIGHT_EDGE,&
            LEFT_EDGE:RIGHT_EDGE,&
            LEFT_EDGE:RIGHT_EDGE),intent(OUT) :: atSameLevel
       
     end subroutine gr_neghAtSameLevel
  end interface

  interface
     subroutine gr_sanitizeDataAfterInterp(blkList,count, info, layers)
       integer,intent(IN) :: count
       integer, dimension(count), intent(IN) :: blkList
       character(len=*), intent(IN) :: info
       integer, dimension(MDIM), intent(IN):: layers
     end subroutine gr_sanitizeDataAfterInterp
  end interface


  interface
     subroutine gr_findWhichChildren(numNegh,Negh,whichChildren)
       integer,intent(IN) :: numNegh
       integer, dimension(MDIM),intent(IN) :: Negh
       integer, intent(OUT) :: whichChildren(numNegh)
     end subroutine gr_findWhichChildren
  end interface


  interface
     subroutine gr_findAllNeghID(blockID,surrBlksSummary)
       use gr_interfaceTypeDecl, ONLY: AllBlockRegions_t
       integer, intent(IN) :: blockID
       type (AllBlockRegions_t), intent(OUT) :: surrBlksSummary
     end subroutine gr_findAllNeghID
  end interface


  interface
     subroutine gr_checkGridState()
       implicit none
     end subroutine gr_checkGridState
  end interface

  interface
     subroutine gr_updateRefinement( gridChanged)
       implicit none
       logical, intent(out),OPTIONAL :: gridChanged
     end subroutine gr_updateRefinement
  end interface
  
  interface 
     subroutine gr_getInteriorBlkPtr(blockID,dataPtr,gridDataStruct)
       implicit none
       integer, intent(IN) :: blockID
       real,dimension(:,:,:,:),pointer :: dataPtr
       integer, intent(IN) :: gridDataStruct
     end subroutine gr_getInteriorBlkPtr
  end interface

  interface 
     subroutine gr_releaseInteriorBlkPtr(blockID,dataPtr,gridDataStruct)
       implicit none
       integer, intent(IN) :: blockID
       real,dimension(:,:,:,:),pointer :: dataPtr
       integer, intent(IN) :: gridDataStruct
     end subroutine gr_releaseInteriorBlkPtr
  end interface

  interface
     subroutine gr_GCAllocScratch(gridDataStruct,blkCnt,blkList,&
          indCnt,indList, gcCnt)
       integer, intent(IN) :: gridDataStruct, blkCnt, indCnt
       integer, dimension(blkCnt), intent(IN) :: blkList
       integer, dimension(indCnt), intent(IN) :: indList
       integer, dimension(NDIM), intent(IN) :: gcCnt
     end subroutine gr_GCAllocScratch
  end interface

  interface
     subroutine gr_GCReleaseScratch(gridDataStruct)
       integer, intent(IN) :: gridDataStruct
     end subroutine gr_GCReleaseScratch
  end interface

  interface
     subroutine gr_GCTransferOneBlk(mode,indCnt,indList,offset,&
          blkLimits,blkLimitsGC,&
          flatArray,blkArray)
       integer, intent(IN) :: indCnt,offset
       logical, intent(IN) :: mode
       integer,dimension(indCnt),intent(IN) :: indList
       integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits,blkLimitsGC
       real, pointer, dimension(:) :: flatArray
       real, pointer, dimension(:,:,:,:) :: blkArray
     end subroutine gr_GCTransferOneBlk
  end interface

  interface
     subroutine gr_setBlockType(blockID,type)
       integer, intent(IN) :: blockID, type
     end subroutine gr_setBlockType
  end interface


  interface
     subroutine gr_hypreGridStatus (blockCount, blockList, nvars)
       implicit none  
       integer, intent(IN):: blockCount
       integer, dimension(blockCount),intent(IN):: blockList
       integer, intent(IN),OPTIONAL :: nvars
     end subroutine gr_hypreGridStatus
  end interface
  
  interface
     subroutine gr_hypreSetupGrid (blockCount, blockList, nvars)
       implicit none 
       integer,                      intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       integer,OPTIONAL,             intent(IN) :: nvars
     end subroutine gr_hypreSetupGrid
  end interface
  
  interface
     subroutine gr_hypreGetFaceB (direction, iFactorB, blkLimits, blkLimitsGC, solnVec, flux, numVars)
       implicit none
       integer, intent(IN) :: direction
       integer, intent(IN) :: iFactorB
       integer, intent(IN) :: blkLimits (2,MDIM) 
       integer, intent(IN) :: blkLimitsGC (2,MDIM)
       real, intent(IN)    :: solnVec(NUNK_VARS, blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(HIGH,KAXIS))   
       real, intent(INOUT) :: flux(NFLUXES,blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(HIGH,KAXIS))  
       integer, intent(IN) :: numVars
     end subroutine gr_hypreGetFaceB
  end interface

  interface
     subroutine gr_hypreComputeB (blockCount, blockList, iVar, iFactorA, &
          iFactorB, dt, theta, bcTypes, bcValues, iFactorD)
       integer, intent(IN) :: iVar
       integer, intent(IN) :: iFactorB
       integer, intent(IN) :: iFactorA
       integer, OPTIONAL, intent(IN) :: iFactorD
       real, intent(IN) :: dt
       real, intent(IN) :: theta
       integer, intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       integer, intent(IN) :: bcTypes(6)
       real,    intent(IN) :: bcValues(2,6)  
     end subroutine gr_hypreComputeB
  end interface
  
  interface 
     subroutine gr_hypreCreateMatrix(iVar, iFactorB, iFactorA, bcTypes, bcValues, &
          dt, alpha, blockCount, blockList, JacobiMatrix)
       integer, intent(IN) :: iVar
       integer, intent(IN) :: iFactorB
       integer, intent(IN) :: iFactorA
       integer, intent(IN) :: bcTypes(6)
       real,    intent(IN) :: bcValues(2,6)
       real,    intent(IN) :: dt
       real,    intent(IN) :: alpha
       integer, intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       logical, intent(IN) :: JacobiMatrix
     end subroutine gr_hypreCreateMatrix
     subroutine gr_hypreCreateMatrixFcB(iVar, iFactorB, iFactorA, bcTypes, bcValues, &
          dt, alpha, blockCount, blockList, JacobiMatrix)
       integer, intent(IN) :: iVar
       integer, intent(IN) :: iFactorB
       integer, intent(IN) :: iFactorA
       integer, intent(IN) :: bcTypes(6)
       real,    intent(IN) :: bcValues(2,6)
       real,    intent(IN) :: dt
       real,    intent(IN) :: alpha
       integer, intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       logical, intent(IN) :: JacobiMatrix
     end subroutine gr_hypreCreateMatrixFcB
     subroutine gr_hypreCreateMatrixAnisoCond(iVar, iFactorB, iFactorA, bcTypes, bcValues, &
          dt, alpha, blockCount, blockList, JacobiMatrix)
       integer, intent(IN) :: iVar
       integer, intent(IN) :: iFactorB
       integer, intent(IN) :: iFactorA
       integer, intent(IN) :: bcTypes(6)
       real,    intent(IN) :: bcValues(2,6)
       real,    intent(IN) :: dt
       real,    intent(IN) :: alpha
       integer, intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       logical, intent(IN) :: JacobiMatrix
     end subroutine gr_hypreCreateMatrixAnisoCond
     subroutine gr_hypreCreateMatrixFcBAnisoCond(iVar, iFactorB, iFactorA, bcTypes, bcValues, &
          dt, alpha, blockCount, blockList, JacobiMatrix)
       integer, intent(IN) :: iVar
       integer, intent(IN) :: iFactorB
       integer, intent(IN) :: iFactorA
       integer, intent(IN) :: bcTypes(6)
       real,    intent(IN) :: bcValues(2,6)
       real,    intent(IN) :: dt
       real,    intent(IN) :: alpha
       integer, intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       logical, intent(IN) :: JacobiMatrix
     end subroutine gr_hypreCreateMatrixFcBAnisoCond
  end interface


  interface
     subroutine gr_hypreApplyBcToFaceAnisoCond(blkLimits,blkLimitsGC,part,var,iFactorB,bcType,direction, &
          bcValue, dt, theta, del, Lower, scalefactor, faceArea, solnVec, blockID)
       integer, intent(IN) :: blkLimits (2,MDIM) 
       integer, intent(IN) :: blkLimitsGC (2,MDIM)
       integer, intent(IN) :: part
       integer, intent(IN) :: var
       integer, intent(IN) :: iFactorB
       integer, intent(IN) :: bcType
       integer, intent(IN) :: direction
       real,    intent(IN) :: bcValue(2)
       real,    intent(IN) :: dt
       real,    intent(IN) :: theta
       real,    intent(IN) :: del
       integer, intent(IN) :: Lower(MDIM), blockID 
       real,    intent(IN) :: scalefactor
       real,    intent(IN) :: faceArea(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
       real,    intent(IN) :: solnVec(NUNK_VARS, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))   
     end subroutine gr_hypreApplyBcToFaceAnisoCond

  end interface


  
  interface 
     subroutine gr_hypreCreateMatrix1Blk(iFactorB, dt, &
          theta, &
          blkLimits, datasize, solnVec, nentries, faces, mylevel, var, faceAreas, &
          del, negh_del, xflux,yflux,zflux, &
          lb, blockID)
       implicit none
       integer, intent(IN) :: iFactorB
       real,    intent(IN) :: dt
       real,    intent(IN) :: theta
       real, intent(IN), dimension(MDIM)     :: del
       real, intent(IN), dimension(2*MDIM, MDIM) :: negh_del
       integer, intent(IN), dimension(2,MDIM):: blkLimits 
       integer, intent(IN) :: datasize(MDIM)
       integer, intent(IN) ::  mylevel
       integer, intent(IN) ::  var
       integer, intent(IN) ::  nentries
       integer, dimension(2,MDIM), intent(IN) :: faces 
       real, intent(IN), dimension(:,:,:,:) :: xflux, yflux, zflux
       real, intent(IN) :: faceAreas  (:,:,:,:)  
       integer, intent(IN) :: lb, blockID
       real, POINTER, DIMENSION(:,:,:,:) :: solnVec
     end subroutine gr_hypreCreateMatrix1Blk
  end interface

  interface
     subroutine gr_xyzToBlockLevel(lev, xyz, ijk)
       integer, intent(in) :: lev
       real, intent(in) :: xyz(NDIM)
       integer, intent(out) :: ijk(NDIM)
     end subroutine gr_xyzToBlockLevel
  end interface

  interface
     Subroutine gr_xyzToBlock(xyz, procID, blkID)
       real, dimension(MDIM),intent(IN) :: xyz
       integer, intent(OUT) :: procID
       integer, intent(OUT) :: blkID
     End Subroutine gr_xyzToBlock
  end interface

end module gr_interface
