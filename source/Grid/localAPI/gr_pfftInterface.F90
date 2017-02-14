!!****ih* source/Grid/localAPI/gr_pfftInterface
!!
!! NAME 
!!   gr_pfftInterface
!!
!! SYNOPSIS
!!   
!!   use gr_pfftInterface
!!
!! DESCRIPTION
!! 
!!  This is the header file for the Pfft solver that defines its
!!  interfaces.
!!
!!***

module gr_pfftinterface
#include "constants.h"
  implicit none

  interface
     subroutine gr_pfftInit()
     end subroutine gr_pfftInit
  end interface

  interface
     subroutine gr_pfftFinalize()
     end subroutine gr_pfftFinalize
  end interface

  interface
     subroutine gr_pfftDcftForward(inArray,outArray,trig,len,lda,numVec,&
          transformType,scale)
       integer, intent(IN) :: len, lda, numVec, transformType
       real, intent(IN) :: scale
       real, dimension(:),intent(IN) :: trig,inArray
       real,dimension(:), intent(OUT) :: outArray
     end subroutine gr_pfftDcftForward
  end interface
  
  interface
     subroutine gr_pfftDcftInverse(inArray,outArray,trig,len,lda,numVec,&
          transformType,scale)
       integer, intent(IN) :: len, lda, numVec, transformType
       real, intent(IN) :: scale
       real, dimension(:),intent(IN) :: trig,inArray
       real,dimension(:), intent(OUT) :: outArray
     end subroutine gr_pfftDcftInverse
  end interface

  interface
     subroutine gr_pfftDisperse(inArray,outArray,nx,ny,nx1)
       integer,intent(IN) :: nx,ny,nx1
       real,intent(IN),dimension(:) :: inArray
       real,intent(OUT), dimension(:) :: outArray
     end subroutine gr_pfftDisperse
  end interface

  interface
     subroutine gr_pfftGenMap()
     end subroutine gr_pfftGenMap
  end interface

  interface
     subroutine gr_pfftInitMetaData(ndim)
       integer, intent(IN) :: ndim
     end subroutine gr_pfftInitMetaData
  end interface
  
  interface
     subroutine gr_pfftIntersperse(inArray,outArray,nx,ny,nx1)
       integer, intent(IN) :: nx,ny,nx1
       real,dimension(:), intent(IN) :: inArray
       real,dimension(:), intent(OUT) :: outArray
     end subroutine gr_pfftIntersperse
  end interface

  interface
     subroutine gr_pfftLocalTranspose(inArray,outArray,nx,nx1,baseDatType)
       integer,intent(IN) :: nx,nx1,baseDatType
       real,intent(IN),dimension(:) :: inArray
       real,intent(OUT), dimension(:) :: outArray
     end subroutine gr_pfftLocalTranspose
  end interface

  interface
     subroutine gr_pfftSetupDim(length, transformType, trig, scale, inDatSize, outDatSize, datInc, factor)
       integer,intent(IN) :: length,transformType
       real,dimension(:), pointer :: trig
       real,intent(OUT) :: scale
       integer,intent(OUT) :: inDatSize, outDatSize, datInc, factor
     end subroutine gr_pfftSetupDim
  end interface
  
  interface gr_pfftTranspose
     subroutine gr_pfftTranspose(dir,baseDatType,inArray,outArray,inlen,outlen,pr,comm)
       integer, intent(IN) :: pr,comm,dir,baseDatType
       integer, dimension(MDIM), intent(IN) :: inlen,outlen
       real,dimension(:),intent(INOUT) :: inArray,outArray
     end subroutine gr_pfftTranspose
     subroutine gr_pfftTranspose3DArr(dir,baseDatType,inArray,outArray,inlen,outlen,pr,comm)
       integer, intent(IN) :: pr,comm,dir,baseDatType
       integer, dimension(MDIM), intent(IN) :: inlen,outlen
       real,dimension(:,:,:),intent(INOUT) :: inArray
       real,dimension(:,:,:),intent(INOUT) :: outArray
     end subroutine gr_pfftTranspose3DArr
  end interface

  interface
     subroutine gr_pfftGetLocalLimits(axis1,axis2)
       integer,intent(IN) :: axis1, axis2
     end subroutine gr_pfftGetLocalLimits
  end interface


  interface
     subroutine gr_pfftWave()
     end subroutine gr_pfftWave
  end interface

  interface
     subroutine gr_pfftDerivs(transArray)
       real,dimension(:),intent(INOUT) :: transArray
     end subroutine gr_pfftDerivs
  end interface
  
  interface
     subroutine gr_pfftGetProcGrid(dims, pfftMyPE, pfftNumProcs, &
          pfftGlobalLen, pfftProcGrid)
       implicit none
       integer, intent(IN) :: dims, pfftMyPE, pfftNumProcs
       integer, dimension(MDIM), intent(IN) :: pfftGlobalLen
       integer, dimension(MDIM), intent(OUT) :: pfftProcGrid
     end subroutine gr_pfftGetProcGrid
  end interface


  interface
     subroutine gr_pfftGenMapHelper(axis, fragmentPtr, maxSingleProcData)
       integer, intent(IN) :: axis
       integer, dimension(:), intent(IN) :: fragmentPtr
       integer, intent(OUT) :: maxSingleProcData
     end subroutine gr_pfftGenMapHelper
  end interface

  !Unable to specify intent for "buffer" and "pfftArray".  This is because
  !sometimes the procedure is called with "pfftArray" which is an intent(IN)
  !in the calling subroutine, and sometimes "pfftArray" is an intent(OUT) in the
  !calling subroutine.  Specifying intent(INOUT) does not work.
  interface
     subroutine gr_pfftBufferTransfer(direction, axis, buffer, pfftArray)
       integer, intent(IN) :: direction, axis
       real, dimension(:,:) :: buffer
       real, dimension(:) :: pfftArray       
     end subroutine gr_pfftBufferTransfer
  end interface

  interface
     subroutine gr_pfftGetDestPfftCoords(startPos, pfftProcCoords)       
       implicit none
       integer, dimension(1:MDIM), intent(IN) :: startPos
       integer, dimension(1:MDIM), intent(OUT) :: pfftProcCoords
     end subroutine gr_pfftGetDestPfftCoords
  end interface

  interface
     subroutine gr_pfftCopyToSendMap(blockID, blockStartPos, blockEndPos, bsize, axis)
       implicit none
       integer, intent(IN) :: blockID
       integer, dimension(1:MDIM), intent(IN) :: blockStartPos, blockEndPos
       integer, intent(IN) :: bsize, axis
     end subroutine gr_pfftCopyToSendMap
  end interface

  interface
     subroutine gr_pfftHandleJaxisFragments(direction, gridVar)
       implicit none
       integer, intent(IN) :: direction, gridVar
     end subroutine gr_pfftHandleJaxisFragments
  end interface

  interface
     subroutine gr_pfftHandleKaxisFragments(direction)
       implicit none
       integer, intent(IN) :: direction
     end subroutine gr_pfftHandleKaxisFragments
  end interface

  interface
     subroutine gr_pfftPrintCommBuffers(buffer, map, logUnit)       
       implicit none
       real, dimension(:,:), intent(IN) :: buffer
       integer, dimension(:,:), intent(IN) :: map
       integer, intent(IN) :: logUnit
     end subroutine gr_pfftPrintCommBuffers
  end interface

  interface
     subroutine gr_pfftGroupUsableProcesses(myPE, globalProcs, originalComm, newComm)
       implicit none
       integer, intent(IN) :: myPE, globalProcs, originalComm
       integer, intent(OUT) :: newComm
     end subroutine gr_pfftGroupUsableProcesses
  end interface

 interface
     subroutine gr_pfftGetLocalLimitsAnytime(axis1,axis2,meaxis,&
          currentGridShape,baseDatType,currentLocalLimits)
       implicit none
       integer,intent(IN), target :: axis1
       integer,intent(IN) :: axis2
       integer,intent(IN), OPTIONAL, target ::  meaxis
       integer, dimension(MDIM), intent(IN), OPTIONAL, target :: currentGridShape
       integer, intent(IN), OPTIONAL, target :: baseDatType
       integer, intent(INOUT), OPTIONAL, target :: currentLocalLimits(LOW:HIGH,MDIM)
     end subroutine gr_pfftGetLocalLimitsAnytime
  end interface

  interface
     subroutine gr_pfftPoissonDirect (iDirection, solveflag, inSize, &
          localSize, globalSize, transformType, inArray, outArray)
       integer, intent(in)    :: iDirection, solveflag, inSize  
       integer, dimension(MDIM),intent(in) :: localSize,globalSize,transformType
       real,  dimension(inSize),intent(in) :: inArray
       real,dimension(inSize), intent(out) :: outArray
     end subroutine gr_pfftPoissonDirect
  end interface

  interface
     subroutine gr_pfftMapToInput(gridVar, pfftInputArray)
       implicit none
       integer, intent(IN) :: gridVar
       real, dimension(:), target :: pfftInputArray
     end subroutine gr_pfftMapToInput
  end interface

  interface
     subroutine gr_pfftMapFromOutput(gridVar, pfftOutputArray)
       implicit none
       integer, intent(IN) :: gridVar
       real, dimension(:), target :: pfftOutputArray
     end subroutine gr_pfftMapFromOutput
  end interface

  interface
     subroutine gr_pfftInitialiseStorage()
       implicit none
     end subroutine gr_pfftInitialiseStorage
  end interface

  interface
     subroutine gr_pfftFinaliseStorage()
       implicit none
     end subroutine gr_pfftFinaliseStorage
  end interface

  interface
     subroutine gr_pfftCreateSendNode(flashProcID, flashBlockID, flashStartPos, &
          flashEndPos, pfftProcID, pfftStartPos, pfftEndPos)
       implicit none
       integer, intent(IN) :: flashProcID, flashBlockID
       integer, dimension(1:MDIM), intent(IN) :: flashStartPos, flashEndPos
       integer, intent(IN) :: pfftProcID
       integer, dimension(1:MDIM), intent(IN) :: pfftStartPos, pfftEndPos
     end subroutine gr_pfftCreateSendNode
  end interface

  interface
     subroutine gr_pfftCommunicateNodeMetaData()
       implicit none
     end subroutine gr_pfftCommunicateNodeMetaData
  end interface

  interface
     subroutine gr_pfftGridPointTable(pfft_inLen)
       integer, dimension(1:MDIM), intent(IN) :: pfft_inLen
     end subroutine gr_pfftGridPointTable
  end interface


  interface
     function gr_pfftMakePencilIn3dSpace &
          (pencilGlobalLen, totalProcs, gr_pfftFnArgConstraint)
       use gr_pfftInterfaceTypeDecl, ONLY : PossibleGrid_t
       implicit none
       integer, dimension(MDIM), intent(IN) :: pencilGlobalLen
       integer, intent(IN) :: totalProcs
       interface     
          logical function gr_pfftFnArgConstraint &
               (pencilGlobalLen, totalProcs, iProcs, jProcs, kProcs)
            implicit none
            integer, dimension(MDIM), intent(IN) :: pencilGlobalLen
            integer, intent(IN) :: totalProcs, iProcs, jProcs, kProcs
          end function gr_pfftFnArgConstraint
       end interface
       type (PossibleGrid_t) :: gr_pfftMakePencilIn3dSpace  !Return.
     end function gr_pfftMakePencilIn3dSpace
  end interface


  interface
     logical function gr_pfftFnArgHardConstraint &
          (pencilGlobalLen, totalProcs, iProcs, jProcs, kProcs)
       implicit none
       integer, dimension(MDIM), intent(IN) :: pencilGlobalLen
       integer, intent(IN) :: totalProcs, iProcs, jProcs, kProcs
     end function gr_pfftFnArgHardConstraint
  end interface


  interface
     logical function gr_pfftFnArgMediumConstraint &
          (pencilGlobalLen, totalProcs, iProcs, jProcs, kProcs)
       implicit none
       integer, dimension(MDIM), intent(IN) :: pencilGlobalLen
       integer, intent(IN) :: totalProcs, iProcs, jProcs, kProcs
     end function gr_pfftFnArgMediumConstraint
  end interface


  interface
     logical function gr_pfftFnArgEasyConstraint &
          (pencilGlobalLen, totalProcs, iProcs, jProcs, kProcs)
       implicit none
       integer, dimension(MDIM), intent(IN) :: pencilGlobalLen
       integer, intent(IN) :: totalProcs, iProcs, jProcs, kProcs
     end function gr_pfftFnArgEasyConstraint
  end interface

  interface
     subroutine gr_pfftSpecifyTransform(transformType,baseDatType, bcTypes)
       implicit none
       integer, dimension(MDIM), intent(OUT) :: transformType
       integer, dimension(0:MDIM), intent(OUT), OPTIONAL :: baseDatType
       integer, dimension(2*MDIM), intent(IN),  OPTIONAL :: bcTypes
     end subroutine gr_pfftSpecifyTransform
  end interface

  interface
     subroutine gr_pfftGenSingleMap(solveLevel, leafMapMode)
       implicit none
       integer, intent(IN) :: solveLevel
       logical, optional, intent(IN) :: leafMapMode
     end subroutine gr_pfftGenSingleMap
  end interface

  interface
     subroutine gr_pfftValidateSelectedLevel(inOutLevel)
       implicit none
       integer, intent(INOUT) :: inOutLevel
     end subroutine gr_pfftValidateSelectedLevel
  end interface

  interface
     subroutine gr_pfftCreateSendFragment(flashProcID, flashBlockID, &
          lBlockFragStart, lBlockFragEnd, lActualBlockFragStart, lActualBlockFragEnd, &
          blkType, blkRefLev, solveLevel, pfftProc, &
          lPencilFragStart, lPencilFragEnd)
       implicit none
       integer, intent(IN) :: flashProcID, flashBlockID
       integer, dimension(1:MDIM), intent(IN) :: lBlockFragStart, lBlockFragEnd, &
            lActualBlockFragStart, lActualBlockFragEnd
       integer, intent(IN) :: blkType, blkRefLev, solveLevel, pfftProc
       integer, dimension(1:MDIM), intent(IN) :: lPencilFragStart, lPencilFragEnd 
     end subroutine gr_pfftCreateSendFragment
  end interface

end module gr_pfftinterface
