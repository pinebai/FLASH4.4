!!****ih* source/Grid/localAPI/gr_ptInterface
!!
!! NAME 
!!
!!   gr_ptInterface
!!
!! SYNOPSIS
!!   
!!   use gr_ptInterface
!!
!! DESCRIPTION
!! 
!! This is the header file for the GridParticles
!! subunit that defines its interfaces.
!!
!!***

module gr_ptInterface
  implicit none
#include "constants.h"
#include "Flash.h"

  interface 
     subroutine gr_ptInit()    
     end subroutine gr_ptInit
  end interface

  interface 
     subroutine gr_ptFinalize()
     end subroutine gr_ptFinalize
  end interface

  interface 
     subroutine gr_ptMarkRefineDerefine ()
     end subroutine gr_ptMarkRefineDerefine
  end interface

  interface
     subroutine gr_ptLocalMatch(particles,localNumParticles,&
          part_props,maxParticlesPerProc,&
          sourceBuf,numSource,destBuf,numDest)
       integer, intent(IN) :: part_props
       integer, intent(IN) :: numSource
       integer, intent(OUT) :: numDest
       integer, intent(INOUT) :: localNumParticles
       integer, intent(IN) :: maxParticlesPerProc
       real,dimension(part_props,maxParticlesPerProc),intent(INOUT) :: particles
       real,dimension(part_props,numSource),intent(INOUT) :: sourceBuf,destBuf
     end subroutine gr_ptLocalMatch
  end interface

  interface

     subroutine gr_ptMoveSieve(particles,localNumParticles,part_props,&
          maxParticlesPerProc, numDest, matchProp)
 
       integer,intent(INOUT) :: localNumParticles
       integer,intent(IN) :: part_props
       integer,intent(IN) :: maxParticlesPerProc
       integer,intent(INOUT) :: numDest
       integer,optional, intent(IN) :: matchProp
       
       real, dimension(part_props, maxParticlesPerProc),intent(INOUT) :: particles
     end subroutine gr_ptMoveSieve
  end interface

  interface 
     subroutine gr_ptMoveOffBlk(particles,part_props,localNumParticles,bufferDim2,destBuf,numDest)
       integer, intent(IN) :: part_props,bufferDim2
       integer,intent(INOUT)::localNumParticles,numDest
       real,dimension(part_props,bufferDim2),intent(INOUT)::particles,destBuf
     end subroutine gr_ptMoveOffBlk
  end interface

  interface
     subroutine gr_ptMoveOffProc(face,axis,index, propCount, maxPerProc,boundary,&
          lnegh,rnegh,corner,localNum,particles)
       integer, intent(IN) :: face, axis, index, propCount, maxPerProc
       logical,intent(IN) :: boundary
       integer, intent(IN) :: lnegh,rnegh
       real,dimension(LOW:HIGH),intent(IN) :: corner
       integer,intent(INOUT) :: localNum
       real,dimension(propCount,maxPerProc),intent(INOUT) :: particles
     end subroutine gr_ptMoveOffProc
  end interface

  interface
     subroutine gr_ptOneFaceBC(particle,propCount,axis,face,blockID,lostParticles)
!! ,moved)
       integer, intent(IN) :: propCount
       real,dimension(propCount),intent(INOUT)::particle
       integer, intent(IN) :: axis, face, blockID
       integer, intent(INOUT) :: lostParticles
 !!      logical, intent(INOUT) :: moved
     end subroutine gr_ptOneFaceBC
  end interface

  interface 
     subroutine gr_ptHandleExcess(particles,propCount,localNum,maxPerProc)
       integer,intent(IN) :: propCount
       integer,intent(IN) :: maxPerProc
       integer,intent(INOUT) :: localNum
       real,dimension(propCount,maxPerProc),intent(INOUT) :: particles
     end subroutine gr_ptHandleExcess
  end interface

  interface 
     subroutine gr_ptApplyBC(particle, propCount, blockID, lostParticles, moved, negh, bndBox)

       integer, intent(IN)    :: blockID,propCount
       real,dimension(propCount),intent(INOUT)::particle
       integer, intent(INOUT) :: lostParticles
       logical, intent(INOUT) :: moved
       integer, dimension(MDIM),intent(INOUT) :: negh
       real, dimension(LOW:HIGH, MDIM), intent(IN) :: bndBox
     end subroutine gr_ptApplyBC
  end interface

  interface
     subroutine gr_ptMapInit()
     end subroutine gr_ptMapInit
  end interface

  interface
     subroutine gr_ptApplyBCsOneBlk (blkLimits, blkLimitsGC, blockID)
       
       implicit none
       integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: blkLimits, blkLimitsGC
       integer, intent(IN) :: blockID
     end subroutine gr_ptApplyBCsOneBlk
  end interface

  interface
     subroutine gr_ptStoreOffBlockCells(particlesPerBlk, blockList, blockCount, blkLimitsGC, blkSize, guard, BufferSize)
       
       integer,dimension(MAXBLOCKS), intent(IN) :: particlesPerBlk, blockList
       integer,intent(IN) :: blockCount
       integer,dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
       integer,dimension(MDIM), intent(IN) :: blkSize, guard
       integer, intent(OUT) :: BufferSize

     end subroutine gr_ptStoreOffBlockCells
  end interface

  interface 
     subroutine gr_ptFindNegh(srcBlkID,guardCellID,negh,neghCornerID,numNegh)

       integer, intent(IN) :: srcBlkID
       integer, dimension(MDIM), intent(IN) :: guardCellID
       integer, dimension(BLKNO:TYPENO,ABSMAXNEGH),intent(OUT):: negh
       integer, dimension(MDIM,ABSMAXNEGH),intent(OUT) :: neghCornerID
       integer, intent(OUT) :: numNegh
       
     end subroutine gr_ptFindNegh
  end interface


  interface 
     subroutine gr_ptSearchBlk (cornerID,negh)

       integer,dimension(MDIM), intent(IN) :: cornerID
       integer,dimension(BLKNO:TYPENO), intent(INOUT) :: negh

     end subroutine gr_ptSearchBlk
  end interface


  interface
     subroutine gr_ptGetSrcDestCoords(blkSize, guard, guardCellID, srcCornerID, srcStride, destCornerID, &
          guardCoords, negh, srcCoords, destCoords)  

       integer,dimension(MDIM), intent(IN)  :: blkSize, guard, guardCellID, srcCornerID, srcStride, destCornerID
       integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: guardCoords
       integer,dimension(BLKNO:TYPENO), intent(IN):: negh
       integer,dimension(LOW:HIGH,MDIM), intent(OUT)  :: srcCoords, destCoords

     end subroutine gr_ptGetSrcDestCoords
  end interface


  interface 
     subroutine gr_ptSameProcMap (srcCoords, destCoords, negh, varGrid)

       integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: srcCoords, destCoords
       integer,dimension(BLKNO:TYPENO), intent(IN):: negh
       integer,intent(IN) :: varGrid

     end subroutine gr_ptSameProcMap
  end interface


  interface 
     subroutine gr_ptOffProcMap(srcCoords, destCoords, bufferSize, sendBuf, &
          sendSize, sendBufPtr, negh, neghCornerID)

       integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: srcCoords, destCoords
       integer,intent(IN) :: bufferSize
       real,dimension(bufferSize),intent(INOUT) :: sendBuf
       integer,intent(INOUT) :: sendSize
       integer,intent(INOUT) :: sendBufPtr
       integer,dimension(BLKNO:TYPENO), intent(IN):: negh
       integer,dimension(MDIM), intent(IN) :: neghCornerID

     end subroutine gr_ptOffProcMap
  end interface


  interface 
     subroutine gr_ptProlongSmear(ioff,joff,koff,prolongedSection)
       
       integer, INTENT(in) :: ioff, joff, koff
       real,  dimension(1:2, 1:2, 1:2), INTENT(out) :: prolongedSection

     end subroutine gr_ptProlongSmear
  end interface


  interface 
     subroutine gr_ptMoveMappedData(varGrid,bufferSize,sendBuf,sendCount,recvBuf)

       integer,intent(IN) :: varGrid
       integer,intent(IN) :: bufferSize
       real,dimension(bufferSize),intent(INOUT) :: sendBuf
       integer,intent(INOUT) :: sendCount
       real,dimension(bufferSize),intent(INOUT) :: recvBuf

     end subroutine gr_ptMoveMappedData
  end interface


  interface
     subroutine gr_ptPackUnpackData(varGrid, bufferSize, sendBuf, sendCount, recvBuf, recvCount)

       integer,intent(IN) :: varGrid
       integer,intent(IN) :: bufferSize
       real,dimension(bufferSize),intent(OUT) :: sendBuf
       integer,intent(OUT) :: sendCount
       real,dimension(bufferSize),intent(IN) :: recvBuf
       integer,intent(IN) :: recvCount
       
     end subroutine gr_ptPackUnpackData
  end interface


  interface
     subroutine gr_ptParticleAtFcBdry(partID, fcBdry)

       integer, intent(IN) :: partID
       logical, intent(OUT) :: fcBdry
     
     end subroutine gr_ptParticleAtFcBdry
  end interface


  interface
     subroutine gr_ptExchangePartialMap(blkLimits,blkLimitsGC,bufSize,axis,face,&
          sendBuf,recvBuf)

       implicit none
       integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits,blkLimitsGC
       integer,intent(IN) :: bufSize,axis,face
       real,intent(inout),dimension(bufSize):: sendBuf,recvBuf

     end subroutine gr_ptExchangePartialMap
  end interface


  interface
     subroutine gr_ptVerifyBlock(particles,propCount,localNumParticles,maxParticlesPerProc)

       implicit none
       integer,intent(IN) :: propCount
       integer,intent(IN) :: localNumParticles
       integer,intent(IN) :: maxParticlesPerProc
       real, dimension(propCount, maxParticlesPerProc),intent(IN) :: particles

     end subroutine gr_ptVerifyBlock
  end interface


  interface
     subroutine gr_ptParseMetadata(bufferSize, dataBuffer, headerPtr, &
          negh, neghCornerID, sectionCoords, numbElements)

       implicit none
       integer, intent(IN) :: bufferSize
       real, dimension(bufferSize), intent(IN) :: dataBuffer
       integer, intent(IN) :: headerPtr
       integer, dimension(MDIM), intent(OUT) :: negh, neghCornerID
       integer, dimension(LOW:HIGH,MDIM), intent(OUT) :: sectionCoords
       integer, intent(OUT) :: numbElements

     end subroutine gr_ptParseMetadata
  end interface


  interface
     subroutine gr_ptDumpState(bufferSize, dataBuffer, bufferContentSize)

       implicit none
       integer, intent(IN) :: bufferSize
       real, dimension(bufferSize),intent(IN) :: dataBuffer
       integer,intent(IN) :: bufferContentSize

     end subroutine gr_ptDumpState
  end interface


  interface
     subroutine gr_ptGetChildData(guardCellID, blkSize, parentCornerID, &
          srcStride, neghBlkID, neghProcID, negh, neghCornerID, numNegh)
       implicit none
       integer, dimension(MDIM), intent(IN) :: guardCellID, blkSize, &
            parentCornerID, srcStride
       integer, intent(IN) :: neghBlkID, neghProcID
       integer, dimension(BLKNO:TYPENO,ABSMAXNEGH), intent(OUT) :: negh
       integer, dimension(MDIM,ABSMAXNEGH), intent(OUT) :: neghCornerID
       integer, intent(OUT) :: numNegh
     end subroutine gr_ptGetChildData
  end interface

  interface
     subroutine gr_ptMove(databuf,propCount, localCount,maxCount, moveDone)

       integer,intent(INOUT) :: localCount
       integer,intent(IN) :: maxCount, propCount
       
       real, dimension(propCount, maxCount),intent(INOUT) :: databuf
       logical, intent(INOUT) :: moveDone
       
     end subroutine gr_ptMove
  end interface
  
  interface
     subroutine gr_ptMovePttoPt(dataBuf,propCount,maxCount,localCount, numDest)
       integer,intent(INOUT) :: localCount
       integer,intent(IN) :: propCount, maxCount, numDest
       real, dimension(propCount, maxCount),intent(INOUT) :: dataBuf
     end subroutine gr_ptMovePttoPt
  end interface


  interface

     subroutine gr_ptSetIndices(index_list,count)
       
       integer, intent(IN) :: count
       integer,dimension(count), intent(IN) :: index_list

     end subroutine gr_ptSetIndices
  end interface

  interface

     subroutine gr_ptResetIndices(index_list,count)
       
       integer, intent(IN) :: count
       integer,dimension(count), intent(IN) :: index_list

     end subroutine gr_ptResetIndices
  end interface

end module gr_ptInterface
