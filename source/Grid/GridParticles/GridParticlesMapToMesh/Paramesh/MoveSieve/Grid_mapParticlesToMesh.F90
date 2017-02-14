!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/MoveSieve/Grid_mapParticlesToMesh
!!
!! NAME
!!  Grid_mapParticlesToMesh
!!
!! SYNOPSIS
!!
!!  call Grid_mapParticlesToMesh(real(INOUT) :: particles(part_props,numParticles),
!!                               integer(IN) :: part_props,
!!                               integer(IN) :: numParticles,
!!                               integer(IN) :: maxParticlesPerProc,
!!                               integer(IN) :: propPart,
!!                               integer(IN) :: varGrid, 
!!                      OPTIONAL,integer(IN) :: mode,
!!                      OPTIONAL,integer(IN) :: ptInfo)
!!
!! DESCRIPTION
!!
!! Routine to map a quantity defined for particles onto the  mesh.
!! The approximate steps are:
!! 1.  Determine the size of the send / receive buffer.  This will eventually 
!!     contain the guard cells that need to be delivered to an off-processor block.
!!     (Whilst determining this size, store valuable information in a data-structure
!!     named "gr_ptDomain".)
!! 2.  Loop over each block on this processor.  If it contains at least one particle 
!!     then smear these particles over a temporary block.  The details of the smear 
!!     depend on the particle mapping scheme, e.g. NGP, CIC, TSC.  The temporary block 
!!     is contiguous in memory to improve performance.
!! 3.  Apply external boundary conditions to the temporary block.  
!! 4.  Copy internal section of temporary block into the source block.  
!! 5.  The boundaries of the temporary block belong to another block on the source 
!!     processor or a block on another processor.  Perform any prolongation / 
!!     restriction for when the touching blocks are at different resolutions.  
!!     Copy the guard cells directly to another block on this processor, or 
!!     place the guard cells in the send / receive buffer.
!! 6.  Send the send buffer to another processor, and receive from another processor.
!!
!! ARGUMENTS
!!  particles : List of particles. It is two dimensional real array, the first dimension
!!              represents each particle's properties, and second dimension is index to
!!              particles.
!!
!!  part_props : number of properties in the particles datastructure
!!
!!  numParticles : While coming in it contains the current number of particles mapped to
!!                      this processor. After all the data structure movement, the number
!!                      of local particles might change, and the new value is put back into it
!!  maxParticlesPerProc : The size of the buffer allocated for particles at compile time. This
!!                        is the largest number of particles that a simulation can have on one
!!                        processor at any time.
!!  varGrid:   Index of gridded variable to receive interpolated
!!                           quantity
!!  propPart:   Index of particle attribute to interpolate onto the mesh
!!  mode:       (Optional) If zero (default), zero varGrid first; if
!!                              nonzero, do not zero varGrid first
!!  ptInfo:     (Optional) additional info to pass to the Particles unit
!!                              in some calls, for debugging
!!
!! NOTES
!!   This subroutine assumes that each particle is associated with the correct
!!   block.  To verify whether this is the case, the application should be
!!   resetup with -defines=DEBUG_GRIDPARTICLES and then rebuilt.  This ensures
!!   that the particle movement in Grid_moveParticles is correct.
!!***

!!REORDER(4): solnVec


subroutine Grid_mapParticlesToMesh (particles,part_props,numParticles,&
     maxParticlesPerProc,propPart, varGrid, mode, ptInfo)
  
  use gr_ptData, ONLY : gr_ptBlkList, gr_ptBlkCount, gr_ptBuf
  use gr_ptMapData, ONLY : gr_ptDomain, NUMBGUARDREGIONS, gr_ptSmearLen
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY: Logfile_stampMessage
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_meshMe
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getListOfBlocks, Grid_sortParticles, &
       Grid_getBlkCornerID
  use gr_ptInterface, ONLY :  gr_ptStoreOffBlockCells, gr_ptSameProcMap, &
       gr_ptOffProcMap, gr_ptMoveMappedData,gr_ptApplyBCsOneBlk
  use Particles_interface, ONLY : Particles_mapToMeshOneBlk 
  use tree, ONLY : lrefine

  implicit none

#include "constants.h"
#include "Flash.h"
#include "gr_ptMapToMesh.h"

  integer,intent(IN) :: maxParticlesPerProc,numParticles, part_props
  real,dimension(part_props,numParticles),intent(INOUT) :: particles
  integer, INTENT(in) :: propPart, varGrid
  integer, INTENT(in), optional :: mode
  integer, INTENT(in), optional :: ptInfo

  integer,parameter :: particleTypes=1
  integer,dimension(MAXBLOCKS,particleTypes) :: particlesPerBlk
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer :: blockID,pStart,pEnd, vmode
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC,srcCoords,destCoords
  integer,dimension(MDIM) :: blkSize,blkSizeGC,guard,srcCornerID,srcStride
  real, allocatable, dimension(:) :: sendBuf, recvBuf
  integer,dimension(BLKID:REFLEVELDIF):: negh
  integer,dimension(MDIM) :: neghCornerID
  integer :: sendBufPtr, blkNo, numNegh, n, sendCount
  integer :: BufferSize, regionIter, ib, jb, kb, ie, je, ke, error
  integer :: localNumParticles

  localNumParticles=numParticles

  if (present(mode)) then
     vmode = mode
  else
     vmode = 0
  end if

  !print *, "Processor", gr_meshMe, "in Grid_mapParticlesToMesh"

  !! We are doing this because in paramesh all blocks have the
  !! same dimension, and this is paramesh specific implementation
  !! so this calculation can be kept out of the loop.
  blockID=1
  call Grid_getBlkIndexLimits(blockID,blkLimits, blkLimitsGC)

  guard=blkLimits(LOW,:)-blkLimitsGC(LOW,:)
  blkSizeGC=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
  blkSize=blkLimits(HIGH,:)-blkLimits(LOW,:)+1


  !! All the blocks are zeroed initially, rather than one at 
  !! a time within the blkNo loop.  This is to prevent solnVec being
  !! zeroed when it contains valid values from an earlier iteration 
  !! of blkNo.
  call Grid_getListOfBlocks(LEAF,gr_ptBlkList,gr_ptBlkCount)

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
  print *, "Processor", gr_meshMe, "LEAF min refinement:", &
       minval(lrefine(gr_ptBlkList(1:gr_ptBlkCount))), "LEAF max refinement:", &
       maxval(lrefine(gr_ptBlkList(1:gr_ptBlkCount))), "no.particles:", numParticles
#endif

  do blkNo = 1, gr_ptblkCount, 1
     call Grid_getBlkPtr(gr_ptblkList(blkNo),solnVec,CENTER)
     if(vmode /= 0) then
        solnVec(varGrid,blkLimitsGC(LOW,IAXIS):blkLimitsGC(LOW,IAXIS)+guard(IAXIS)-1,:,:) = 0.0
        solnVec(varGrid,blkLimitsGC(HIGH,IAXIS)-guard(IAXIS)+1:blkLimitsGC(HIGH,IAXIS),:,:) = 0.0

        if(NDIM >= 2) then
           solnVec(varGrid,:,blkLimitsGC(LOW,JAXIS):blkLimitsGC(LOW,JAXIS)+guard(JAXIS)-1,:) = 0.0
           solnVec(varGrid,:,blkLimitsGC(HIGH,JAXIS)-guard(JAXIS)+1:blkLimitsGC(HIGH,JAXIS),:) = 0.0
        end if

        if(NDIM == 3) then
           solnVec(varGrid,:,:,blkLimitsGC(LOW,KAXIS):blkLimitsGC(LOW,KAXIS)+guard(KAXIS)-1) = 0.0
           solnVec(varGrid,:,:,blkLimitsGC(HIGH,KAXIS)-guard(KAXIS)+1:blkLimitsGC(HIGH,KAXIS)) = 0.0
        end if
     else
        solnVec(varGrid,:,:,:) = 0.0
     end if
     call Grid_releaseBlkPtr(gr_ptblkList(blkNo),solnVec,CENTER)
  end do



  ! The Grid_mapParticlesToMesh routine requires that
  ! the particles are sorted in block order
  call Grid_sortParticles(particles, part_props, localNumParticles, &
       particleTypes, numParticles, particlesPerBlk,BLK_PART_PROP)
#ifdef DEBUG_ALL
  print*,'On', gr_meshMe,', sum(particlesPerBlk)=',sum(particlesPerBlk)
#endif

  call gr_ensureValidNeighborInfo(10) ! We want valid grid information as after a guardcell-filling operation

  !---------------------------------------------------------------------------------------
  ! In order to determine the size of the send / receive buffer, we will count the number 
  ! of cells that exist off-processor.  
  ! Whilst doing this, we will keep track of all the neighbors in a temporary data structure
  ! named gr_ptDomain.
  !---------------------------------------------------------------------------------------
  if (gr_ptSmearLen > 0) then
     allocate(gr_ptDomain(gr_ptBlkCount), STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("Severe error. Memory cannot be allocated!")
     end if

     !Subroutine modifies module data structure named gr_ptDomain.
     call gr_ptStoreOffBlockCells(particlesPerBlk, gr_ptBlkList, gr_ptBlkCount, &
          blkLimitsGC, blkSize, guard, BufferSize)


     !Each process has the same value for BufferSize (global maximum buffer 
     !size).  If it is zero then no data needs communicating.  Set it to 1 
     !so there is buffer space for gr_ptMoveMappedData to send / recv a single 
     !NULL element.
     if (BufferSize <= 0) then
        BufferSize = 1
        if (gr_meshMe == 0) then
           print *, "[Grid_mapParticlesToMesh]: No off-processor mapping required"
        endif
     end if  


     allocate(sendBuf(BufferSize), recvBuf(BufferSize), STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("Severe error. Memory cannot be allocated!")
     end if
  end if



  ! ---------------------------------------------------------------------------------
  ! Now we can actually map the particles to the mesh.
  ! ---------------------------------------------------------------------------------
  allocate(gr_ptBuf(blkSizeGC(IAXIS),blkSizeGC(JAXIS),blkSizeGC(KAXIS)), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("Severe error. Memory cannot be allocated!")
  end if


  sendBufPtr=1  !Pointer to the first element of sendBuf.
  sendCount=0
  pEnd=0


  do blkNo = 1, gr_ptBlkCount

     blockID = gr_ptBlkList(blkNo)
     if (particlesPerBlk(blockID,particleTypes) > 0) then

        call Grid_getBlkPtr(blockID,solnVec,CENTER)
        gr_ptBuf = 0.0

        pStart=pEnd+1
        pEnd=pEnd+particlesPerBlk(blockID,particleTypes)

        call Particles_mapToMeshOneBlk(blkLimitsGC,guard,blockID,&
             particles(:,pStart:pEnd),pEnd-pStart+1,propPart,gr_ptBuf,&
             particleOffset=ptInfo+pStart-1)

        !Apply boundary conditions (BCs) to the guard cells of this block.
        !If for example, reflecting BCs are in place, a mass accumulated 
        !guard cell may map onto an internal cell of the SAME block.  This 
        !is NOT concerned with the global block configuration.  That is 
        !handled in gr_ptSameProcMap & gr_ptOffProcMap, where by periodic 
        !BCs will map a guard cell onto an internal cell of a DIFFERENT block.

        call gr_ptApplyBCsOneBlk(blkLimits,blkLimitsGC,blockID)

        !! All smearing within the block copied in 
        !! We have to add the gr_ptBuf contribution to what is already 
        !! there because an earlier iteration of blkNo could have written
        !! to this block.
        !! solnVec(varGrid,:,:,:) = solnVec(varGrid,:,:,:) + gr_ptBuf(:,:,:)
        !! Actually, only need internal region...
        solnVec(varGrid,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
             solnVec(varGrid,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) + &
             gr_ptBuf(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))

        call Grid_releaseBlkPtr(blockID,solnVec,CENTER)



        if (gr_ptSmearLen > 0) then

           !! For each section of the gcells find the negh. If the negh is
           !! on the local processor at the same reflevel, copy the 
           !! values into the block. If the negh is there, but at different
           !! reflevel restrict or prolong the values
           !! If the negh is not on the processor, put it in the sendBuf
           !! Even in sendbuf, there can be two ways of putting info.
           !! The destination might be known before sending, or it may 
           !! not be known, both situations have to be handled.

           ! Initialise certain data for the case when NDIM < MDIM.
           ! Required because gr_ptDomain data structure only contains data for NDIM.
           srcCoords(LOW:HIGH, NDIM:MDIM) = 1
           destCoords(LOW:HIGH, NDIM:MDIM) = 1
           neghCornerID(NDIM:MDIM) = 0

           ! Loop over each guard cell region in a block.
           do regionIter = 1, NUMBGUARDREGIONS

              numNegh = gr_ptDomain(blkNo) % haloRegion(regionIter) % numNegh

              if (numNegh > 0) then
                 do n = 1, numNegh

                    negh(:) = gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % negh(:)
                    neghCornerID(1:NDIM) = &
                         gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % cornerID(1:NDIM)
                    srcCoords(LOW:HIGH,1:NDIM) = &
                         gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % srcCoords(LOW:HIGH,1:NDIM)
                    destCoords(LOW:HIGH,1:NDIM) = & 
                         gr_ptDomain(blkNo) % haloRegion(regionIter) % neighbor(n) % destCoords(LOW:HIGH,1:NDIM)

                    ib=srcCoords(LOW,IAXIS); ie=srcCoords(HIGH,IAXIS)
                    jb=srcCoords(LOW,JAXIS); je=srcCoords(HIGH,JAXIS)
                    kb=srcCoords(LOW,KAXIS); ke=srcCoords(HIGH,KAXIS)

                    !Move onto the next guard cell section if the current sections contains no data.
                    !Assuming each guard cell region contains at least some data, then
                    !at least one processor will use all of the space in sendBuf.
                    if(maxVal(abs(gr_ptBuf(ib:ie,jb:je,kb:ke)))/=0.0) then

                       if(negh(BLKPROC)==gr_meshMe) then
                          call gr_ptSameProcMap(srcCoords,destCoords,negh,varGrid)
                       else
                          !print *, "Processor", gr_meshMe, "calling gr_ptOffProcMap.  CornerID=", &
                          !neghCornerID(:,n)
                          call gr_ptOffProcMap(srcCoords,destCoords,BufferSize,&
                               sendBuf,sendCount,sendBufPtr,negh,neghCornerID)
                          !print *, "Processor", gr_meshMe, "packed", sendCount, & 
                          !"reals out of", BufferSize, "space, which will be sent to processor", negh(BLKPROC,n)
                       end if
                    end if

                 end do
              end if
           end do

        end if  !End -> if (gr_ptSmearLen > 0).
        

     end if !Test value in particlesPerBlk
  end do !End loop over blocks.



  deallocate(gr_ptBuf)

  if (gr_ptSmearLen > 0) then

     if(sendCount > BufferSize) then
        call Driver_abortFlash("Severe error. Communication buffer too small!!!!")
     end if
     
     call gr_ptMoveMappedData(varGrid,BufferSize,sendBuf,sendCount,recvBuf)

     deallocate(gr_ptDomain)  !Used in gr_ptMoveMappedData when it calls gr_ptDumpState.
     deallocate(sendBuf)
     deallocate(recvBuf)
  end if

  !print *, "Processor", gr_meshMe, "finished in Grid_mapParticlesToMesh"

end subroutine Grid_mapParticlesToMesh
