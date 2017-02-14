!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/UG/Grid_mapParticlesToMesh
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
!!  propPart:   Index of particle attribute to interpolate onto the mesh
!!  varGrid:    Index of gridded variable to receive interpolated
!!                         quantity
!!  mode:       (Optional) If zero (default), zero varGrid first; if
!!                              nonzero, do not zero varGrid first
!!  ptInfo:     (Optional) additional info to pass to the Particles unit
!!                              in some calls, for debugging
!!
!! PARAMETERS
!! 
!!***

!!REORDER(4): solnVec


#include "constants.h"

subroutine Grid_mapParticlesToMesh (particles,part_props,numParticles,&
                                    maxParticlesPerProc,&
                                    propPart, varGrid, mode, ptInfo)
  
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_sortParticles, Grid_getBlkBoundBox
#include "Flash.h"

  use gr_ptData, ONLY : gr_ptDestBuf, gr_ptSourceBuf,gr_ptBlkList, gr_ptBlkCount,&
       gr_ptBuf
  use Particles_interface, ONLY : Particles_mapToMeshOneBlk
  use Grid_data, ONLY :  gr_axisNumProcs, gr_meshMe
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY: Logfile_stampMessage
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_ptInterface, ONLY : gr_ptApplyBCsOneBlk,gr_ptExchangePartialMap

  implicit none
 
#include "Flash_mpi.h"

  integer,intent(IN) :: numParticles, part_props,maxParticlesPerProc
  real,dimension(part_props,maxParticlesPerProc),intent(INOUT) :: particles
  integer, INTENT(in) :: propPart, varGrid
  integer, INTENT(in), optional :: mode
  integer, INTENT(in), optional :: ptInfo

  integer,dimension(MAXBLOCKS)::ParticlesPerBlk
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer :: blockID,pStart,pEnd, vmode
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer,dimension(MDIM) :: blkSize,blkSizeGC,guard
  logical :: fcBdry=.false.
  real, allocatable, dimension(:) :: sendBuf,recvBuf
  integer :: axis, bufSize
  integer :: particleTypes=1
  real,dimension(LOW:HIGH,MDIM) :: bndBox
  integer :: localNumParticles

  localNumParticles=numParticles

  if (present(mode)) then
     vmode = mode
  else
     vmode = 0
  end if

  call Grid_sortParticles(particles, part_props, localNumParticles, &
       particleTypes,maxParticlesPerProc, ParticlesPerBlk,BLK_PART_PROP)

  blockID=1

  call Grid_getBlkPtr(blockID,solnVec,CENTER)
  call Grid_getBlkIndexLimits(blockID,blkLimits, blkLimitsGC)
  blkSizeGC=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
  blkSize=blkLimits(HIGH,:)-blkLimits(LOW,:)+1
  allocate(gr_ptBuf(blkSizeGC(IAXIS),blkSizeGC(JAXIS),blkSizeGC(KAXIS)))
  gr_ptBuf = 0.0

  if(vmode/=0.0) then
     !Copy the central section into gr_ptBuf.
     !This is equivalent to zeroing the guard cells.
     gr_ptBuf(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
          solnVec(varGrid,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) 
  end if

  pStart=1
  pEnd=particlesPerBlk(blockID)
  guard=blkLimits(LOW,:)-blkLimitsGC(LOW,:)


  call Particles_mapToMeshOneBlk(blkLimitsGC,guard,blockID,&
       particles(:,pStart:pEnd),pEnd-pStart+1,propPart,gr_ptBuf,&
       particleOffset=ptInfo)


  bufSize=minval(blkSizeGC)
  bufSize=(blkSizeGC(IAXIS)*blkSizeGC(JAXIS)*blkSizeGC(KAXIS))/bufSize
  bufSize=bufSize*2
  allocate(recvBuf(bufSize))
  allocate(sendBuf(bufSize))


  !Apply boundary conditions to the cells in this block.
  !i.e. consider the data exchange with respect to this single block.
  call gr_ptApplyBCsOneBlk(blkLimits, blkLimitsGC, blockID)


  !Now consider the data exchange between the single block and the rest of the 
  !computational domain.
  do axis = IAXIS,NDIM

     !There must be more than one processor in a single dimension for the 
     !guard cell exchange to make sense.  (N.B. If there is 
     !a single block in a domain with periodic BCs, the data exchange  
     !was handled previously by gr_ptApplyBCsOneBlk.)

     if (gr_axisNumProcs(axis) > 1) then
        call gr_ptExchangePartialMap(blkLimits,blkLimitsGC,bufSize,axis,LOW,sendBuf,recvBuf)
        call gr_ptExchangePartialMap(blkLimits,blkLimitsGC,bufSize,axis,HIGH,sendBuf,recvBuf)
     end if

  end do


  !We overwrite rather than accumulate because we have already 
  !accumulated the central section when vmode/=0.
  solnVec(varGrid,:,:,:) = gr_ptBuf
  deallocate(gr_ptBuf)
  deallocate(sendBuf)
  deallocate(recvBuf)
  call Grid_releaseBlkPtr(blockID,solnVec,CENTER)

end subroutine Grid_mapParticlesToMesh
