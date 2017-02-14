!!****if* source/Grid/GridParticles/GridParticlesMapFromMesh/Grid_mapMeshToParticles
!!
!! NAME
!!
!!  Grid_mapMeshToParticles
!!
!! SYNOPSIS
!!
!!  Grid_mapMeshToParticles(real(inout) :: particles(part_props,numParticles),
!!                          integer(in) :: part_props,
!!                          integer(in) :: part_blkID,
!!                          integer(in) :: numParticles,
!!                          integer(in) :: posAttrib(POS),
!!                          integer(in) :: numAttrib,
!!                          integer(in) :: attrib(:,numAttrib),
!!                          integer(in) :: mapType,
!!               optional,  integer(in) :: gridDataStruct)
!!
!! DESCRIPTION
!!
!!  Routine to map a quantity defined on the mesh onto the particle positions.
!!
!! ARGUMENTS
!!
!!     particles:     Data structure containing particles information
!!     part_props : number of particle attributes
!!     part_blkID : the index of particle attributes that carries the block number
!!     numParticles : the number of particles on my proc
!!     posAttrib           : particles data structure indices pointing
!!                           to particle positions
!!     numAttrib           : number of attributes that need to be mapped
!!     attrib              : list of attributes and their corresponding
!!                           mesh data structure indices
!!                           processor
!!     mapType : method for mapping grid quantities to particles
!!     gridDataStruct: The Grid data structure that varGrid refers to; one of
!!                    CENTER, FACEX, FACEY, FACEZ, GRIDVAR.  
!!                    If this argument is not present, the default is CENTER, 
!!                    so that attrib indicates a list of variable in UNK.
!!
!!  
!!***
!*******************************************************************************

subroutine Grid_mapMeshToParticles (particles, part_props,part_blkID,&
                                    numParticles,posAttrib,&
                                    numAttrib, attrib,&
                                    mapType,gridDataStruct)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr,&
       Grid_getLocalNumBlks, Grid_getBlkBoundBox, Grid_getDeltas
  use Particles_interface, ONLY : Particles_mapFromMesh

  implicit none

#include "Flash.h"
#include "constants.h"
#include "GridParticles.h"
#include "Particles.h"


  integer, INTENT(in) :: part_props, numParticles, numAttrib, part_blkID
  real, INTENT(inout),dimension(part_props,numParticles) :: particles
  integer,dimension(MDIM), intent(IN) :: posAttrib
  integer, dimension(PART_ATTR_DS_SIZE,numAttrib),INTENT(in) :: attrib
  integer, INTENT(IN) :: mapType
  integer, optional, intent(IN) :: gridDataStruct

  integer :: i,j,k,currentBlk,blkCount,inBlk, prevBlk

  real, pointer, dimension(:,:,:,:) :: solnVec
  real, dimension(LOW:HIGH,MDIM) :: bndBox
  real,dimension(MDIM) :: delta, pos
  integer :: gDataStruct
  real, dimension(numAttrib) :: partAttribVec

  if(present(gridDataStruct)) then
     gDataStruct=gridDataStruct
  else
     gDataStruct=CENTER
  end if

  if(numParticles>0) then

     call Grid_getLocalNumBlks(blkCount)

     currentBlk=int(particles(part_blkID,1))
     prevBlk=currentBlk
     call Grid_getBlkPtr(currentBlk,solnVec,gDataStruct)
     call Grid_getBlkBoundBox(currentBlk,bndBox)
     call Grid_getDeltas(currentBlk,delta)
     do i = 1, numParticles
#ifdef DEBUG_GRIDPARTICLES
        if((particles(part_blkID, i) < 0) .or. (particles(part_blkID, i) > blkCount)) then
           call Driver_abortFlash("BLK_PART_PROP out of bounds")
        end if
#endif
        
        currentBlk=int(particles(part_blkID,i))
        if(currentBlk /= prevBlk)then
           call Grid_releaseBlkPtr(prevBlk,solnVec,gridDataStruct)
           call Grid_getBlkPtr(currentBlk,solnVec,gridDataStruct)
           call Grid_getBlkBoundBox(currentBlk,bndBox)
           call Grid_getDeltas(currentBlk,delta)
        end if
        do j = 1,MDIM
           pos(j)=particles(posAttrib(j),i)
        end do
        
        call Particles_mapFromMesh (mapType, numAttrib, attrib,&
             pos, bndBox,delta,solnVec, partAttribVec)
        do j = 1,numAttrib
           particles(attrib(PART_DS_IND,j),i)=partAttribVec(j)
        end do
        prevBlk=currentBlk
     enddo
     call Grid_releaseBlkPtr(currentBlk,solnVec,gridDataStruct)
  end if

end subroutine Grid_mapMeshToParticles

