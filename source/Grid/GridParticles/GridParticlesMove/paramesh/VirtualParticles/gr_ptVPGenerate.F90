!!****if* source/Grid/GridParticles/GridParticlesMove/paramesh/VirtualParticles/gr_ptVPGenerate
!!
!! NAME
!!
!!  gr_ptVPGenerate
!!
!! SYNOPSIS
!!
!!  gr_ptVPGenerate(real(INOUT)   :: particle(:),
!!                    integer(IN)   :: propCount,
!!                   integer(OUT)   :: destCount, 
!!                      real(OUT)   :: destParticles,
!!                    integer(IN)   :: blockID,
!!                    logical(IN)   :: newBlkID)
!!
!! DESCRIPTION
!!     
!!    This routine is used in moving the non stationaly data elements 
!!    associated with structures like particles and ray, when a data element
!!    moves off a block without re-gridding. Here every element currently 
!!    on the processor is examined to see if it still belongs to the same block.
!!    If it does not, it is further examimned to see if it has moved out of the physical boundary.
!!    If is out of physical boundary, it may either leave the domain, stay on the same block
!!    or be moved  to destBuf, which holds elements to be passed to the next processor, depending
!!    on the boundary conditions. If it is still in the physical domain, it may have
!!    moved to another block on the same processor, in which case only its BLK
!!    needs to change, otherwise it is moved to destBuf.
!!
!! ARGUMENTS
!!
!!     particle -  single element with all its attributes
!!                           values
!!     propCount - number of element attributes
!!
!!     destCount - the count of virtual particles created
!!     destParticles - the newly created virtual particles
!!     blockID   - the block to which the current particle is attached
!!     newBlkID  - indicates if the blockid is new, and therefore
!!                 boundblock and surrblks are needed
!! 
!!***

#ifdef DEBUG_ALL
#define DEBUG_PARTICLES
#endif

subroutine gr_ptVPGenerate(particle,propCount,destCount,destParticles,&
                           blockID,newBlkID)
#include "constants.h"
#include "Flash.h"
#include "gr_ptMapToMesh.h"

  use gr_ptData, ONLY : gr_ptProc,gr_ptBlk,gr_ptPosx,gr_ptPosy,gr_ptPosz
  use gr_ptVPData, ONLY : gr_ptVPMaxCount
  use Grid_data,ONLY : gr_globalDomain, gr_domainBC

  implicit none
  integer, intent(IN) :: propCount
  real,dimension(propCount),intent(INOUT):: particle
  integer,intent(OUT) :: destCount
  real,dimension(propCount, gr_ptVPMaxCount), intent(OUT) :: destParticles
  integer,intent(IN) :: blockID
  logical, intent(IN) :: newBlkID

!!  integer :: blockID, procID
  integer,parameter :: MAXCOUNT=7
  integer,dimension(BLKNO:PROCNO,MAXCOUNT) :: neghID
  integer,dimension(BLKNO:PROCNO)::destNeghID
  integer,dimension(BLKID:REFLEVELDIF,ABSMAXNEGH):: negh
  integer,dimension(MDIM,ABSMAXNEGH) :: neghCornerID
  integer :: numNegh
  real,dimension(MDIM) ::  pos
  integer,dimension(MDIM) ::  onBoundary
  logical :: leftDomain
  integer :: vpartCount
  integer,dimension(IAXIS:VP_LEAVE,MAXCOUNT)::vparticles
  integer i,j
  real :: auxTag

  destCount = 0 
  call gr_ptVPBC(particle,propCount, leftDomain, onBoundary)

  if(leftDomain)then !! nothing needs to be done except 
                     !!letting the particle cease to exist
     particle(gr_ptBlk)=LOST
     destCount=-1
  else
     vpartCount=MAXCOUNT
     pos(1:MDIM)=particle(gr_ptPosx:gr_ptPosz)

     call gr_ptVPMatchCondition(pos, blockID, newBlkID, vpartCount,vparticles)
    
     if(vpartCount>0) then !! the particle needs handling 
        auxTag=particle(TAG_PART_PROP)
        !destCount = 0         !! DestCount initialized
        do i = 1,vpartCount
           destNeghID=0
           if(vparticles(VP_LEAVE,i)==1) then
              !! This part is necessary when the particle is in the center of a face
              !! and the neighbor along that face is refined. Then potentially there
              !! there is a possibility of four copies of the particles, and if it is
              !! leaving the current block we need to determine which of the four 
              !! neghbors is the destination
             
              call gr_findNeghID(blockID,pos,vparticles(IAXIS:KAXIS,i),destNeghID)

           end if

           !! And now we find out if we missed any neghbors is call to match condition
           !! because of fine-coarse boundary at this neighbor
           call gr_ptFindNegh(blockID,vparticles(IAXIS:KAXIS,i),neghid,neghCornerID,numNegh)

           !! If there was a fine-coarse boundary and the particle is sitting close to
           !! to the middle of that boundary numNegh will be greater than 1
           !! destCount = 0      
           do j = 1,numNegh
              destCount=destCount+1
              destParticles(:,destCount)=particle(:)

              if (onBoundary(IAXIS)==vparticles(IAXIS,i)) then
                 if ((onBoundary(IAXIS)==LEFT_EDGE).and.&
                      (gr_domainBC(LOW,IAXIS)==PERIODIC)) then
                    destParticles(gr_ptPosx,destCount)=particle(gr_ptPosx)+(gr_globalDomain(HIGH,IAXIS)-gr_globalDomain(LOW,IAXIS))
                 elseif ((onBoundary(IAXIS)==RIGHT_EDGE).and.&
                      (gr_domainBC(HIGH,IAXIS)==PERIODIC)) then
                    destParticles(gr_ptPosx,destCount)=particle(gr_ptPosx)-(gr_globalDomain(HIGH,IAXIS)-gr_globalDomain(LOW,IAXIS))
                 end if
              endif
              if (onBoundary(JAXIS)==vparticles(JAXIS,i)) then
                  if ((onBoundary(JAXIS)==LEFT_EDGE).and.&
                      (gr_domainBC(LOW,JAXIS)==PERIODIC)) then
                     destParticles(gr_ptPosy,destCount)=particle(gr_ptPosy)+(gr_globalDomain(HIGH,JAXIS)-gr_globalDomain(LOW,JAXIS))
                 elseif  ((onBoundary(JAXIS)==RIGHT_EDGE).and.&
                      (gr_domainBC(HIGH,JAXIS)==PERIODIC)) then
                     destParticles(gr_ptPosy,destCount)=particle(gr_ptPosy)-(gr_globalDomain(HIGH,JAXIS)-gr_globalDomain(LOW,JAXIS))
                 end if
              endif
#if NDIM==3
              if (onBoundary(KAXIS)==vparticles(KAXIS,i)) then
                  if ((onBoundary(KAXIS)==LEFT_EDGE).and.&
                      (gr_domainBC(LOW,KAXIS)==PERIODIC)) then
                     destParticles(gr_ptPosz,destCount)=particle(gr_ptPosz)+(gr_globalDomain(HIGH,KAXIS)-gr_globalDomain(LOW,KAXIS))
                 elseif ((onBoundary(KAXIS)==RIGHT_EDGE).and.&
                      (gr_domainBC(HIGH,KAXIS)==PERIODIC)) then
                     destParticles(gr_ptPosz,destCount)=particle(gr_ptPosz)-(gr_globalDomain(HIGH,KAXIS)-gr_globalDomain(LOW,KAXIS))
                  end if
              end if
#endif

              destParticles(PROC_PART_PROP,destCount)=neghid(PROCNO,j)
              destParticles(BLK_PART_PROP,destCount)=neghid(BLKNO,j)
              !! If the neghid found by gr_ptFindNegh matches with
              !! that found by gr_findNeghID then this is the real
              !! copy of the particle
              if((neghid(BLKNO,j)==destNeghID(BLKNO)).and.((neghid(PROCNO,j)==destNeghID(PROCNO))))then
                 particle(TAG_PART_PROP)=-particle(TAG_PART_PROP)
                 
#ifdef DEBUG_VPARTICLES
                 write(*,*) 'Real particle i=',i,'new TAG=',int(particle(TAG_PART_PROP))
#endif
                 
              else
                 destParticles(TAG_PART_PROP,destCount)=-auxTag

#ifdef DEBUG_VPARTICLES
                 write(*,*) destCount,'on blk=',int(destParticles(BLK_PART_PROP,destCount)),' TAG=', & 
                      destParticles(TAG_PART_PROP,destCount)
#endif
              end if
           end do
           
        end do
       
     end if
  end if
  !write(*,*)'Final no. of virtual particles created is= ',destcount
end subroutine gr_ptVPGenerate
