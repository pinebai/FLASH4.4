!!****if* source/Grid/GridParticles/GridParticlesMove/UG/Directional/gr_ptMoveOffProc
!!
!! NAME
!!
!!  gr_ptMoveOffProc
!!
!! SYNOPSIS
!!
!!   gr_ptMoveOffProc(integer(IN)   :: face, 
!!                   integer(IN)   :: axis,
!!                   integer(IN)   :: index, 
!!                   integer(IN)   :: propCount,
!!                   integer(IN)   :: maxPerProc,
!!                   logical(IN)   :: boundary,
!!                   integer(IN)   :: lnegh,
!!                   integer(IN)   :: rnegh,
!!                   real(IN)      :: corner(LOW:HIGH),
!!                 integer(inout)  :: localNum,
!!                   real,(inout)  :: particles(:,:))
!!
!! DESCRIPTION
!!     
!!    This routine is used in moving the particles data when a particle moves
!!    off a block because of time integration. The routine examines a single face
!!    of the block at a time. Here every particle currently 
!!    on the processor is examined to see if it has moved out of the block across
!!    the specific face and axis under consideration.
!!    If the face was along physical boundary, then it may either leave the domain, 
!!    stay on the same block or be moved  to gr_ptDestBuf, 
!!    which holds particles to be passed to the next processor, depending
!!    on the boundary conditions. Once all particles have been examined, the contents
!!    of gr_ptDestBuf are sent to the neighbor along the face, and the contents
!!    of opposite face neighbors gr_ptDestBuf are received. The newly received particles
!!    are added to the particles data structure. Note that it is assumed that no particle
!!    can move across a full block in a single time step.
!!
!! ARGUMENTS
!!
!!     face - valid values are LOW and HIGH, indicating whether it is lower face or
!!            the upperface along the specified dimension
!!     axis - the physical dimension under consideration
!!     index - index into the particles data structure indicating the particle position 
!!             specified axis
!!     propCount - the count of particle properties
!!     maxPerProc - the maximum number of particles allowed on one processor
!!     boundary - indicates whether the face under consideration is on physical boundary
!!     lnegh     - neghbor to send data to
!!     rnegh     - neghbor to receive data from
!!     corner   - the coordinate of the lower upper face along the axis
!!     localNum -   The number of particles in the particle array that
!!                           are stored on this processor.
!!     particles -           A 2 dimensional array of particles and their property
!!                           values
!!
!!***

subroutine gr_ptMoveOffProc(face,axis,index, propCount, maxPerProc,boundary,&
     lnegh,rnegh,corner,localNum,particles)
#include "constants.h"
#include "Flash.h"

  use Grid_data, ONLY : gr_axisComm,gr_domainBC
  use gr_ptData, ONLY : gr_ptDestBuf, gr_ptSourceBuf,gr_ptBlk
  use gr_ptInterface, ONLY : gr_ptOneFaceBC
  use Driver_interface, ONLY: Driver_abortFlash
  implicit none
  include "Flash_mpi.h"


  integer, intent(IN) :: face, axis, index, propCount, maxPerProc
  logical,intent(IN) :: boundary
  integer, intent(IN) :: lnegh,rnegh
  real,dimension(LOW:HIGH),intent(IN) :: corner
  integer,intent(INOUT) :: localNum
  real,dimension(propCount,maxPerProc),intent(INOUT) :: particles

!!  real, allocatable, dimension(:,:) :: sendBuf,recvBuf

  integer :: j,k,numDest, lostParticles, count

  integer :: pid, pend,bufsize
  integer :: recvCount, tag, ierr
  integer,dimension(MPI_STATUS_SIZE) :: status
  real :: nothing
  logical :: moved

  numDest = 0
  lostParticles=0
  count=0
  tag=20*axis
  nothing = 0.0
  if(face>HIGH)call Driver_abortFlash("Grid_moveParticles: face value is not LOW/HIGH")
  
  pend=localNum
  j=1
  do pid = 1,localNum
     if(particles(gr_ptBlk,j)>0.0)then  !! verify that it is a valid
        !! particle before processing it
        
        moved = (face==LOW).and.(particles(index,j)<corner(face))
        moved = moved.or.((face==HIGH).and.(particles(index,j)>corner(face)))
        
        if(moved) then !! find if this particle moves

           if(boundary) then                       !! if it is on physical boundary
              lostParticles=0
              call gr_ptOneFaceBC(particles(:,j),propCount, axis,face,int(particles(gr_ptBlk,j)), &
                   &              lostParticles)
           end if

           moved = (particles(index,j)<corner(face)).or.(particles(index,j)>corner(face))

           if(moved) then
              numDest=numDest+1
              gr_ptDestBuf(:,numDest)=particles(:,j)
              particles(gr_ptBlk,j)=NONEXISTENT  !! make the particle invalid
           end if
        end if
        if(particles(gr_ptBlk,j)==NONEXISTENT) then
           if(j<pend)particles(:,j)=particles(:,pend)
           pend=pend-1
        else
           j=j+1
        end if
     end if
  end do

  localNum=pend
  count=numDest*propCount
  if(count==0)then   !! if there were no particles to send to the neighbor, 
     count=1         !! send a single data item so sendreceive doesn't hang
     gr_ptDestBuf(1,1)=nothing
  end if

  bufsize=propCount*maxPerProc

  recvCount=bufSize


  call MPI_Sendrecv(gr_ptDestBuf,count,FLASH_REAL,lnegh,tag,&  !! send to left neghbor, receive from 
       gr_ptSourceBuf,recvCount,FLASH_REAL,rnegh,tag,&           !! the right one
       gr_axisComm(axis),status,ierr)
  call MPI_Get_count(status,FLASH_REAL,recvCount,ierr)    !! find out how many particle received

  if(recvCount>1)then       !! if some particles were received, then append them to
     k=recvCount/propCount
     do j = 1,k
        particles(:,localNum+j)=gr_ptSourceBuf(:,j)
     end do
     localNum=localNum+k
  end if
end subroutine gr_ptMoveOffProc
