!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/UG/gr_ptExchangePartialMap
!!
!! NAME
!!  gr_ptExchangePartialMap
!!
!! SYNOPSIS
!!
!!   gr_ptExchangePartialMap(integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits
!!                           integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
!!                           integer,intent(IN) :: bufSize
!!                           integer,intent(IN) :: axis
!!                           integer,intent(IN) :: face
!!                           real,dimension(bufSize),intent(INOUT) :: sendBuf
!!                           real,dimension(bufSize),intent(INOUT) :: recvBuf
!!
!! DESCRIPTION
!!
!! Routine to send the guard cells on each block to its neighboring block 
!! for a particular axis and face.  All processors along a certain axis must 
!! participate.  Periodic boundary conditions are assumed, so as to prevent deadlock. 
!! This does not pose a problem to non-periodic simulations because the guard cells 
!! are zero in these simulations.  As such, accumulation of guard cells 
!! of zero value does not change the original cell value.
!!
!! ARGUMENTS
!!
!!                           blkLimits: A block's internal cell range.
!!                           blkLimitsGC: A block's complete cell range including guard cells.
!!                           bufSize: Size of the send / receive buffer.
!!                           axis: Current axis of communication.
!!                           face: Current direction of communication along a certain axis.
!!                           sendBuf: Send buffer for communication.
!!                           recvBuf: Receive buffer for communication.
!!
!!***

subroutine gr_ptExchangePartialMap(blkLimits,blkLimitsGC,bufSize,axis,face,&
     sendBuf,recvBuf)

  use Grid_data, ONLY : gr_axisMe, gr_axisNumProcs, gr_axisComm, gr_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_ptData, ONLY : gr_ptBuf
  use gr_ptMapData, ONLY : gr_ptSmearLen
  use Grid_interface, ONLY : Grid_getBlkBC
  use ut_conversionInterface, ONLY : ut_convertToMemoryOffset, ut_convertToArrayIndicies

  implicit none  

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits,blkLimitsGC
  integer,intent(IN) :: bufSize,axis,face
  real,intent(inout),dimension(bufSize):: sendBuf,recvBuf

  integer,dimension(MDIM) :: blkSize,blkSizeGC
  integer, dimension(LOW:HIGH, MDIM) :: faces
  integer,dimension(MDIM) :: beg,fin,addFact,elementCoord
  integer :: i,j,k,n, snegh,rnegh,tag,ierr,nRecd
  integer,dimension(MPI_STATUS_SIZE) :: status
  integer :: i1, j1, k1, smearLength, memoryOffset

  smearLength = gr_ptSmearLen
  blkSizeGC=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
  blkSize=blkLimits(HIGH,:)-blkLimits(LOW,:)+1

#ifdef DEBUG_GRIDMAPPARTICLES
  if((NDIM==1).and.(axis/=IAXIS))then
     print*,"[gr_ptExchangePartialMap]: Axis=",axis
     call Driver_abortFlash("[gr_ptExchangePartialMap]: problem is one dimensional")
  end if
  if((NDIM==2).and.(axis==KAXIS))then
     print*,"[gr_ptExchangePartialMap]: Axis=",axis
     call Driver_abortFlash("[gr_ptExchangePartialMap]: problem is two dimensional")
  end if
  print*,"The blocksize sent in and face is", blkSize, face
#endif


  if (smearLength < 0) then
     call Driver_abortFlash("[gr_ptExchangePartialMap]: smearLength is less than 0!")
  else if (smearLength == 0) then
     !Nothing to exchange.
     return
  end if


  !-----------------------------------------------------------------------
  !Determine the range of guard cells whose value needs to be tested.
  !Remember, there are NDIM rounds of communication.
  !-----------------------------------------------------------------------
  beg = 1
  fin = 1
  addfact = 0

  if(axis==IAXIS) then

     do i=JAXIS,NDIM
        beg(i)=blkLimits(LOW,i)-smearLength
        fin(i)=blkLimits(HIGH,i)+smearLength
     end do

  elseif(axis==JAXIS) then

     beg(IAXIS)=blkLimits(LOW,IAXIS)
     fin(IAXIS)=blkLimits(HIGH,IAXIS)
     if(NDIM>2) then
        beg(KAXIS)=blkLimits(LOW,KAXIS)-smearLength
        fin(KAXIS)=blkLimits(HIGH,KAXIS)+smearLength
     end if

  else

     do i=IAXIS,JAXIS
        beg(i)=blkLimits(LOW,i)
        fin(i)=blkLimits(HIGH,i)
     end do

  end if


  !We consider smearLength guard cells for the axis of interest.
  if(face==LOW) then
     beg(axis)=blkLimits(LOW,axis)-smearLength
     fin(axis)=blklimits(LOW,axis)-1
     addfact(axis)= blkSize(axis)
  else
     beg(axis)=blkLimits(HIGH,axis)+1
     fin(axis)=blklimits(HIGH,axis)+smearLength
     addfact(axis)= -blkSize(axis)
  end if


  !-----------------------------------------------------------------------
  !Determine whether any data must be packed into the send buffer.
  !-----------------------------------------------------------------------
  n=0

  call Grid_getBlkBC(1, faces)

  !There is only a valid neighbor under the following circumstances.
  if( (faces(face,axis)==NOT_BOUNDARY) .or. (faces(face,axis)==PERIODIC) ) then

     do k = beg(KAXIS),fin(KAXIS)
        k1 = k+addfact(KAXIS)
        do j = beg(JAXIS),fin(JAXIS)
           j1 = j+addfact(JAXIS)
           do i = beg(IAXIS),fin(IAXIS)
              i1=i+addfact(IAXIS)
              if(gr_ptBuf(i,j,k)/=0.0) then

                 elementCoord(IAXIS) = i1
                 elementCoord(JAXIS) = j1
                 elementCoord(KAXIS) = k1

                 call ut_convertToMemoryOffset(MDIM, elementCoord, blkLimitsGC(LOW,:), blkLimitsGC(HIGH,:), memoryOffset)
                 sendBuf(n+1) = memoryOffset
                 n=n+2
                 sendBuf(n)=gr_ptBuf(i,j,k)

              end if
           end do
        end do
     end do
  end if


  !This fragment is executed when there is a neighbor we want to send info. to, but no data.
  !It is also executed when we must involve a neighbor to prevent deadlock, even 
  !if we do not want to communicate.
  if(n==0)then
     n=1
     sendBuf(n)=0.0
  end if



  !-----------------------------------------------------------------------
  !Determine which processors we will communicate with .
  !The code is designed for PERIODIC boundaries.  When we use 
  !different boundary conditions no data can end up in the
  !guard cells, so we are just communicating to prevent deadlock.
  !-----------------------------------------------------------------------
  tag=20
  if(face==LOW) then
     snegh=gr_axisMe(axis)-1     !! determine the negh to send
     if(snegh<0) snegh=gr_axisNumProcs(axis)-1    !! left boundary

     rnegh=gr_axisMe(axis)+1    !! do the same for right neighbor 
     if(rnegh==gr_axisNumProcs(axis)) rnegh=0
  else
     rnegh=gr_axisMe(axis)-1     !! determine the negh to send
     if(rnegh<0) rnegh=gr_axisNumProcs(axis)-1    !! left boundary

     snegh=gr_axisMe(axis)+1    !! do the same for right neighbor 
     if(snegh==gr_axisNumProcs(axis)) snegh=0
  end if


#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
  print *, "Processor", gr_meshMe, "sending to processor", snegh, "receiving from processor", rnegh, &
       "along communicator axis:", axis, "sending", n, "elements."
#endif


  call MPI_Sendrecv(sendBuf,n,FLASH_REAL,snegh,tag,&  !! send to one negh and receive from the
       recvBuf,bufSize,FLASH_REAL,rnegh,tag,&   !! other one on the current axis
       gr_axisComm(axis),status,ierr)
  call MPI_Get_count(status,FLASH_REAL,nRecd,ierr) !! find number received



  !-----------------------------------------------------------------------
  !Extract any data from the receive buffer.
  !-----------------------------------------------------------------------
  if(nRecd>1) then

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
    print *, "Processor", gr_meshMe, "received", nRecd, "elements"
#endif

     do n = 1, nRecd, 2

        memoryOffset = int(recvBuf(n))
        call ut_convertToArrayIndicies(MDIM, memoryOffset, blkLimitsGC(LOW,:), blkLimitsGC(HIGH,:), elementCoord)

        i = elementCoord(IAXIS)
        j = elementCoord(JAXIS)
        k = elementCoord(KAXIS)

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
        print *, "Processor", gr_meshMe, "mapping to central element", i, j, k
#endif
        gr_ptBuf(i,j,k) = gr_ptBuf(i,j,k) + recvBuf(n+1)

     end do

  end if

  return
end subroutine gr_ptExchangePartialMap
