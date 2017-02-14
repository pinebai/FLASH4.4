!!****if* source/Particles/ParticlesMain/Particles_dump
!!
!! NAME
!!    Particles_dump
!!
!! SYNOPSIS
!!    Particles_dump()
!!                   integer(IN) :: blockCount,
!!                   integer(IN), dimension(:)  :: blockList,
!!                   integer(IN) :: nstep,
!!                   real(IN)    :: time)
!!
!! DESCRIPTION
!!    Dump particle information to a plain old file
!!    Quick and dirty for blue gene testing
!!
!! ARGUMENTS
!!
!!  
!!  blockCount:             integer(IN)     number of blocks within this processor
!!  blockList(blockCount):  integer(IN)     block IDs within this processor
!!  nstep:                  integer(IN)     current time step index
!!  time:                   real(IN)        current simulation time
!!
!!
!! PARAMETERS
!!
!!***

!**********************************************************************
subroutine Particles_dump(blockCount,blockList,nstep,time,dt)

  use Particles_data, ONLY:  pt_numLocal, particles,pt_meshMe, pt_meshComm, pt_meshNumProcs
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY :  Grid_getBlkBoundBox

  implicit none
  
#include "Flash.h"
#include "constants.h"
 include "Flash_mpi.h"
  

  
  
  integer, intent(IN) :: blockCount
  integer, intent(IN) :: blockList(blockCount)
  integer, intent(IN) :: nstep
  real, intent(IN)    :: time, dt

  integer ::  fileUnit = 2
  character(len=23)   :: filename
  integer,dimension(4) :: prNum
  integer :: temp, i, p,blockID
  real,dimension(2,MDIM) :: limit

  integer, allocatable           :: partsToLeft(:)
  integer           :: particleOffset, globalNumParticles, ierr,pos, fh
  real, allocatable, dimension(:) :: partBuf
  integer ::    blkID
  integer ::  lsize, gsize, type_subarray, count, disp, status(MPI_STATUS_SIZE)
  integer :: sizeSubArray, sizeReal

  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: localnp
  integer, save :: recNo=0

  allocate(partsToLeft(0:pt_meshNumProcs))
  recNo=recNo+1
  ! some trick here to get a filename that reflects the processor or time step
  temp = pt_meshMe
  call Grid_getBlkBoundBox(1,limit)
  filename = "test_ParticlesDump"

  localnp = pt_numLocal

  if(pt_numLocal > 0) then
     allocate(partBuf(0:(pt_numLocal*4 -1)))
  end if

  ! remember that you're writing out the velocity from the OLD time step
  !write(fileUnit,890)pt_meshMe,nstep,time,dt,pt_numLocal
  pos = 0

  do p = 1, pt_numLocal
     
     partBuf(pos) = particles(TAG_PART_PROP,p)
     partBuf(pos+1) = particles(BLK_PART_PROP,p)
     partBuf(pos+2) = particles(POSX_PART_PROP,p)
     partBuf(pos+3) = particles(POSY_PART_PROP,p)
     
!!     if(recNo==1)print*,int(particles(TAG_PART_PROP,p))

     blkID = int(particles(BLK_PART_PROP,p))
     call Grid_getBlkBoundBox(blkID,limit)
!!$     if(partBuf(pos)==1.2201220e+12) then
!!$        print*,'dump',pt_meshMe,int(particles(BLK_PART_PROP,p)),p
!!$        print*,'dump',particles(POSX_PART_PROP,p),particles(POSY_PART_PROP,p)
!!$        print*,'dump',limit(LOW:HIGH,IAXIS)
!!$        print*,'dump',limit(LOW:HIGH,JAXIS)
!!$     end if

     if(particles(POSX_PART_PROP,p)<limit(1,1))print*,'particle p = ',p,' out of limit lowerface iaxis', limit(1,1)
     if(particles(POSX_PART_PROP,p)>limit(2,1))print*,'particle p = ',p,' out of limit upperface iaxis',limit(2,1)
     if(NDIM > 1) then
        if(particles(POSY_PART_PROP,p)<limit(1,2))print*,'particle p = ',p,' out of limit lowerface jaxis',limit(1,2)
        if(particles(POSY_PART_PROP,p)>limit(2,2))print*,'particle p = ',p,' out of limit upperface jaxis',limit(2,2)
     end if
     if(NDIM>2)then
        if(particles(POSZ_PART_PROP,p)<limit(1,3))print*,'particle p = ',p,' out of limit lowerface kaxis',limit(1,3)
        if(particles(POSZ_PART_PROP,p)>limit(2,3))print*,'particle p = ',p,' out of limit upperface kaxis',limit(2,3)
     end if
     pos=pos+4
     
  enddo
  
!  close(fileUnit)

  globalNumParticles = 0

  !GET the particles offset - this could be a subroutine of its own
  !gather all local particles from all procs and put them into an array partsToLeft
  call MPI_Allgather(pt_numLocal, 1, MPI_INTEGER, &
       partsToLeft, 1, MPI_INTEGER, &
       pt_meshComm, ierr)
  

  do i = 0, pt_meshNumProcs-1
     globalNumParticles = globalNumParticles + partsToLeft(i)
  end do
  
  do i = pt_meshNumProcs-1,1,-1
     partsToLeft(i) = partsToLeft(i-1)
  end do
  
  partsToLeft(0) = 0

  do i = 2,pt_meshNumProcs-1
     partsToLeft(i) = partsToLeft(i) + partsToLeft(i-1)
  end do
  
  particleOffset = partsToLeft(pt_meshMe)

  offset = (particleOffset * 4  + ((nstep -1)*globalNumParticles*4))*8

  !print *, "info = ", particleOffset, offset, nstep, globalNumParticles, pt_meshMe


  call MPI_File_open(pt_meshComm,filename,&
          MPI_MODE_CREATE+MPI_MODE_RDWR,MPI_INFO_NULL,fh,ierr)
  if(ierr < 0) then
     call Driver_abortFlash("Particles_dump: Error opening the file")
  end if
  
  call MPI_File_set_view(fh,offset,FLASH_REAL,FLASH_REAL,"native",&
       MPI_INFO_NULL,ierr)
  if(ierr < 0) then
     call Driver_abortFlash("Particles_dump: Error setting the file view")
  end if

  if(pt_numLocal > 0) then
     call MPI_File_write(fh,partBuf,4*pt_numLocal,FLASH_REAL,status,ierr)
  end if
  if(ierr < 0) then
     call Driver_abortFlash("Particles_dump: Error writing the file")
  end if

  call MPI_File_close(fh,ierr)
  if(ierr < 0) then
     call Driver_abortFlash("Particles_dump: Error closing the file")
  end if

  if(pt_numLocal > 0) then
     deallocate(partBuf)
  end if

  deallocate(partsToLeft)

!890 format(I5,5X,I5,5X,G12.4,3X,G15.7,3X,I5)
!900 format(I5,3X,I6,2X,I5,3X,6(G12.5,3X))

  
end subroutine Particles_dump

