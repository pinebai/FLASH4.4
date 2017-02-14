!!****if* source/Grid/GridSolvers/Pfft/gr_pfftTranspose
!!
!! NAME
!!
!!  gr_pfftTranspose
!!
!! SYNOPSIS
!!
!!  call gr_pfftTranspose (integer(IN) :: dir,
!!                         integer(IN) :: baseDatType,
!!                         real(INOUT) :: inArray(:),
!!                         real([IN]OUT)   :: outArray(:),
!!                         integer(IN) :: len(MDIM) , 
!!                         integer(IN) :: lda,
!!                         integer(IN) :: pr,
!!                         integer(IN) :: comm)
!!
!!  call gr_pfftTranspose[3DArr] (integer(IN) :: dir,
!!                                integer(IN) :: baseDatType,
!!                                real(INOUT) :: inArray(:),
!!                                real([IN]OUT)   :: outArray(:),
!!                                integer(IN) :: len(MDIM) , 
!!                                integer(IN) :: lda,
!!                                integer(IN) :: pr,
!!                                integer(IN) :: comm)
!!
!! DESCRIPTION
!!
!! This routine carries out the distributed transpose of the data as follows:
!! the original local array size is <len(IAXIS),len(JAXIS),len(KAXIS)>
!! The correspoding global size of the domain is <len(IAXIS),len(JAXIS)*nproc(JAXIS),
!! len(KAXIS)*nproc(KAXIS)>, where nproc(IAXIS:JAXIS) contains the configuration
!! of the processor grid. nproc(IAXIS) is assumed to be 1.
!!
!! If the transpose direction is forward, then the resulting local and global sizes
!! respectively are : <lda, len(KAXIS),newlen(IAXIS)/nproc(JAXIS)> and 
!! <lda, len(KAXIS)*nproc(KAXIS),newlen(IAXIS)>, where lda is input giving the 
!! actual global size of the data along x-axis. if lda is not a multiple of 
!! nproc(JAXIS), then len(JAXIS) is calculated by rounding lda to the nearest higher
!! multiple of nproc(JAXIS), and then dividing by nproc(JAXIS). Similarly newlen(IAXIS)
!! is calculated by rounding len(IAXIS) to nearest higher multiple of nproc(JAXIS).
!!
!! If the transpose is in inverse direction then the resulting local and global 
!! distributions are : <lda,newlen(IAXIS)/nproc(KAXIS),len(JAXIS)> and
!! <lda,newlen(IAXIS),len(JAXIS)>.
!!
!! ARGUMENTS
!!
!!    dir      -  whether forward or inverse transpose
!!    baseDatType     -  whether real or complex data
!!    inArray  -  array containing the data to be transposed. It is 
!!                overwritten during the transpose
!!    outArray -  array containing transposed data
!!    len      -  integer array containing the dimensions of the dataset
!!    lda      -  datasize in first dimension of outArray.
!!    pr       -  The number of processors in the communicator
!!    comm     -  The communicator that includes all the processors along
!!                the row or the column of the grid, that participate in the
!!                distributed transpose. 
!!
!!
!! NOTES
!!  This routine can handle real and complex data. If the data are complex
!!  then the routine assumes that len(IAXIS) is actually the complex vector
!!  size, and multiplies it by 2 internally. Similarly the returned value
!!  of len(IAXIS) is also that of a complex vector.
!!
!!***

#include "Pfft.h"
#include "constants.h"

subroutine gr_pfftTranspose(dir,baseDatType,inArray,outArray,inlen,outlen,pr,comm)
  implicit none

  integer, intent(IN) :: pr,comm,dir,baseDatType
  integer, dimension(MDIM), intent(IN) :: inlen,outlen
  real, dimension(:),intent(INOUT) :: inArray,outArray

  call gr_pfftTransposeImplementation(dir,baseDatType,&
       inArray,size(inArray),&
       outArray,size(outArray),&
       inlen,outlen,pr,comm)

  return
end subroutine gr_pfftTranspose



subroutine gr_pfftTranspose3DArr(dir,baseDatType,inArray,outArray,inlen,outlen,pr,comm)
  implicit none

  integer, intent(IN) :: pr,comm,dir,baseDatType
  integer, dimension(MDIM), intent(IN) :: inlen,outlen
  real,dimension(:,:,:),intent(INOUT) :: inArray
  real,dimension(:,:,:),intent(INOUT) :: outArray

  call gr_pfftTransposeImplementation(dir,baseDatType,&
       inArray,size(inArray),&
       outArray,size(outArray),&
       inlen,outlen,pr,comm)

  return
end subroutine gr_pfftTranspose3DArr





! Do not make an explicit interface for the following subroutine visible
! to the wrappers above - falling back to Fortran77-like level of knowledge
! is the whole point here! - KW
subroutine gr_pfftTransposeImplementation(dir,baseDatType,&
     inArray,inArraySize,&
     outArray,outArraySize,&
     inlen,outlen,pr,comm)
  use gr_pfftInterface, ONLY : gr_pfftLocalTranspose, gr_pfftIntersperse,&
        gr_pfftDisperse
  use gr_pfftData, ONLY : pfft_ndim,pfft_me
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Flash_mpi.h"

  integer, intent(IN) :: pr,comm,dir,baseDatType
  integer, dimension(MDIM), intent(IN) :: inlen,outlen
  integer, intent(IN) :: inArraySize, outArraySize
  real, dimension(inArraySize),intent(INOUT) :: inArray
  real, dimension(outArraySize),intent(INOUT) :: outArray
  
  integer ::  nx,nx1,npx,nPerProc
  integer ::  ierr, factor, numbCommProcs

!! If the transform is real to complex, then the data arriving
!! into the transpose is complex, but if it is sine or cosine transform
!! thent the data is real. Data manipulation needs to know if it should
!! treat the data as complex or real
!! nx1 is the real length of the data size in the second dimension
!! It is an input to account for the fact, that ny*pr could be > nx1
!! if nx1 is not an exact multiple of pr
  
  factor=1
  if((baseDatType==PFFT_PCLDATA_COMPLEX) .or. &
       (baseDatType==PFFT_PCLDATA_COMPLEX_STUFFED) .or. &
       (baseDatType==PFFT_PCLDATA_COMPLEX_EXTENDED)) &
       factor=2

  if(dir==PFFT_FORWARD) then

     !! Please see pfft_init for explanation of inLen and outLen sizes

     nx=inLen(IAXIS)
     nx1=inLen(JAXIS)*inLen(KAXIS)

     call gr_pfftLocalTranspose(inArray,outArray,nx,nx1,baseDatType)
     !! Here the shape of outArray is effectively
     !! inLen(JAXIS), inLen(KAXIS), inLen(IAXIS)

     npx=outLen(pfft_ndim) !! it is always the last dimension that is divided.

     if(pr>1) then

        nPerProc = nx1*npx*factor
        call MPI_Comm_size(comm, numbCommProcs, ierr)       
        
        if((nPerProc * numbCommProcs > size(inArray)).or.&
             (nPerProc * numbCommProcs > size(outArray))) then
 
           print *, "Chunk size:", nPerProc * numbCommProcs, "number of processors:", numbCommProcs, &
                "size inArray:", size(inArray), "size(outArray):", size(outArray) 
           call Driver_abortFlash("Severe error: array will be accessed out of bounds(1)... exiting")           
        end if
           

        call MPI_ALLTOALL(outArray,nPerProc,FLASH_REAL,inArray,nPerProc,&
                          FLASH_REAL,comm,ierr)
        !!   At this point the array is of the form 
        !!   (inLen(JAXIS),inLen(KAXIS),outLen(KAXIS),pr) 
        !!   so that what was second dimension is now all in processor
        !!   and what was the first dimension gets distributed. 
        !!   The third dimension becomes the second dimension,
        !!   but its distributions remains unchanged.

        nx=inLen(JAXIS)*factor
        nx1=inLen(KAXIS)*npx
        npx = outLen(IAXIS)*factor
        call gr_pfftIntersperse(inArray,outArray,nx,nx1,npx)
        !! Now the outArray shape is outLen(IAXIS),outLen(JAXIS),
        !! outLen(KAXIS)
     end if
  else
     nx=inLen(IAXIS)*factor
     npx=outLen(JAXIS)*factor
     nx1=inLen(JAXIS)*inLen(KAXIS)

     if(pr > 1) then
        call gr_pfftDisperse(inArray,outArray,nx,nx1,npx)

        !! the shape of outArray is (outLen(JAXIS),outLen(KAXIS),
        !! inLen(KAXIS),pr)

        nPerProc=npx*nx1

        call MPI_Comm_size(comm, numbCommProcs, ierr)       
        
        if((nPerProc * numbCommProcs > size(inArray)).or.&
             (nPerProc * numbCommProcs > size(outArray))) then
 
           print *, "Chunk size:", nPerProc * numbCommProcs, "number of processors:", numbCommProcs, &
                "size inArray:", size(inArray), "size(outArray):", size(outArray) 
           call Driver_abortFlash("Severe error: array will be accessed out of bounds(2)... exiting")
           
        end if


        call MPI_ALLTOALL(outArray,nPerProc,FLASH_REAL,inArray,nPerProc,&
                          FLASH_REAL,comm,ierr)
        !! Now the shape of outArray is (outLen(JAXIS),outLen(KAXIS),
        !! outLen(KAXIS))
     end if
     nx=outLen(JAXIS)*outLen(KAXIS)
     nx1=outLen(IAXIS)
     call gr_pfftLocalTranspose(inArray,outArray,nx,nx1,baseDatType)
  end if
  return
end subroutine gr_pfftTransposeImplementation
