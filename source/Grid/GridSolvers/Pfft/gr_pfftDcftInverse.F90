!!****if* source/Grid/GridSolvers/Pfft/gr_pfftDcftInverse
!!
!! NAME
!!
!!  gr_pfftDcftInverse
!!
!! SYNOPSIS
!!  
!!  gr_pfftDcftInverse(real(IN)   :: inArray(:),
!!                     real(OUT)  :: outArray(:),
!!                     real(IN)   :: trig(:),
!!                     integer(IN):: len,
!!                     integer(IN):: lda,
!!                     integer(IN):: numVec,
!!                     integer(IN):: transformType,
!!                     real(IN)   :: scale)
!! 
!! DESCRIPTION
!!
!!   This routine calculates inverse transforms of numVec single
!!   dimensional vectors in non-parallel mode. The allowed tranforms 
!!   are listed in the arguments
!!
!! ARGUMENTS
!!  
!!  inArray - Array containing input data
!!  outArray - Array for storing transformed data
!!  trig  - this array contains precomputed trignometric tables
!!  len - length of individual vectors
!!  lda - Sometime the vectors have some trailing spaces, so the
!!        actual storage for the vectors is larger than their size
!!        lda is the actual storage allocated for each vector
!! numVec - number of vectors being transformed.
!! transformType - the allowed value for this argument are:
!!                 PFFT_REAL2C  (real to complex)
!!                 PFFT_REAL  (real to real)
!!                 PFFT_COMPLEX (complex to compex)
!!                 PFFT_COS (cosine)
!!                 PFFT_SIN (sine)
!!                 PFFT_COSQ (cos with a phase shift)
!!                 PFFT_SINQ (sin with a phasee shift)
!!  scale - the scaling factor
!!
!! NOTES
!!
!!  Grid_pfftInit should have been called before this routine
!!
!!***

subroutine gr_pfftDcftInverse(inArray,outArray,trig,len,lda,numVec,&
     transformType,scale)
#include "Pfft.h"
  implicit none
  integer, intent(IN) :: len, lda, numVec, transformType
  real, intent(IN) :: scale
  real, dimension(:),intent(IN) :: trig,inArray
  real,dimension(:), intent(OUT) :: outArray

  integer :: i,j,n,iout,jout,nout, skip
  real :: temp
  real,allocatable :: tempOut(:)

  select case(transformType)
  case(PFFT_REAL)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=inArray(j:n)
        call RFFTB(len,outArray(j),trig)
        temp = maxval(outArray(j:n))
        j=j+lda
     end do
  case(PFFT_REAL2C)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j)=inArray(j)
        outArray(j+1:n-1)=inArray(j+2:n)
        outArray(n)=0.0
        call RFFTB(len,outArray(j),trig)
        temp = maxval(outArray(j:n))
        j=j+lda
     end do
  case(PFFT_REAL2C_STUFF)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j)=inArray(j)
        outArray(j+1:n-1)=inArray(j+2:n)
        outArray(n)=inArray(j+1)
        call RFFTB(len,outArray(j),trig)
        temp = maxval(outArray(j:n))
        j=j+lda
     end do
  case(PFFT_REAL2C_EXTEND)
     j=1; jout=1
     do i=1,numVec
        n=j+len; nout=jout+len-1
        outArray(jout)=inArray(j) !skip inArray(j+1)
        outArray(jout+1:nout)=inArray(j+2:n) !drop inArray(n+1)
        call RFFTB(len,outArray(jout),trig)
        j=j+lda+2; jout=jout+lda
     end do
  case(PFFT_COMPLEX)
     j=1
     do i=1,numVec
        n=j+2*len-1
        outArray(j:n)=inArray(j:n)
        call CFFTB(len,outArray(j),trig)
        j=j+2*lda
     end do
  case(PFFT_COS)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=inArray(j:n)
        !call COST(len,outArray(j),trig)  !f77 fftpack
        call RCOST(len,outArray(j),trig)  !f90 fftpack
        j=j+lda
     end do
  case(PFFT_SIN)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=inArray(j:n)
        !call SINT(len,outArray(j),trig)  !f77 fftpack
        call RSINT(len,outArray(j),trig)  !f90 fftpack
        j=j+lda
     end do
  case(PFFT_COSQ)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=inArray(j:n)
        call COSQB(len,outArray(j),trig)
        j=j+lda
     end do
  case(PFFT_SINQ)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=inArray(j:n)
        call SINQB(len,outArray(j),trig)
        j=j+lda
     end do
  case(PFFT_COS_CC)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=inArray(j:n)
        call COSQF(len,outArray(j),trig)
        j=j+lda
     end do
  case(PFFT_SIN_CC)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=inArray(j:n)
        call SINQF(len,outArray(j),trig)
        j=j+lda
     end do
  case(PFFT_COS_IV)
     j=1
     do i=1,numVec-1
        n=j+len-1
        outArray(j:n)=inArray(j:n)
        outArray(j+len:n+len)=-inArray(n:j:-1)
        call COSQB(2*len,outArray(j),trig)
        skip=1
        do iout=j,n
           outArray(iout)=outArray(iout+skip)
           skip=skip+1
        end do
        j=j+lda
     end do
     n=j+len-1
     if(j>1 .AND. lda .GE. 2*len) then
        outArray(j-len:n-len)=inArray(j:n)
        outArray(j:n)=-inArray(n:j:-1)
        call COSQB(2*len,outArray(j-len),trig)
        skip=1
        do iout=j,n
           outArray(iout)=outArray(iout-len+skip)
           skip=skip+1
        end do
     else
        allocate(tempOut(2*len))
        tempOut(1:len)=inArray(j:n)
        tempOut(1+len:2*len)=-inArray(n:j:-1)
        call COSQB(2*len,tempOut,trig)
        outArray(j:n)=tempOut(2:2*len:2)
        deallocate(tempOut)
     end if
  case(PFFT_SIN_IV)
     j=1
!     print*,'numVec,len are',numVec,len
     do i=1,numVec-1
        n=j+len-1
!        if(i==1) print*,'Inv inA:',j,inArray(j:n)
        outArray(j:n)=inArray(j:n)
        outArray(j+len:n+len)=inArray(n:j:-1)
        call SINQB(2*len,outArray(j),trig)
        skip=0
        do iout=j,n
           outArray(iout)=outArray(iout+skip)
           skip=skip+1
        end do
        j=j+lda
     end do
     n=j+len-1
     if(j>1 .AND. lda .GE. 2*len) then
        outArray(j-len:n-len)=inArray(j:n)
        outArray(j:n)=inArray(n:j:-1)
        call SINQB(2*len,outArray(j-len),trig)
        skip=0
        do iout=j,n
           outArray(iout)=outArray(iout-len+skip)
           skip=skip+1
        end do
     else
        allocate(tempOut(2*len))
        tempOut(1:len)=inArray(j:n)
        tempOut(1+len:2*len)=inArray(n:j:-1)
        call SINQB(2*len,tempOut,trig)
        outArray(j:n)=tempOut(1:2*len-1:2)
        deallocate(tempOut)
     end if
  end select
  return
end subroutine gr_pfftDcftInverse
