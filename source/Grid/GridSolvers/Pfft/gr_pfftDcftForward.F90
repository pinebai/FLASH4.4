!!****if* source/Grid/GridSolvers/Pfft/gr_pfftDcftForward
!!
!! NAME
!!
!!  gr_pfftDcftForward
!!
!! SYNOPSIS
!!  
!!  gr_pfftDcftForward(real(IN)   :: inArray(:),
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
!!   This routine calculates forward transforms of numVec single
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
!!                 PFFT_COMPLEX (complex to complex)
!!                 PFFT_COS (cosine)
!!                 PFFT_SIN (sine)
!!                 PFFT_COSQ (cos with a phase shift)
!!                 PFFT_SINQ (sin with a phase shift)
!!  scale - the scaling factor
!!
!! NOTES
!!
!!  Grid_pfftInit should have been called before this routine
!!  DEV: Currently, only PFFT_REAL2C, PFFT_REAL, and PFFT_COMPLEX are somewhat tested
!!
!!***

subroutine gr_pfftDcftForward(inArray,outArray,trig,len,lda,numVec,&
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
        outArray(j:n)=scale*inArray(j:n)
        call RFFTF(len,outArray(j),trig)
       j=j+lda
     end do

  case(PFFT_REAL2C)
     j=1
     do i=1,numVec-1
        n=j+len-1
        outArray(j+1:n+1)=scale*inArray(j:n)
        call RFFTF(len,outArray(j+1),trig)
        outArray(j)=outArray(j+1)
        outArray(j+1)=0.0
       j=j+lda
     end do
     n=j+len-1
     outArray(j:n)=scale*inArray(j:n)
     call RFFTF(len,outArray(j),trig)
     do i=0,len-3
        outArray(n-i)=outArray(n-i-1)
     end do
     outArray(j+1)=0.0
  case(PFFT_REAL2C_STUFF)
     ! like PFFT_REAL2C, but stuffing Nyquist component into the
     ! imaginary part of the DC component (instead of dropping
     ! the former and zeroing the latter). - KW
     j=1
     do i=1,numVec-1
        n=j+len-1
        outArray(j+1:n+1)=scale*inArray(j:n)
        call RFFTF(len,outArray(j+1),trig)
        outArray(j)=outArray(j+1)
        outArray(j+1)=outarray(n+1)
       j=j+lda
     end do
     n=j+len-1
     outArray(j:n)=scale*inArray(j:n)
     call RFFTF(len,outArray(j),trig)
     temp = outarray(n)
     do i=0,len-3
        outArray(n-i)=outArray(n-i-1)
     end do
     outArray(j+1)=temp
  case(PFFT_REAL2C_EXTEND)
     ! like PFFT_REAL2C, but extending length of output by one
     ! complex number to store Nyquist component (with zero
     ! imaginary part) as well as DC component (also with zero
     ! imaginary part). - KW
     j=1; jout=1
     do i=1,numVec
        n=j+len-1; nout=jout+len+1
        outArray(jout+1:jout+len)=scale*inArray(j:n)
        call RFFTF(len,outArray(jout+1),trig)
        outArray(jout)=outArray(jout+1)
        outArray(jout+1)=0.0
        outarray(nout)=0.0
        j=j+lda; jout=jout+lda+2
     end do
  case(PFFT_COMPLEX)
     j=1
     do i=1,numVec
        n=j+2*len-1
        outArray(j:n)=scale*inArray(j:n)
        call CFFTF(len,outArray(j),trig)
        j=j+2*lda
     end do
  case(PFFT_COS)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=scale*inArray(j:n)
!        call COST(len,outArray(j),trig)  !f77 fftpack
        call RCOST(len,outArray(j),trig)  !f90 fftpack
        j=j+lda
     end do
  case(PFFT_SIN)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=scale*inArray(j:n)
!        call SINT(len,outArray(j),trig)  !f77 fftpack
        call RSINT(len,outArray(j),trig)  !f90 fftpack
        j=j+lda
     end do
  case(PFFT_COSQ)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=scale*inArray(j:n)
        call COSQF(len,outArray(j),trig)
        j=j+lda
     end do
  case(PFFT_SINQ)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=scale*inArray(j:n)
        call SINQF(len,outArray(j),trig)
        j=j+lda
     end do
  case(PFFT_COS_CC)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=scale*inArray(j:n)
        call COSQB(len,outArray(j),trig)
        j=j+lda
     end do
  case(PFFT_SIN_CC)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=scale*inArray(j:n)
        call SINQB(len,outArray(j),trig)
        j=j+lda
     end do
  case(PFFT_COS_IV)
     j=1
     do i=1,numVec-1
        n=j+len-1
        outArray(j:n)=scale*inArray(j:n)
        outArray(j+len:n+len)=-scale*inArray(n:j:-1)
!        if(i==1) print*,'Fwd-ouA:',j,outArray(j:n+len)
        call COSQB(2*len,outArray(j),trig)
!        if(i==1) print*,'Fwd+ouA:',j,outArray(j:n+len)
        skip=1
        do iout=j,n
           outArray(iout)=outArray(iout+skip)
           skip=skip+1
        end do
        j=j+lda
     end do
     n=j+len-1
     if (j>1 .AND. lda .GE. 2*len) then
        outArray(j-len:n-len)=scale*inArray(j:n)
        outArray(j:n)=-scale*inArray(n:j:-1)
        call COSQB(2*len,outArray(j-len),trig)
        skip=1
        do iout=j,n
           outArray(iout)=outArray(iout-len+skip)
           skip=skip+1
        end do
     else
        allocate(tempOut(2*len))
        tempOut(1:len)=scale*inArray(j:n)
        tempOut(1+len:2*len)=-scale*inArray(n:j:-1)
        call COSQB(2*len,tempOut,trig)
        outArray(j:n)=tempOut(2:2*len:2)
        deallocate(tempOut)
     end if
  case(PFFT_SIN_IV)
     j=1
     do i=1,numVec-1
        n=j+len-1
        outArray(j:n)=scale*inArray(j:n)
        outArray(j+len:n+len)=scale*inArray(n:j:-1)
!        if(i==1) print*,'Fwd-ouA:',j,outArray(j:n+len)
        call SINQB(2*len,outArray(j),trig)
!        if(i==1) print*,'Fwd+ouA:',j,outArray(j:n+len)
        skip=0
        do iout=j,n
           outArray(iout)=outArray(iout+skip)
           skip=skip+1
        end do
        j=j+lda
     end do
     n=j+len-1
     if(j>1 .AND. lda .GE. 2*len) then
        outArray(j-len:n-len)=scale*inArray(j:n)
        outArray(j:n)=scale*inArray(n:j:-1)
        call SINQB(2*len,outArray(j-len),trig)
        skip=0
        do iout=j,n
           outArray(iout)=outArray(iout-len+skip)
           skip=skip+1
        end do
     else
        allocate(tempOut(2*len))
        tempOut(1:len)=scale*inArray(j:n)
        tempOut(1+len:2*len)=scale*inArray(n:j:-1)
        call SINQB(2*len,tempOut,trig)
        outArray(j:n)=tempOut(1:2*len-1:2)
        deallocate(tempOut)
     end if

  end select

  return
end subroutine gr_pfftDcftForward
