!!****if* source/Grid/GridSolvers/Pfft/DirectSolver/Generic_Direct/gr_pfftDcftForward
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
!!                 PFFT_REAL  (real to complex)
!!                 PFFT_COMPLEX (complex to compex)
!!  scale - the scaling factor
!!
!! NOTES
!!
!!  Grid_pfftInit should have been called before this routine
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

  integer :: i,j,n
  real :: temp 
 
  select case(transformType)
  case(PFFT_REAL)                 ! PFFT_REAL Contiguous Complex coeffs For IN Plane Solver
     j=1
     do i=1,numVec 
        n=j+len-1
        outArray(j:n)=scale*inArray(j:n)
        call RFFTF(len,outArray(j),trig)
       j=j+lda
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
        call COSQB(len,outArray(j),trig)
        j=j+lda
     end do
  case(PFFT_SINQ)
     j=1
     do i=1,numVec
        n=j+len-1
        outArray(j:n)=scale*inArray(j:n)
        call SINQB(len,outArray(j),trig)
        j=j+lda
     end do
  end select

  return
end subroutine gr_pfftDcftForward
