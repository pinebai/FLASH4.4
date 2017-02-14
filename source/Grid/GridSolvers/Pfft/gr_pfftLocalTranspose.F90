!!****if* source/Grid/GridSolvers/Pfft/gr_pfftLocalTranspose
!!
!! NAME
!! 
!!  gr_pfftLocalTranspose
!!
!! SYNOPSIS
!!
!!  gr_pfftLocalTranspose(real(IN)   :: inArray(:),
!!                        real(OUT)  :: outArray(:),
!!                        integer(IN):: nx,
!!                        integer(IN):: nx1,
!!                        integer(IN):: type)
!! 
!! DESCRIPTION
!!
!!   This routine is for local transpositions of 2 dimensional data
!!   the routine can handle real and complex data
!!   
!! ARGUMENTS
!!  
!!  inArray - Array containing input data
!!  outArray - Array for storing transposed data
!!  nx - length of the first dimension
!!  nx1 - length of the second dimension
!!  type - the datatype (real or complex)
!! 
!!***

subroutine gr_pfftLocalTranspose(inArray,outArray,nx,nx1,baseDatType)

#include "Pfft.h"
  implicit none
  integer,intent(IN) :: nx,nx1,baseDatType
  real,intent(IN),dimension(:) :: inArray
  real,intent(OUT), dimension(:) :: outArray
  integer :: i,j,k,ind1,ind2


  if( .NOT.((baseDatType==PFFT_PCLDATA_COMPLEX).or.&
       (baseDatType==PFFT_PCLDATA_COMPLEX_STUFFED).or.&
       (baseDatType==PFFT_PCLDATA_COMPLEX_EXTENDED)) ) then  !! If the data is real, then simple transpose
     do j = 1,nx1
        do i = 1,nx
           outArray(j+nx1*(i-1))=inArray(i+nx*(j-1))
        end do
     end do
  else   !! if data is complex then two consecutive data must stay together
     do j = 1,nx1
        do i=1,nx
           ind1=2*j-1+2*nx1*(i-1)
           ind2=2*i-1+2*nx*(j-1)
           outArray(ind1) = inArray(ind2)
           outArray(ind1+1) = inArray(ind2+1)
        end do
     end do
  end if
  return
end subroutine gr_pfftLocalTranspose
