!!****if* source/Grid/GridSolvers/Pfft/DirectSolver/SecondOrder_FD3p/gr_pfftDerivs
!!
!! NAME 
!!
!!   gr_pfftDerivs
!!
!! SYNOPSIS
!!
!!   gr_pfftDerivs(real(INOUT) :: transArray)
!!
!! DESCRIPTION 
!!
!!  Calculate derivatives of the transformed data inplace
!! 
!! ARGUMENTS
!!
!!  transArray : At input it contains transformed data, at 
!!               output it has the derivatives of the transformed data
!!
!! NOTES
!!   In Numerical Recipes in Fortran (1986) this is supposed to be equation 17.4.5
!!   In Antia, this is supposed to be equation 13.170
!!   However, the current implementation of gr_pfftWave does NOT have any cosines in it.
!!   And there are no cosines in this routine.  So LBR'll be damned if she can figure out
!!   how this routine is supposed to be solving
!!     \hat{u} = \hat{rho} / 2/delx**2 *(cos (2*pi*m/M)-1))
!!
!!***
subroutine gr_pfftDerivs(transArray)

#include "constants.h"  
#include "Pfft.h"

  use gr_pfftData, ONLY : pfft_myPE, pfft_wave, pfft_localLimits, pfft_outLen,&
       pfft_transformType,pfft_dimOrder,pfft_ndim, pfft_usableProc

  use Grid_interface, ONLY : Grid_getDeltas

  implicit none
  
  real,dimension(:),intent(INOUT) :: transArray

  integer :: jfactor, kfactor
  real,dimension(MDIM) :: delta
  integer :: blockID=1
  real :: partialJ,partialK,partial
  integer :: i,j,k,ii,jj,kk,n
  logical :: isComplex
  integer :: multiplier,jOffset,kOffset

  if (.not.pfft_usableProc) return

  multiplier=2  ! Data entering this particular implementation is always complex

  jfactor=max(1,pfft_outLen(IAXIS))*multiplier
  kfactor=max(1,jfactor*pfft_outLen(JAXIS))

!!$ print *, "pfft_myPE:", pfft_myPE, &
!!$       "pfft_localLimits(IAXIS):", pfft_localLimits(:,IAXIS), &
!!$       "pfft_localLimits(JAXIS):", pfft_localLimits(:,JAXIS), &
!!$       "pfft_localLimits(KAXIS):", pfft_localLimits(:,KAXIS)

  kk=0
  kOffset=0
  do k = pfft_localLimits(LOW,KAXIS),pfft_localLimits(HIGH,KAXIS)
     if(pfft_ndim==MDIM) then
        kk=kk+1
        partialK=pfft_wave(k+1,JAXIS) ! Warning J
     else
        partialK=0.0
     end if
     jj=0
     jOffset=0
     do j = pfft_localLimits(LOW,JAXIS),pfft_localLimits(HIGH,JAXIS)
        if(pfft_ndim>1) then
           jj=jj+1
           partialJ=pfft_wave(j+1,IAXIS)+partialK  ! Warning I
        else
           partialJ=0.0
        end if
        ii=0
        n=kOffset+jOffset+1
        do i= pfft_localLimits(LOW,IAXIS),pfft_localLimits(HIGH,IAXIS)
           ii=ii+1
           partial=pfft_wave(i+1,KAXIS)+partialJ  ! Warning K
           if(partial/=0.0)partial= -1.0/partial
          
              transArray(n:n+1)=transArray(n:n+1)*partial
              if (i+j+k .eq. 0) transArray(n:n+1)= 0.
              n=n+2
        end do
        jOffset=jOffset+jfactor
     end do
     kOffset=kOffset+kfactor
  end do
  
end subroutine gr_pfftDerivs
