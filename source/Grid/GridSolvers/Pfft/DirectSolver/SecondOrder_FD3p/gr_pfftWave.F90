!!****if* source/Grid/GridSolvers/Pfft/DirectSolver/SecondOrder_FD3p/gr_pfftWave
!!
!! NAME 
!!
!!   gr_pfftWave
!!
!! SYNOPSIS
!!
!!   gr_pfftWave()
!!
!! DESCRIPTION 
!!
!!  Calculate the right wave numbers for derivatives 
!! 
!! ARGUMENTS
!!
!!
!!***
subroutine gr_pfftWave()
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"

  use gr_pfftData, ONLY : pfft_wave, pfft_localLimits, pfft_transformType,&
       pfft_dimOrder, pfft_globalLen
  use Grid_data,ONLY : gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax
 

  implicit none

  real :: factor,tmp
  integer :: i,j,n, nq,beg,fin,l
  integer :: n2mh
  real :: Lx, Ly, Lz, dx, dy, dz, dx2q, dy2q, dz2q

  Lx = gr_imax - gr_imin
  Ly = gr_jmax - gr_jmin   
  Lz = gr_kmax - gr_kmin

  dx = Lx / real(pfft_globalLen(IAXIS))
  dy = Ly / real(pfft_globalLen(JAXIS))
  dz = Lz / real(pfft_globalLen(KAXIS))

  !...............................................................
  !     Modified wave numbers
  !...............................................................

  ! X wavenumbers:
  dx2q = 1.0/(dx*dx);
  n2mh = pfft_globalLen(IAXIS)/2;            
  do l=1 , n2mh
     pfft_wave(l+1,IAXIS) = cos(2.*PI*real(l)/real(pfft_globalLen(IAXIS)));
     pfft_wave(l+1,IAXIS) = 2.*(1.-pfft_wave(l+1,IAXIS))*dx2q;
  enddo
  do l=n2mh+1 , pfft_globalLen(IAXIS)-1;
     pfft_wave(l+1,IAXIS) = pfft_wave(pfft_globalLen(IAXIS)+1-l,IAXIS);
  enddo


#if NDIM >= 2
  ! Y wavenumbers:
  dy2q = 1.0/(dy*dy);
  n2mh = pfft_globalLen(JAXIS)/2;            
  do l=1 , n2mh
     pfft_wave(l+1,JAXIS) = cos(2.*PI*real(l)/real(pfft_globalLen(JAXIS)));
     pfft_wave(l+1,JAXIS) = 2.*(1.-pfft_wave(l+1,JAXIS))*dy2q;
  enddo
  do l=n2mh+1 , pfft_globalLen(JAXIS)-1
     pfft_wave(l+1,JAXIS) = pfft_wave(pfft_globalLen(JAXIS)+1-l,JAXIS);
  enddo
#endif

#if NDIM == 3
  ! Z wavenumbers:
  dz2q = 1.0/(dz*dz);
  n2mh = pfft_globalLen(KAXIS)/2;     
  do l=1 , n2mh
     pfft_wave(l+1,KAXIS) = cos(2.*PI*real(l)/real(pfft_globalLen(KAXIS)));
     pfft_wave(l+1,KAXIS) = 2.*(1.-pfft_wave(l+1,KAXIS))*dz2q;
  enddo
  do l=n2mh+1 , pfft_globalLen(KAXIS)-1
     pfft_wave(l+1,KAXIS) = pfft_wave(pfft_globalLen(KAXIS)+1-l,KAXIS);
  enddo
#endif

  return
end subroutine gr_pfftWave
