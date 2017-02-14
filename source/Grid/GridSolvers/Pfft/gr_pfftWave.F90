!!****if* source/Grid/GridSolvers/Pfft/gr_pfftWave
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
       pfft_dimOrder, pfft_globalLen, gr_pfftDiffOpDiscretize
  use Grid_data,ONLY : gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax

  implicit none

  real :: factor,tmp
  integer :: i,j,jj,n, nq,beg,fin
  real :: Lx, Ly, Lz, domainPhysLen(MDIM), domainPhysLenOrig(MDIM)

  Lx = gr_imax - gr_imin
  Ly = gr_jmax - gr_jmin   
  Lz = gr_kmax - gr_kmin
  domainPhysLenOrig = (/Lx, Ly, Lz/)
  domainPhysLen(1:NDIM) = domainPhysLenOrig(pfft_dimOrder(1:NDIM))

  if (gr_pfftDiffOpDiscretize==1) then
    do i = 1,NDIM
!     print *, 'gr_pfftWave, i=',i,' transformType=',pfft_transformType(pfft_dimOrder(i))
     beg=pfft_localLimits(LOW,i)
     fin=pfft_localLimits(HIGH,i)
     if(fin.GE.beg) then
!        print*,'gr_pfftWave 1:',i,pfft_transformType(pfft_dimOrder(i))
        select case(pfft_transformType(pfft_dimOrder(i)))
        case(PFFT_SIN_CC)
           factor=PI/pfft_globalLen(pfft_dimOrder(i))
           do j = beg,fin
              tmp=(j+1)*factor
              pfft_wave(j-beg+1,i)=gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
           end do
        case(PFFT_COS_CC)
           factor=PI/pfft_globalLen(pfft_dimOrder(i))
           do j = beg,fin
              tmp=j*factor
              pfft_wave(j-beg+1,i)=gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
           end do
        case(PFFT_SIN_IV)
           factor=PI/(pfft_globalLen(pfft_dimOrder(i))*2.0)
           do j = beg,fin
              tmp=(j+j+1)*factor      !DEV: ??
              pfft_wave(j-beg+1,i)=gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
           end do
        case(PFFT_COS_IV)
           factor=PI/(pfft_globalLen(pfft_dimOrder(i))*2.0)
           do j = beg,fin
              tmp=(j+j+1)*factor      !DEV: ??
              pfft_wave(j-beg+1,i)=gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
           end do
        case(PFFT_COSQ)
           factor=PI/(pfft_globalLen(pfft_dimOrder(i))*2.0)
           do j = beg,fin
              tmp=sin((j+1)*factor)*(pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              pfft_wave(j-beg+1,i)=4.0*gr_twice(tmp)
           end do
        case(PFFT_COS)
           factor=PI/(2.0*(pfft_globalLen(pfft_dimOrder(i))-1))
           do j = beg,fin
              tmp=sin(j*factor)*(pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              pfft_wave(j-beg+1,i)=4.0*gr_twice(tmp)
           end do
        case(PFFT_SINQ)
           factor=PI/(pfft_globalLen(pfft_dimOrder(i))*2.0)
           do j = beg,fin
              tmp=sin((j+1)*factor)*(pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              pfft_wave(j-beg+1,i)=4.0*gr_twice(tmp)
           end do
        case(PFFT_SIN)
           factor=PI/(2.0*(pfft_globalLen(pfft_dimOrder(i))+1))
           do j = beg,fin
              tmp=sin((j+1)*factor)*(pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              pfft_wave(j-beg+1,i)=4.0*gr_twice(tmp)
           end do
        case(PFFT_REAL2C,PFFT_REAL2C_STUFF,PFFT_REAL2C_EXTEND)
           factor = 2.0*PI / domainPhysLen(i)
           do j = beg,fin
              n=j-beg+1
              tmp=j * factor
              pfft_wave(n,i)= gr_twice(tmp)
              n=n+1
           end do
        case(PFFT_REAL)         !tested - KW
           factor=2.0*PI/pfft_globalLen(pfft_dimOrder(i))
           do j = beg,fin
              n=j-beg+1
              jj = (j+1)/2
              tmp=jj * factor
              pfft_wave(n,i)= gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              n=n+1
           end do
        case(PFFT_COMPLEX)
           factor = 2.0*PI / domainPhysLen(i)
           nq=(pfft_globalLen(pfft_dimOrder(i))+1)/2
           do j = beg,fin
              n=j-beg+1
              if(j>nq)then
                 tmp=j-pfft_globalLen(pfft_dimOrder(i))
                 tmp = tmp * factor
                 pfft_wave(n,i)= gr_twice(tmp)
              elseif(j==nq)then
                 pfft_wave(n,i)= 0.0
              else
                 tmp=j * factor
                 pfft_wave(n,i)= gr_twice(tmp)
              end if
              n=n+1
           end do
        end select
     end if
    end do
  else if (gr_pfftDiffOpDiscretize==2) then

    do i = 1,NDIM
     !print *, 'gr_pfftWave, i=',i,'means phys.',pfft_dimOrder(i),', transformType=',pfft_transformType(pfft_dimOrder(i))
     beg=pfft_localLimits(LOW,i)
     fin=pfft_localLimits(HIGH,i)
     if(fin.GE.beg) then
!        print*,'gr_pfftWave 2:',i,pfft_transformType(pfft_dimOrder(i))
        select case(pfft_transformType(pfft_dimOrder(i)))
        case(PFFT_SIN_CC)
           factor=PI/pfft_globalLen(pfft_dimOrder(i))
!!$           print*,'tt(_CC),i,beg,fin:',pfft_transformType(pfft_dimOrder(i)),i,beg,fin
           do j = beg,fin
              tmp=(j+1)*factor
              pfft_wave(j-beg+1,i)=gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
           end do
        case(PFFT_COS_CC)
           factor=PI/pfft_globalLen(pfft_dimOrder(i))
!!$           print*,'tt(_CC),i,beg,fin:',pfft_transformType(pfft_dimOrder(i)),i,beg,fin
           do j = beg,fin
              tmp=j*factor
              pfft_wave(j-beg+1,i)=gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
           end do
        case(PFFT_SIN_IV)
           factor=PI/(pfft_globalLen(pfft_dimOrder(i))*2.0)
!!$           print*,'tt(_IV),i,beg,fin:',pfft_transformType(pfft_dimOrder(i)),i,beg,fin
           do j = beg,fin
              tmp=(j+j+1)*factor      !DEV: ??
              pfft_wave(j-beg+1,i)=gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
           end do
        case(PFFT_COS_IV)
           factor=PI/(pfft_globalLen(pfft_dimOrder(i))*2.0)
!!$           print*,'tt(_CC),i,beg,fin:',pfft_transformType(pfft_dimOrder(i)),i,beg,fin
           do j = beg,fin
              tmp=(j+j+1)*factor      !DEV: ??
              pfft_wave(j-beg+1,i)=gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
           end do
        case(PFFT_COSQ)
           factor=PI/(pfft_globalLen(pfft_dimOrder(i))*2.0)
           do j = beg,fin
              tmp=sin((j+1)*factor)*(pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              pfft_wave(j-beg+1,i)=4.0*tmp*tmp
           end do
        case(PFFT_COS)
           factor=PI/(2.0*(pfft_globalLen(pfft_dimOrder(i))-1))
           do j = beg,fin
              tmp=sin(j*factor)*(pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              pfft_wave(j-beg+1,i)=4.0*tmp*tmp
           end do
        case(PFFT_SINQ)
           factor=PI/(pfft_globalLen(pfft_dimOrder(i))*2.0)
           do j = beg,fin
              tmp=sin((j+1)*factor)*(pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              pfft_wave(j-beg+1,i)=4.0*tmp*tmp
           end do
        case(PFFT_SIN)
           factor=PI/(2.0*(pfft_globalLen(pfft_dimOrder(i))+1))
           do j = beg,fin
              tmp=sin((j+1)*factor)*(pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              pfft_wave(j-beg+1,i)=4.0*tmp*tmp
           end do
        case(PFFT_REAL2C,PFFT_REAL2C_STUFF,PFFT_REAL2C_EXTEND)
           factor = 2.0*PI / pfft_globalLen(pfft_dimOrder(i))
           do j = beg,fin
              n=j-beg+1
              tmp=j * factor
              pfft_wave(n,i)= gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              n=n+1
           end do
        case(PFFT_REAL)         !tested - KW
           factor=2.0*PI/pfft_globalLen(pfft_dimOrder(i))
!!$           print*,'tt(REAL),i,beg,fin:',pfft_transformType(pfft_dimOrder(i)),i,beg,fin
           do j = beg,fin
              n=j-beg+1
              jj = (j+1)/2
              tmp=jj * factor
              pfft_wave(n,i)= gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              n=n+1
           end do
        case(PFFT_COMPLEX)
           factor = 2.0*PI / pfft_globalLen(pfft_dimOrder(i))
           nq=(pfft_globalLen(pfft_dimOrder(i))+1)/2
           do j = beg,fin
              n=j-beg+1
              if(j>nq)then
                 tmp=j-pfft_globalLen(pfft_dimOrder(i))
                 tmp = tmp * factor
                 pfft_wave(n,i)= gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              elseif(j==nq)then
                 pfft_wave(n,i)= 0.0
              else
                 tmp=j * factor
                 pfft_wave(n,i)= gr_twiceScaled(tmp, pfft_globalLen(pfft_dimOrder(i))/domainPhysLen(i))
              end if
              n=n+1
           end do
        end select
     end if
    end do

  end if

  return

contains

  real function gr_once(waveno)
    real, intent(in) :: waveno
    select case (gr_pfftDiffOpDiscretize)
    case(1)
       gr_once = waveno
    case(2)
       gr_once = sin(0.5*waveno)*2
!!$    case default
!!$       gr_once = waveno
    end select
  end function gr_once

  real function gr_twice(waveno)
    real, intent(in) :: waveno
    select case (gr_pfftDiffOpDiscretize)
    case(1)
       gr_twice = waveno*waveno
    case(2)
       gr_twice = gr_once(waveno)
       gr_twice = gr_twice*gr_twice
!!$    case default
!!$       gr_twice = waveno*waveno
    end select
  end function gr_twice

  real function gr_twiceScaled(waveno,scale)
    real, intent(in) :: waveno, scale
    select case (gr_pfftDiffOpDiscretize)
    case(1)
       gr_twiceScaled = (waveno*scale)**2
    case(2)
       gr_twiceScaled = gr_once(waveno)*scale
       gr_twiceScaled = gr_twiceScaled*gr_twiceScaled
!!$    case default
!!$       gr_twiceScaled = gr_twice(waveno)*scale*scale
    end select
  end function gr_twiceScaled
end subroutine gr_pfftWave
