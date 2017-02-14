!!****if* source/Grid/GridSolvers/Multipole/gr_mpoleSphABTerms
!!
!! NAME
!!
!!  gr_mpoleSphABTerms
!! 
!!
!! SYNOPSIS
!!
!!  gr_mpoleSphABTerms(integer, intent(in) :: blockID)
!!
!!
!! DESCRIPTION
!! 
!!  Compute intermediate terms for the given block to be used in calculation
!!  of moments. These terms are used in the recurrence relations.
!!
!! 
!! ARGUMENTS
!!
!!  blockID - my block number
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleSphABTerms (blockID)

  use gr_mpoleData, ONLY : G_2DSPHERICAL,&
                         gbnd, xmin, lstep, Moment, dsinv, r2, yzn, xzn,&
                         pint, pleg, mpole_lmax

  use Grid_interface, ONLY : Grid_getCellCoords, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  
#include "constants.h"
#include "Flash.h"

  implicit none
  
  integer, intent(in) :: blockID
  integer,dimension(LOW:HIGH,MDIM):: blkLimits,blkLimitsGC
  integer   :: i, j, k,l,k2
  integer :: m, isize, jsize
  integer :: startingPos(MDIM), dataSize(MDIM)
  real,pointer,dimension(:,:,:,:):: solnData
  real :: tmp,tmp1

!------------------------------------------------------------------------

  ! interface coordinates
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  isize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jsize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
 
  call Grid_getCellCoords(IAXIS, blockID, RIGHT_EDGE, .true., xzn, isize)
  call Grid_getCellCoords(JAXIS, blockID, RIGHT_EDGE, .true., yzn, jsize)


! generate Legendre polynomials

  pleg(:,0) = 1.e0
  pleg(:,1) = cos(yzn) 
     
  do k = 2,mpole_lmax+1   !! diagnostic determine whether this is lmax or mmax
     tmp1 = 1.e0 - 1.e0/float(k)   
     tmp  = 1.e0 + tmp1
     pleg(:,k) = tmp*pleg(:,1)*pleg(:,k-1) - tmp1*pleg(:,k-2)
  end do

! integrals of Legendre polynomials:
!                               int_{th_{i-1}}^{th_i} sin th d th pleg(cos th)

  do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
     pint(j,0) = pleg(j,1)-pleg(j-1,1)
  end do
  
  do k = lstep, mpole_lmax, lstep
     tmp1 = 1.e0/float(k)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        pint(j,k) = (pleg(j,k+1) - pleg(j,1) * pleg(j,k)&
             & - pleg(j-1,k+1) + pleg(j-1,1) * pleg(j-1,k))*tmp1
     end do
  end do
  
!-------------------------------------------------------
!
! to avoid differences of nearly identical terms build :
!     
! sum of the inner mass contribution from center to surface
! sum of the outer mass contribution from surface to center
!
! Description of the Moment array
!
! Moment(:,1,1,l,0), contributions on all possible position vectors
!                    for the lth moment (term atrm)
! Moment(:,2,1,l,0), same as above, btrm
! Moment(:,1,2,l,0), fac1 ((xznl/xznr)**(l+1))
! Moment(:,2,2,l,0), fac2 ((xznr/xznr+1)**l
!
!-------------------------------------------------------


  call Grid_getBlkPtr(blockID,solnData,CENTER)
  ! identify the innermost zone of global grid
  
  i = blkLimits(LOW,IAXIS)
  if ( xzn(NGUARD) == xmin ) then
     
     gbnd = 0.e0
     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        gbnd = gbnd + 0.5e0*solnData(DENS_VAR,i,j,1) *pint(i,0)
     end do

     gbnd = gbnd*(r2(1)-r2(0))

  end if

! traverse global radial grid

  k = dsinv*(xzn(i)-xzn(i-1)) - 1

  do l = 0, mpole_lmax, lstep
     tmp1 = 1.e0/float(l+3)
     
     !    i-doloop is the outer loop now
     
     do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
        k2 = dsinv*(xzn(i)-xmin)

!       all factors should be later integrated over global grid in j

!       precompute atrm(1) factor

        tmp = 0.e0

        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           tmp = tmp + solnData(DENS_VAR,i,j,1)*pint(j,l)
        end do

!       get atrm

        do m = k2,k2+k
           Moment(m,1,1,l,0) = Moment(m,1,1,l,0) &
                              +tmp1*tmp*(r2(m)-r2(m-1)*Moment(m,1,2,l,0))
        end do

!       get btrm

        if ( l == 2 ) then

           do m = k2,k2+k
              Moment(m-1,2,1,l,0) = Moment(m-1,2,1,l,0)                   &
                                   +tmp*r2(m)*log(1.e0/Moment(m,1,2,0,0)) &
                                   *Moment(m-1,2,2,l,0)
           end do

        else

           do m = k2,k2+k
              Moment(m-1,2,1,l,0) = Moment(m-1,2,1,l,0)                     &
                                   +tmp*(r2(m)*Moment(m-1,2,2,l,0)-r2(m-1)) &
                                   /float(2-l)
           end do

        end if

     end do ! i

  end do ! l
  call Grid_releaseBlkPtr(blockID,solnData)
end subroutine gr_mpoleSphABTerms



