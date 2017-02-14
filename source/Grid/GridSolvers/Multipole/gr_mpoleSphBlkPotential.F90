!!****if* source/Grid/GridSolvers/Multipole/gr_mpoleSphBlkPotential
!!
!!  NAME
!!
!!    gr_mpoleSphBlkPotential
!! 
!!
!! SYNOPSIS
!!
!!   gr_mpoleSphBlkPotential(integer, intent(in) :: blockID)
!!
!! DESCRIPTION
!! 
!! Compute gravitational potential for the given block.
!!
!! ARGUMENTS
!!
!!   blockID - my block number
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleSphBlkPotential(blockID)
  
  use gr_mpoleData, ONLY : point_mass, point_mass_rsoft, gpot, twopi,newton,&
                         gbnd, xmin, lstep, Moment, dsinv, yzn, xzn,&
                         pint, pleg, mpole_lmax

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_getCellCoords

#include "constants.h"
#include "Flash.h"

  implicit none
  
  integer,intent(IN) :: blockID
  
  integer   :: i, j, k, l, k2, i_i
  integer :: isize,jsize
  real :: tmp, tmp1, gbnd_, gfac
  real,pointer,dimension(:,:,:,:) :: solnData
  integer,dimension(LOW:HIGH,MDIM)::blkLimits, blkLimitsGC
!------------------------------------------------------------------------------

! symmetry-dependent constant

  gfac = float(lstep)*twopi*Newton

! interface coordinates
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  isize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jsize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
 
  call Grid_getCellCoords(IAXIS, blockID, RIGHT_EDGE, .true., xzn, isize)
  call Grid_getCellCoords(JAXIS, blockID, RIGHT_EDGE, .true., yzn, jsize)


! generate Legendre polynomials

  pleg(:,0) = 1.e0
  pleg(:,1) = cos(yzn) 

  do k = 2,mpole_lmax+1   ! diagnostic determine whether this is lmax or mmax
     tmp1 = 1.e0 - 1.e0/float(k)   
     tmp  = 1.e0 + tmp1
     pleg(:,k) = tmp*pleg(:,1)*pleg(:,k-1) - tmp1*pleg(:,k-2)
  enddo

! integrals of Legendre polynomials:
!                          int_{th_{i-1}}^{th_i} sin th d th pleg(cos th)
  
  do j = blkLimits(LOW,JAXIS)-1,blkLimits(HIGH,JAXIS)+1
     pint(j,0) = pleg(j,1) - pleg(j-1,1)
  end do
  
  do k = lstep, mpole_lmax, lstep  !! ! diagnostic whether mmax/lmax
     tmp1 = 1.e0/float(k)
     do j = blkLimits(LOW,JAXIS)-1,blkLimits(HIGH,JAXIS)+1
        pint(j,k) = (pleg(j,k+1)-pleg(j,k)*pleg(j,k)                &
                   - pleg(j-1,k+1)+pleg(j-1,k)*pleg(j-1,k))*tmp1
     end do
  end do

! initialize potential

  if ( point_mass.ne.0.e0 ) then

     tmp1 = point_mass_rsoft*(xzn(NGUARD+1)-xzn(NGUARD))

     do i = blkLimits(LOW,IAXIS)-1,blkLimits(HIGH,IAXIS)
        if ( abs(xzn(i)) < tmp1 ) then
           tmp = tmp1
        else
           tmp = xzn(i)
        end if

        gpot(i,:,1) = -Newton * point_mass / tmp
     end do

   else

      gpot(:,:,1) = 0.e0

   end if

! all interfaces except for the innermost one

  if ( xzn(NGUARD) == xmin ) then
     i_i = NGUARD+1
  else
     i_i = NGUARD
  end if

  do i = i_i,blkLimits(HIGH,IAXIS)
     k2 = dsinv*(xzn(i)-xmin)

     do l = 0,mpole_lmax
        do j = blkLimits(LOW,JAXIS)-1,blkLimits(HIGH,JAXIS)
           gpot(i,j,1) = gpot(i,j,1)                               &
                        +gfac*pleg(j,l)                            &
                        *(Moment(k2,1,1,l,0)+Moment(k2,2,1,l,0))
        end do
     end do
  end do

! deal with the innermost interface

  if ( xzn(NGUARD) == xmin ) then

     gbnd_ = (Moment(1,2,1,0,0)+gbnd)*gfac

     do j = blkLimits(LOW,JAXIS)-1,blkLimits(HIGH,JAXIS)
        gpot(NGUARD,j,1) = gpot(NGUARD,j,1) + gbnd_
     end do

  end if

! interpolate from the interfaces to the zone centers
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
     do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
        solnData(GPOT_VAR,i,j,1) = 0.25e0                             &
             *( gpot(i-1,j-1,1)+gpot(i,j-1,1)    &
             +gpot(i-1,j  ,1)+gpot(i,j  ,1))
     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData)

  return
end subroutine gr_mpoleSphBlkPotential
  
