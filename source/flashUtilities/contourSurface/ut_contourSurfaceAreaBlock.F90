!!****if* source/flashUtilities/contourSurface/ut_contourSurfaceAreaBlock
!!
!! NAME
!!
!!  sim_contour_surface_area
!!
!!
!! SYNOPSIS
!!
!!  ut_contourSurfaceAreaBlock(nlevels,isolevels,ctrData, blkLimits, blkIndex, areas)
!!
!!
!! DESCRIPTION
!!
!!  Extracts contour surface of a given variable using marching
!!  cubes algorithm and computs area of this surface for all blocks
!!  on this processor.
!!
!!  note this is just a wrapper routine for a marching cubes routine
!!  the vertex data and box dimensions are calculated and passed to
!!  a subroutine that actually does the marching cubes.
!!
!!  this function should not be called in a 1d simulation
!!
!!
!! ARGUMENTS
!!
!!***

subroutine ut_contourSurfaceAreaBlock(nlevels,isolevels,ctrData, blkLimits, blkIndex, areas)

use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkBoundBox
use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "Flash.h"
#include "constants.h"

  integer,                     intent(IN)  :: nlevels
  real,dimension(nlevels),     intent(IN)  :: isolevels
  real,dimension(:,:,:),       intent(IN)  :: ctrData
  integer,dimension(HIGH,MDIM),intent(IN)  :: blkLimits
  integer,                     intent(IN)  :: blkIndex
  real,dimension(nlevels),     intent(OUT) :: areas

  real, allocatable, dimension(:,:,:) :: mcData

  integer, dimension(3) :: mcdataSize

  integer :: i,j,k,ii,jj,kk, istat
  real :: datamin, datamax
  real, dimension(MDIM) :: deltas

  real :: xmin,ymin,zmin,xmax,ymax,zmax
  real, dimension(2,MDIM) :: boundBox

  call Grid_getBlkBoundBox(blkIndex, boundBox)

  xmin = boundBox(1,IAXIS)
  ymin = boundBox(1,JAXIS)
  zmin = boundBox(1,KAXIS)
  xmax = boundBox(2,IAXIS)
  ymax = boundBox(2,JAXIS)
  zmax = boundBox(2,KAXIS)
#if NDIM == 2
  zmax = zmin+xmax-xmin
#endif

  call Grid_getDeltas(blkIndex, deltas)

  if (NDIM==2) deltas(KAXIS) = deltas(IAXIS)

  ! allocate space for data to be passed to marching cubes routines
  ! need one layer of guard cells on all sides
  ! always keep three cells in z direction for faking
  mcdataSize(1) = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+3
  mcdataSize(2) = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+3
  mcdataSize(3) = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+3
  allocate(mcData(mcdataSize(1), mcdataSize(2), mcdataSize(3)),STAT=istat)
   if (istat /= 0) call Driver_abortFlash("Cannot allocate celldata in ut_contourSurfaceAreaBlock")

  datamax = -HUGE(1.0)
  datamin = HUGE(1.0)



  if (NDIM==2) then

    kk = blkLimits(LOW,KAXIS)
    do j = 1, mcdataSize(2)
       jj = blkLimits(LOW,JAXIS)-2 + j
       do i = 1, mcdataSize(1)
          ii = blkLimits(LOW,IAXIS)-2 + i
          mcData(i,j,1) = ctrData(ii,jj,kk)
          mcData(i,j,2) = mcData(i,j,1)
          mcData(i,j,3) = mcData(i,j,1)
          datamax = max(datamax,mcData(i,j,1))
          datamin = min(datamin,mcData(i,j,1))
       end do
    end do
  
  else if (NDIM==3) then

    do k = 1, mcdataSize(3)
       kk = blkLimits(LOW,KAXIS)-2 + k
       do j = 1, mcdataSize(2)
          jj = blkLimits(LOW,JAXIS)-2 + j
          do i = 1, mcdataSize(1)
             ii = blkLimits(LOW,IAXIS)-2 + i
             mcData(i,j,k) = ctrData(ii,jj,kk)
             datamax = max(datamax,mcData(i,j,k))
             datamin = min(datamin,mcData(i,j,k))
          end do
       end do
    end do

  else

     call Driver_abortFlash("Surface area measurement only works in 2 or 3 dimensions")

  endif

  do i = 1, nlevels
     if( datamin <= isolevels(i) .AND. isolevels(i) <= datamax ) then
        call interior_iso_surface_area(mcData,mcdatasize(1),mcdataSize(2),mcdataSize(3), &
                         deltas(1),deltas(2),deltas(3),isolevels(i),areas(i))
#if NDIM == 2
        areas(i) = areas(i)/(zmax-zmin)
#endif
     else
        areas(i) = 0.0
     end if
  enddo

  deallocate(mcData)

end subroutine ut_contourSurfaceAreaBlock


