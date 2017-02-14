!!****if* source/physics/sourceTerms/Turb/TurbMain/Turb_laplacian
!!
!! NAME
!!
!!  Turb_laplacian
!!
!! SYNOPSIS
!!
!!  call Turb_laplacian(real, dimension(:,:,:)(out) :: lapl,
!!                      real, dimension(:,:,:)(in) :: scalar,
!!                      integer(in) :: h,
!!                      integer(in) :: bid)
!!
!! DESCRIPTION
!!      
!!      calculate the laplacian of a scalar field
!!       
!!      Dean Townsley 2008
!!       
!!      Calculate the laplacian of the scalar field "scalar".
!!      The block id (bid) is passed in so that we can retrieve
!!      coordinate info. The laplacian is only computed for the
!!      interior cells, although the indices of lapl also run
!!      over the guard cells to simplify indexing. h is the step size
!!      used to calculate the laplacian
!!
!! ARGUMENTS
!!
!!   lapl : laplacian
!!
!!   scalar : scalar field
!!
!!   h : step size
!!
!!   bid : block id
!!
!!
!!
!!***


#include "Flash.h"
#include "constants.h"
subroutine Turb_laplacian(lapl, scalar, h, bid)

  use Turb_interface, only : Turb_calcCompLimits
  use Grid_interface, only : Grid_getGeometry, Grid_getBlkIndexLimits, &
                             Grid_getDeltas, Grid_getCellCoords
  use Driver_interface, only : Driver_abortFlash
  implicit none
  real, dimension(:,:,:), intent(out) :: lapl
  real, dimension(:,:,:), intent(in) :: scalar
  integer, intent(in) :: bid, h

  integer :: geom
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC, compLimits
  real, dimension(MDIM) :: deltas

  real :: inv_12_dx2, inv_12_dy2, inv_12_dz2
  real :: inv_12_dr2, inv_12_dr, inv_12_dtheta, inv_12_dtheta2, inv_12_dphi2
  real :: d2, d1
  real, dimension(:), allocatable :: inv_r, inv_r2, two_over_r, ctan, inv_sin2
  integer :: i,j,k, istat
  integer :: h2

  h2 = 2*h

  call Grid_getGeometry(geom)
  call Grid_getBlkIndexLimits(bid, blkLimits, blkLimitsGC)
  call Turb_calcCompLimits(blkLimits,compLimits, h)
  call Grid_getDeltas(bid, deltas)
  deltas(IAXIS) = h * deltas(IAXIS)
  if (NDIM >= 2) deltas(JAXIS) = h * deltas(JAXIS)
  if (NDIM == 3) deltas(KAXIS) = h * deltas(KAXIS)

  select case (geom)
  case (CARTESIAN)

     ! seems like any self-respecting optimizing compiler would be able to
     ! pull these out of the loop but we won't trust that
     inv_12_dx2 = 1.0/12.0/deltas(IAXIS)**2
     if (NDIM >= 2) inv_12_dy2 = 1.0/12.0/deltas(JAXIS)**2
     if (NDIM == 3) inv_12_dz2 = 1.0/12.0/deltas(KAXIS)**2

     do k = compLimits(LOW,KAXIS), compLimits(HIGH,KAXIS)
        do j = compLimits(LOW,JAXIS), compLimits(HIGH,JAXIS)
           do i = compLimits(LOW,IAXIS), compLimits(HIGH,IAXIS)
           
              lapl(i,j,k) = ( -scalar(i-h2,j,k) + 16*scalar(i-h,j,k) -30*scalar(i,j,k) &
                                             + 16*scalar(i+h,j,k) - scalar(i+h2,j,k) ) * inv_12_dx2
              if (NDIM >= 2) then
                 lapl(i,j,k) = lapl(i,j,k) + &
                            ( -scalar(i,j-h2,k) + 16*scalar(i,j-h,k) -30*scalar(i,j,k) &
                                             + 16*scalar(i,j+h,k) - scalar(i,j+h2,k) ) * inv_12_dy2
              endif
              if (NDIM == 3) then
                 lapl(i,j,k) = lapl(i,j,k) + &
                            ( -scalar(i,j,k-h2) + 16*scalar(i,j,k-h) -30*scalar(i,j,k) &
                                             + 16*scalar(i,j,k+h) - scalar(i,j,k+h2) ) * inv_12_dz2
              endif

           enddo
        enddo
     enddo

  case (CYLINDRICAL)

     inv_12_dr  = 1.0/12.0/deltas(IAXIS)
     inv_12_dr2 = inv_12_dr/deltas(IAXIS)

     allocate(inv_r(blkLimitsGC(HIGH,IAXIS)),STAT=istat)
     if (istat/=0) call Driver_abortFlash("Unable to allocate inv_r in Turb_laplacian")
     call Grid_getCellCoords(IAXIS, bid, CENTER, .true., inv_r, blkLimitsGC(HIGH,IAXIS))
     inv_r(compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS)) = &
              1.0/inv_r(compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS))

     if (NDIM >= 2) inv_12_dz2 = 1.0/12.0/deltas(JAXIS)**2
     if (NDIM == 3) then
        inv_12_dtheta2 = 1.0/12.0/deltas(KAXIS)**2
        allocate(inv_r2(blkLimitsGC(HIGH,IAXIS)),STAT=istat)
        if (istat/=0) call Driver_abortFlash("Unable to allocate inv_r2 in Turb_laplacian")
        inv_r2(compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS)) = &
                  inv_r(compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS))**2
     endif


     do k = compLimits(LOW,KAXIS), compLimits(HIGH,KAXIS)
        do j = compLimits(LOW,JAXIS), compLimits(HIGH,JAXIS)
           do i = compLimits(LOW,IAXIS), compLimits(HIGH,IAXIS)

              d2 = ( -scalar(i-h2,j,k) + 16*scalar(i-h,j,k) -30*scalar(i,j,k) &
                                    + 16*scalar(i+h,j,k) - scalar(i+h2,j,k) ) * inv_12_dr2
              d1 = ( scalar(i-h2,j,k) -8*scalar(i-h,j,k)  +8*scalar(i+h,j,k) -scalar(i+h2,j,k) )*inv_12_dr

              lapl(i,j,k) = d2 + inv_r(i)*d1
                            
              if (NDIM >= 2) then
                 d2 = ( -scalar(i,j-h2,k) + 16*scalar(i,j-h,k) -30*scalar(i,j,k) &
                                       + 16*scalar(i,j+h,k) - scalar(i,j+h2,k) ) * inv_12_dz2
                 lapl(i,j,k) = lapl(i,j,k) +  d2
              endif
              if (NDIM == 3) then
                 d2 = ( -scalar(i,j,k-h2) + 16*scalar(i,j,k-h) -30*scalar(i,j,k) &
                                       + 16*scalar(i,j,k+h) - scalar(i,j,k+h2) ) * inv_12_dtheta2
                 lapl(i,j,k) = lapl(i,j,k) +  inv_r2(i)*d2
              endif

           enddo
        enddo
     enddo

     if (NDIM==3) deallocate(inv_r2)
     deallocate(inv_r)

  case (SPHERICAL)

     inv_12_dr  = 1.0/12.0/deltas(IAXIS)
     inv_12_dr2 = inv_12_dr/deltas(IAXIS)

     allocate(inv_r2(blkLimitsGC(HIGH,IAXIS)),STAT=istat)
     if (istat/=0) call Driver_abortFlash("Unable to allocate inv_r2 in Turb_laplacian")
     allocate(two_over_r(blkLimitsGC(HIGH,IAXIS)),STAT=istat)
     if (istat/=0) call Driver_abortFlash("Unable to allocate two_over_r in Turb_laplacian")
     call Grid_getCellCoords(IAXIS, bid, CENTER, .true., inv_r2, blkLimitsGC(HIGH,IAXIS))
     two_over_r(compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS)) = &
                   2.0/inv_r2(compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS))
     inv_r2(compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS)) = &
                   1.0/inv_r2(compLimits(LOW,IAXIS):compLimits(HIGH,IAXIS))**2

     if (NDIM >= 2) then
        inv_12_dtheta = 1.0/12.0/deltas(JAXIS)
        inv_12_dtheta2 = inv_12_dtheta/deltas(JAXIS)
        allocate(ctan(blkLimitsGC(HIGH,JAXIS)),STAT=istat)
        if (istat/=0) call Driver_abortFlash("Unable to allocate ctan in Turb_laplacian")
        call Grid_getCellCoords(IAXIS, bid, CENTER, .true., ctan, blkLimitsGC(HIGH,IAXIS))
        ctan(compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS)) = &
                1.0/tan(ctan(compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS)))
     endif 

     if (NDIM == 3) then
        inv_12_dphi2 = 1.0/12.0/deltas(KAXIS)**2
        allocate(inv_sin2(blkLimitsGC(HIGH,JAXIS)),STAT=istat)
        if (istat/=0) call Driver_abortFlash("Unable to allocate inv_sin2 in Turb_laplacian")
        inv_sin2(compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS)) = &
                    1.0+ctan(compLimits(LOW,JAXIS):compLimits(HIGH,JAXIS))**2
     endif


     do k = compLimits(LOW,KAXIS), compLimits(HIGH,KAXIS)
        do j = compLimits(LOW,JAXIS), compLimits(HIGH,JAXIS)
           do i = compLimits(LOW,IAXIS), compLimits(HIGH,IAXIS)
                d2 = ( -scalar(i-h2,j,k) + 16*scalar(i-h,j,k) -30*scalar(i,j,k) &
                                      + 16*scalar(i+h,j,k) - scalar(i+h2,j,k) ) * inv_12_dr2
                d1 = ( scalar(i-h2,j,k)-8*scalar(i-h,j,k)+8*scalar(i+h,j,k) -scalar(i+h2,j,k) )*inv_12_dr
                lapl(i,j,k) = d2 + two_over_r(i)*d1

                if (NDIM>=2) then
                   d2 = ( -scalar(i,j-h2,k) + 16*scalar(i,j-h,k) -30*scalar(i,j,k) &
                                         + 16*scalar(i,j+h,k) - scalar(i,j+h2,k) ) * inv_12_dtheta2
                   d1 = ( scalar(i,j-h2,k)-8*scalar(i,j-h,k)+8*scalar(i,j+h,k)-scalar(i,j+h2,k) )*inv_12_dtheta
                   lapl(i,j,k) = lapl(i,j,k) + inv_r2(i)*( d2 + ctan(j)*d1 )
                endif
                if (NDIM==3) then
                   d2 = ( -scalar(i,j,k-h2) + 16*scalar(i,j,k-h) -30*scalar(i,j,k) &
                                         + 16*scalar(i,j,k+h) - scalar(i,j,k+h2) ) * inv_12_dphi2
                   lapl(i,j,k) = lapl(i,j,k) + inv_r2(i)*inv_sin2(j)*d2
                endif
             enddo
          enddo
       enddo

       if (NDIM==3) deallocate(inv_sin2)
       if (NDIM>=2) deallocate(ctan)
       deallocate(inv_r2)
       deallocate(two_over_r)

 
  end select

  return
end subroutine Turb_laplacian
