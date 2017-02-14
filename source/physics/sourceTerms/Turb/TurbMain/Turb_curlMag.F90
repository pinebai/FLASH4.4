!!****if* source/physics/sourceTerms/Turb/TurbMain/Turb_curlMag
!!
!! NAME
!!
!!  Turb_curlMag
!!
!! SYNOPSIS
!!
!!  call Turb_curlMag(real, dimension(:,:,:)(out) :: curl,
!!                    real, dimension(:,:,:)(in) :: velx,
!!                    real, dimension(:,:,:)(in) :: vely,
!!                    real, dimension(:,:,:)(in) :: velz,
!!                    integer(in) :: h,
!!                    integer(in) :: bid)
!!
!! DESCRIPTION
!!
!!  Aaron Jackson 2010
!!  Calculate the curl of the velocity field.
!!  The block id (bid) is passed in so that we can retrieve
!!  coordinate info. The curl is only computed for the 
!!  interior cells, although the indices of curl also run
!!  over the guard cells to simplify indexing.
!!  For this implementation, we only need the magnitude of
!!  the curl, so we store this information as a scalar
!!  instead of a vector in curl. h is the step size used
!!  to calculate the curl
!!
!! ARGUMENTS
!!
!!   curl : clalculated curl 
!!
!!   velx : x velocity 
!!
!!   vely : y velocity
!!
!!   velz : z velocity
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
subroutine Turb_curlMag(curl, velX, velY, velZ, h, bid)

  use Grid_interface, only : Grid_getGeometry, Grid_getBlkIndexLimits, &
                             Grid_getDeltas, Grid_getCellCoords
  use Driver_interface, only : Driver_abortFlash
  implicit none
  real, dimension(:,:,:), intent(out) :: curl
  real, dimension(:,:,:), intent(in) :: velX, velY, velZ
  integer, intent(in) :: bid, h

  integer :: geom
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(MDIM) :: deltas

  real :: inv_12_dx, inv_12_dy, inv_12_dz
  real :: inv_12_dr, inv_12_dtheta, inv_12_dphi
  real :: dvXdY, dvXdZ, dvYdX, dvYdZ, dvZdX, dvZdY
  real, dimension(:), allocatable :: inv_r, ctan, inv_sin
  integer :: i,j,k, istat
  integer :: h2
  logical :: three_dim

  if (NDIM == 1) then
     ! curl is always zero
     curl(:,:,:) = 0.e0
     return
  endif

  h2 = 2*h

  if (h2 > NGUARD) call Driver_abortFlash("Step size in Turb_curlMag is too large for the number of guard cells")

  three_dim = (NDIM == 3)

  call Grid_getGeometry(geom)
  call Grid_getBlkIndexLimits(bid, blkLimits, blkLimitsGC)
  call Grid_getDeltas(bid, deltas)
  deltas(IAXIS:KAXIS) = h * deltas(IAXIS:KAXIS)

  select case (geom)
  case (CARTESIAN)

     ! seems like any self-respecting optimizing compiler would be able to
     ! pull these out of the loop but we won't trust that
     inv_12_dx = 1.0/12.0/deltas(IAXIS)
     inv_12_dy = 1.0/12.0/deltas(JAXIS)
     if (three_dim) inv_12_dz = 1.0/12.0/deltas(KAXIS)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           
              dvXdY = ( velX(i,j-h2,k) - 8*velX(i,j-h,k)   &
                          + 8*velX(i,j+h,k) - velX(i,j+h2,k) ) * inv_12_dy
              dvYdX = ( velY(i-h2,j,k) - 8*velY(i-h,j,k)   &
                          + 8*velY(i+h,j,k) - velY(i+h2,j,k) ) * inv_12_dx
              if (three_dim) then 
                 dvXdZ = ( velX(i,j,k-h2) - 8*velX(i,j,k-h)   &
                          + 8*velX(i,j,k+h) - velX(i,j,k+h2) ) * inv_12_dz
                 dvYdZ = ( velY(i,j,k-h2) - 8*velY(i,j,k-h)   &
                          + 8*velY(i,j,k+h) - velY(i,j,k+h2) ) * inv_12_dz

                 dvZdX = ( velZ(i-h2,j,k) - 8*velZ(i-h,j,k)   &
                          + 8*velZ(i+h,j,k) - velZ(i+h2,j,k) ) * inv_12_dx
                 dvZdY = ( velZ(i,j-h2,k) - 8*velZ(i,j-h,k)   &
                          + 8*velZ(i,j+h,k) - velZ(i,j+h2,k) ) * inv_12_dy

                 curl(i,j,k) = sqrt( (dvZdY - dvYdZ)**2  +  &
                                     (dvXdZ - dvZdX)**2  +  &
                                     (dvYdX - dvXdY)**2 )
              else
                 curl(i,j,k) = dvYdX - dvXdY
              endif

           enddo
        enddo
     enddo

  case (CYLINDRICAL)

     inv_12_dr  = 1.0/12.0/deltas(IAXIS)
     inv_12_dz = 1.0/12.0/deltas(JAXIS)

     if (three_dim) then
        inv_12_dtheta = 1.0/12.0/deltas(KAXIS)

        allocate(inv_r(blkLimitsGC(HIGH,IAXIS)),STAT=istat)
        if (istat/=0) call Driver_abortFlash("Unable to allocate inv_r in Turb_fsCurl")
        call Grid_getCellCoords(IAXIS, bid, CENTER, .true., inv_r, blkLimitsGC(HIGH,IAXIS))
        inv_r(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)) = &
                 1.0/inv_r(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS))
     endif

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              dvXdY = ( velX(i,j-h2,k) - 8*velX(i,j-h,k)   &
                          + 8*velX(i,j+h,k) - velX(i,j+h2,k) ) * inv_12_dy
              dvYdX = ( velY(i-h2,j,k) - 8*velY(i-h,j,k)   &
                          + 8*velY(i+h,j,k) - velY(i+h2,j,k) ) * inv_12_dx

              if (three_dim) then
                 dvXdZ = ( velX(i,j-h2,k) - 8*velX(i,j-h,k)   &
                          + 8*velX(i,j+h,k) - velX(i,j+h2,k) ) * inv_12_dz
                 dvXdZ = dvXdZ * inv_r(i)
                 dvYdZ = ( velY(i-h2,j,k) - 8*velY(i-h,j,k)   &
                          + 8*velY(i+h,j,k) - velY(i+h2,j,k) ) * inv_12_dz 
                 dvYdZ = dvYdZ * inv_r(i)

                 dvZdX = ( velZ(i-h2,j,k) - 8*velZ(i-h,j,k)   &
                          + 8*velZ(i+h,j,k) - velZ(i+h2,j,k) ) * inv_12_dx
                 dvZdX = velZ(i,j,k)*inv_r(i) + dvZdX
                 dvZdY = ( velZ(i,j-h2,k) - 8*velZ(i,j-h,k)   &
                          + 8*velZ(i,j+h,k) - velZ(i,j+h2,k) ) * inv_12_dy
                 
                 curl(i,j,k) = sqrt( (dvZdY - dvYdZ)**2  +  &
                                     (dvXdZ - dvZdX)**2  +  &
                                     (dvYdX - dvXdY)**2 )
              else
                 curl(i,j,k) = dvYdX - dvXdY
              endif

           enddo
        enddo
     enddo

     if (three_dim) deallocate(inv_r)

  case (SPHERICAL)

     inv_12_dr  = 1.0/12.0/deltas(IAXIS)
     inv_12_dtheta = 1.0/12.0/deltas(JAXIS)

     allocate(inv_r(blkLimitsGC(HIGH,IAXIS)),STAT=istat)
     if (istat/=0) call Driver_abortFlash("Unable to allocate inv_r in Turb_fsCurl")
     call Grid_getCellCoords(IAXIS, bid, CENTER, .true., inv_r, blkLimitsGC(HIGH,IAXIS))
     inv_r(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)) = &
                   1.0/inv_r(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS))
     if (three_dim) then
        inv_12_dphi = 1.0/12.0/deltas(KAXIS)
        allocate(ctan(blkLimitsGC(HIGH,JAXIS)),STAT=istat)
        if (istat/=0) call Driver_abortFlash("Unable to allocate ctan in Turb_fsCurl")
        allocate(inv_sin(blkLimitsGC(HIGH,JAXIS)),STAT=istat)
        if (istat/=0) call Driver_abortFlash("Unable to allocate inv_sin in Turb_fsCurl")
        call Grid_getCellCoords(JAXIS, bid, CENTER, .true., ctan, blkLimitsGC(HIGH,JAXIS))
        inv_sin(blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)) = &
                  1.0/sin(ctan(blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)))
        ctan(blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)) = &
                1.0/tan(ctan(blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)))
     endif 

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              dvXdY = ( velX(i,j-h2,k) - 8*velX(i,j-h,k)   &
                          + 8*velX(i,j+h,k) - velX(i,j+h2,k) ) * inv_12_dy
              dvXdY = dvXdY * inv_r(i)
              dvYdX = ( velY(i-h2,j,k) - 8*velY(i-h,j,k)   &
                          + 8*velY(i+h,j,k) - velY(i+h2,j,k) ) * inv_12_dx
              dvYdX = velY(i,j,k) * inv_r(i) + dvYdX

              if (three_dim) then
                 dvXdZ = ( velX(i,j-h2,k) - 8*velX(i,j-h,k)   &
                          + 8*velX(i,j+h,k) - velX(i,j+h2,k) ) * inv_12_dz
                 dvXdZ = dvXdZ * inv_r(i) * inv_sin(j)
                 dvYdZ = ( velY(i-h2,j,k) - 8*velY(i-h,j,k)   &
                          + 8*velY(i+h,j,k) - velY(i+h2,j,k) ) * inv_12_dz 
                 dvYdZ = dvYdZ * inv_r(i) * inv_sin(j)

                 dvZdX = ( velZ(i-h2,j,k) - 8*velZ(i-h,j,k)   &
                          + 8*velZ(i+h,j,k) - velZ(i+h2,j,k) ) * inv_12_dx
                 dvZdX = velZ(i,j,k) * inv_r(i) + dvZdX
                 dvZdY = ( velZ(i,j-h2,k) - 8*velZ(i,j-h,k)   &
                          + 8*velZ(i,j+h,k) - velZ(i,j+h2,k) ) * inv_12_dy
                 dvZdY = ( velZ(i,j,k) * ctan(j) + dvZdY ) * inv_r(i)
                 
                 curl(i,j,k) = sqrt( (dvZdY - dvYdZ)**2  +  &
                                     (dvXdZ - dvZdX)**2  +  &
                                     (dvYdX - dvXdY)**2 )
              else
                 curl(i,j,k) = dvYdX - dvXdY
              endif

             enddo
          enddo
       enddo

       if (three_dim) then
          deallocate(inv_sin)
          deallocate(ctan)
       endif
       deallocate(inv_r)
 
  end select

  return
end subroutine Turb_curlMag
  
