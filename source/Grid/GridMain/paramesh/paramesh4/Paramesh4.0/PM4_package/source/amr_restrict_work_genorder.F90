!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine amr_restrict_work_genorder(datainw,dataoutw,iopt,order)




!------------------------------------------------------------------------
!
! This routine performs interpolation for the restriction operation on
! cell centered data stored in 'work'.  It uses a lagrange polynomial
! interpolation scheme.  Also, the interpolation stencil is automatically
! shifted to avoid using data in guardcells. CAUTION: you must realize that
! use of this routine with 'order' set to values higher than 1 MAY
! cause asymmetric interpolation and your results may loose symmetry.
!
! Data is passed in in the array 'datainw' and returned in the array
! 'dataoutw'.  The order of the interpolating polynomial is also passed
! in the variable 'order' and can take on value ranging from 1 to 5.
! The last argument 'ivar' specifies which variable in 'work' to apply
! the interpolation to.
!
!
! Written :     Kevin Olson          March 2004
!------------------------------------------------------------------------



      use paramesh_dimensions
      use physicaldata
      use workspace

      implicit none

      real, intent(in)    :: datainw(:,:,:)
      real, intent(inout) :: dataoutw(:,:,:)
      integer, intent(in) :: iopt, order

      real    :: xi, xj, www
      real,    save :: weight(5,3,-4:5)

      integer :: i,j,k
      integer :: i0, j0, k0, is, js, ks
      integer :: iw, jw, kw, iii, jjj, kkk
      integer, save :: iparmin,iparmax
      integer, save :: jparmin,jparmax
      integer, save :: kparmin,kparmax
      integer :: istart, jstart, kstart
      integer :: iend, jend, kend
      integer :: order2

      logical, save :: first = .true.

!------------------------------------

      if (first) then

      first = .false.

      do order2 = 1, 5

! left

      xi = 0.-.5
      do i = 0,order2
         weight(order2,1,i) = 1.
         xj = 0.-.5
         do j = 0,order2
            if (i .ne. j) then
               weight(order2,1,i) = & 
     &              weight(order2,1,i)*(0.-xj)/(xi-xj)
            end if
            xj = xj + 1.
         end do
         xi = xi + 1.
      end do

! middle

      istart = -int(order2/2)
      iend = istart + order2
      xi = real(istart)-.5
      do i = istart,iend
         weight(order2,2,i) = 1.
         xj = real(istart)-.5
         do j = istart,iend
            if (i .ne. j) then
               weight(order2,2,i) = & 
     &              weight(order2,2,i)*(0.-xj)/(xi-xj)
            end if
            xj = xj + 1.
         end do
         xi = xi + 1.
      end do

! right

      istart = -order2 + 1
      iend = istart + order2
      xi = real(istart)-.5
      do i = istart,iend
         weight(order2,3,i) = 1.
         xj = real(istart)-.5
         do j = istart,iend
            if (i .ne. j) then
               weight(order2,3,i) = & 
     &              weight(order2,3,i)*(0.-xj)/(xi-xj)
            end if
            xj = xj + 1.
         end do
         xi = xi + 1.
      end do

      end do

      iparmin = 1+nguard_work
      iparmax = nxb+nguard_work
      jparmin = 1+nguard_work*k2d
      jparmax = nyb+nguard_work*k2d
      kparmin = 1+nguard_work*k3d
      kparmax = nzb+nguard_work*k3d

      end if

      do k0 = kparmin,kparmax,2
      do j0 = jparmin,jparmax,2
      do i0 = iparmin,iparmax,2

        dataoutw(i0,j0,k0) = 0.

        if (ndim == 3) then
           if (k0 == kparmin) then
              kstart = 0
              kw = 1
           elseif (k0 == kparmax-1) then
              kstart = -order + 1
              kw = 3
           else
              kstart = -int(order/2) 
              kw = 2
           end if
           ks = k0+kstart
           kend = kstart + order
        else
           ks     = 1
           kstart = 1
           kend   = 1
        end if

        if (ndim >= 2) then
           if (j0 == jparmin) then
              jstart = 0
              jw = 1
           elseif (j0 == jparmax-1) then
              jstart = -order + 1
              jw = 3
           else
              jstart = -int(order/2)
              jw = 2
           end if
           js = j0+jstart
           jend = jstart + order
        else
           js     = 1
           jstart = 1
           jend   = 1
        end if

        if (i0 == iparmin) then
           istart = 0
           iw = 1
        elseif (i0 == iparmax-1) then
           istart = -order + 1
           iw = 3
        else
           istart = -int(order/2)
           iw = 2
        end if
        is = i0+istart
        iend = istart + order

        k = ks
        do kkk = kstart,kend
           j = js
           do jjj = jstart,jend
              i = is
              do iii = istart,iend

                 if (ndim == 1) then
                    www = weight(order,iw,iii)
                 elseif (ndim == 2) then
                    www = weight(order,iw,iii)* & 
     &                    weight(order,jw,jjj)
                 elseif (ndim == 3) then
                    www = weight(order,iw,iii)* & 
     &                    weight(order,jw,jjj)* & 
     &                    weight(order,kw,kkk)
                 end if

                 if (curvilinear_conserve) then
                    dataoutw(i0,j0,k0) = & 
     &                dataoutw(i0,j0,k0) + & 
     &                (datainw(i,j,k))
                 else
                    dataoutw(i0,j0,k0) = & 
     &                dataoutw(i0,j0,k0) + & 
     &                (www*datainw(i,j,k))
                 endif

                 i = i + 1
              end do
              j = j + 1
           end do
           k = k + 1
        end do

      enddo
      enddo
      enddo

      end subroutine amr_restrict_work_genorder
