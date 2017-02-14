!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_bhGenEwaldFieldLevel
!!
!! NAME
!!
!!  grv_bhGenEwaldFieldLevel
!!
!!
!! SYNOPSIS
!!
!!   call grv_bhGenEwaldFieldLevel(
!!                      integer(in) :: nx,
!!                      integer(in) :: ny,
!!                      integer(in) :: nz,
!!                      real(inout) :: field_Ewald(0:,0:,0:,0:)
!!        )
!!
!! DESCRIPTION
!!
!!   Generates the Ewald field and its partial derivatives for given axis 
!!   orientation and symmetry. Called by grv_bhGenerateEwaldField.
!!
!!
!!
!! ARGUMENTS
!!
!!   nx : number of points in ewald field in direction x
!!   ny : number of points in ewald field in direction y
!!   nz : number of points in ewald field in direction z
!!   field_Ewald - reorganised ewald field array.
!!
!!***

#include "Flash.h"

subroutine grv_bhGenEwaldFieldLevel(nx, ny, nz, field_Ewald)

  use Logfile_interface, ONLY : Logfile_stamp
  use Gravity_data, ONLY : grv_bhEwaldSeriesN, &
    grv_meshNumProcs, grv_meshMe, grv_meshComm, &
    grv_bhLx, grv_bhLy, grv_bhLz, grv_bhEwald_periodicity, &
    grv_bhPotConst, grv_bhDxI
  use grv_bhInterface, ONLY : grv_IntSimpson, grv_Coef1P2I, grv_Coef1P2I_der, &
    grv_bhAccShort, grv_bhAccLong3P, grv_bhAccLong1P2I, &
    grv_bhAccLong2P1I

  implicit none
!#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  integer, intent(IN) :: nx, ny, nz
  real, intent(OUT) :: field_Ewald(0:,0:,0:,0:)

  real, parameter :: pi = PI
  integer :: chunk, i1d, i, j, k, ni, nj, nk
  integer :: hi, hj, hk, ierr
  integer :: es_nrx, es_nry, es_nrz, es_nfx, es_nfy, es_nfz, es_radius2
  real :: ewald_alpha, ewald_dzeta, ewald_eta
  real :: ratio_p, ratio_pinv, ratio_p1, ratio_p2, ratio_pinv1, ratio_pinv2
  real :: Linv
  real :: x, y, z, xni, yni, zni, rni, rni2
  real :: cr1, cr2, cr3, cf1, cf2, cf3

  real :: field_EwaldRaw(0:12, 0:nx, 0:ny, 0:nz)
  real :: field_EwaldRawLoc(0:12, 0:nx, 0:ny, 0:nz)
  real :: ewald_row(0:12), ewald_rowPrep(0:12)

  real :: grv_bhDX, ew_rcmax

! prepare constants for generating Ewald field
! default values regardless of problem orientation
! range of coefficients in sums in real and Fourier space
  es_nrx = grv_bhEwaldSeriesN
  es_nry = grv_bhEwaldSeriesN
  es_nrz = grv_bhEwaldSeriesN
  es_nfx = grv_bhEwaldSeriesN
  es_nfy = grv_bhEwaldSeriesN
  es_nfz = grv_bhEwaldSeriesN
  es_radius2 = grv_bhEwaldSeriesN*grv_bhEwaldSeriesN

! axis ratio of ellipsis
  cr1 = 1.0
  cr2 = 1.0
  cr3 = 1.0
  cf1 = 1.0
  cf2 = 1.0
  cf3 = 1.0

! particular values for problem with given kind of boundary conditions and orientation
! note grv_bhEwald_periodicity = 1,2,4 means 1 direction periodic and 2 isolated,
! grv_bhEwald_periodicity = 3,5,6 means 2 directions periodic and 1 isolated,
! grv_bhEwald_periodicity = 7 means boundary conditions isolated in tree directions 
  if (grv_bhEwald_periodicity == 1) then
    Linv = 1.0/grv_bhLx 
    es_nry = 0
    es_nrz = 0
    es_nfy = 0
    es_nfz = 0
  else if (grv_bhEwald_periodicity == 2) then
    Linv = 1.0/grv_bhLy 
    es_nrx = 0
    es_nrz = 0
    es_nfx = 0
    es_nfz = 0
  else if (grv_bhEwald_periodicity == 4) then
    Linv = 1.0/grv_bhLz
    es_nrx = 0
    es_nry = 0
    es_nfx = 0
    es_nfy = 0
  else if (grv_bhEwald_periodicity == 6) then
    Linv = 1.0/grv_bhLy
    ratio_p = grv_bhLz/grv_bhLy
    es_nrx = 0
    es_nrz = ceiling(grv_bhEwaldSeriesN/ratio_p)
    es_nfx = 0
    es_nfz = ceiling(grv_bhEwaldSeriesN*ratio_p)
    cr3 = ratio_p**2
    cf3 = 1.0/(ratio_p**2)
  else if (grv_bhEwald_periodicity == 5) then
    Linv = 1.0/grv_bhLz
    ratio_p = grv_bhLx/grv_bhLz
    es_nry = 0
    es_nrx = ceiling(grv_bhEwaldSeriesN/ratio_p)
    es_nfy = 0
    es_nfx = ceiling(grv_bhEwaldSeriesN*ratio_p)
    cr1 = ratio_p**2
    cf1 = 1.0/(ratio_p**2)
  else if (grv_bhEwald_periodicity == 3) then
    Linv = 1.0/grv_bhLx
    ratio_p = grv_bhLy/grv_bhLx
    es_nry = ceiling(grv_bhEwaldSeriesN/ratio_p)
    es_nrz = 0
    es_nfy = ceiling(grv_bhEwaldSeriesN*ratio_p)
    es_nfz = 0
    cr2 = ratio_p**2
    cf2 = 1.0/(ratio_p**2)
  else if (grv_bhEwald_periodicity == 7) then
! this case we haven't finished yet 
    Linv = 1.0/grv_bhLx 
    ratio_p1 = grv_bhLy/grv_bhLx
    ratio_p2 = grv_bhLz/grv_bhLx
    es_nry = ceiling(grv_bhEwaldSeriesN/ratio_p1)
    es_nrz = ceiling(grv_bhEwaldSeriesN/ratio_p2)
    es_nfy = ceiling(grv_bhEwaldSeriesN*ratio_p1)
    es_nfz = ceiling(grv_bhEwaldSeriesN*ratio_p2)
    cr2 = ratio_p1**2
    cr3 = ratio_p2**2
    cf2 = 1.0/(ratio_p1**2)
    cf3 = 1.0/(ratio_p2**2)
    ratio_pinv1 = 1.0/ratio_p1
    ratio_pinv2 = 1.0/ratio_p2
  endif

  ratio_pinv = 1.0/ratio_p
  !ewald_alpha = 2.0*Linv*ratio_pinv
  ewald_alpha = 2.0*Linv
  ewald_dzeta = pi*pi*Linv*Linv/(ewald_alpha*ewald_alpha)
 
! compute Ewald field
  ! size of a chunk of data for a given processor
  chunk = 1 + ((nz+1)*(ny+1)*(nx+1) / grv_meshNumProcs)

  ! on all CPUs set each element to zero at first
  field_EwaldRawLoc(:,:,:,:) = 0.0D0

!  grv_bhDX = ewald_xmax / nx
  grv_bhDX = 1.0D0/grv_bhDxI

  do k = 0, nz
    do j = 0, ny
      do i = 0, nx

        ! calculate the 1D index: 0..ewald_field_z*ny*ewald_field_x-1
        i1d = k*(ny+1)*(nx+1) + j*(nx+1) + i
        ! check if this point should be calculated on this processor
        if ((i1d >= grv_meshMe*chunk) .and. (i1d < (grv_meshMe+1)*chunk)) then

            ! coordinates of the point
!          x = i * ewald_xmax / nx
!          y = j * ewald_ymax / nx
!          z = k * ewald_xmax / nx

          x = i * grv_bhDX
          y = j * grv_bhDX
          z = k * grv_bhDX

          ! first term - short range interactions
          do ni = -es_nrx,es_nrx
            do nj = -es_nry,es_nry
              do nk = -es_nrz,es_nrz
              ! terms with non-negligible contributions must lie inside ellipse
                if ((cr1*ni*ni+cr2*nj*nj+cr3*nk*nk) <= es_radius2) then
                  xni = x + ni*grv_bhLx
                  yni = y + nj*grv_bhLy
                  zni = z + nk*grv_bhLz

                  field_EwaldRawLoc(:,i,j,k) = field_EwaldRawLoc(:,i,j,k) + & 
                  & grv_bhAccShort(ewald_alpha,xni,yni,zni)

                endif
              enddo
            enddo
          enddo


          ! second term - long range interactions
          do hi = -es_nfx,es_nfx
            do hj = -es_nfy,es_nfy
              do hk = -es_nfz,es_nfz
                ! terms with non-negligible contributions must lie inside ellipse 
                ! (perpendicular to the previous ellipse)
                if ((cf1*hi*hi+cf2*hj*hj+cf3*hk*hk) <= es_radius2) then

                  select case (grv_bhEwald_periodicity)
                    case (1)
                      ewald_eta = 2*pi*Linv*sqrt(y**2+z**2)
                      ewald_row = grv_bhAccLong1P2I(Linv, ewald_dzeta, ewald_eta, hi, x, y, z)

                    case (2)
                      ewald_eta = 2*pi*Linv*sqrt(x**2+z**2)
!                      ewald_rowPrep = grv_bhAccLong1P2I(Linv, ewald_dzeta, ewald_eta, hj, y, -x, z)
                      ewald_rowPrep = grv_bhAccLong1P2I(Linv, ewald_dzeta, ewald_eta, hj, y, x, z)

                      ewald_row(0) = ewald_rowPrep(0)
                      ewald_row(1) = ewald_rowPrep(2)
                      ewald_row(2) = ewald_rowPrep(1)
                      ewald_row(3) = ewald_rowPrep(3)
                      ewald_row(4) = ewald_rowPrep(8)
                      ewald_row(5) = ewald_rowPrep(5)
                      ewald_row(6) = ewald_rowPrep(9)
                      ewald_row(8) = ewald_rowPrep(4)
                      ewald_row(9) = ewald_rowPrep(6)
                      ewald_row(12) = ewald_rowPrep(12)

                    case (4)
                      ewald_eta = 2*pi*Linv*sqrt(x**2+y**2)
!                      ewald_rowPrep = grv_bhAccLong1P2I(Linv, ewald_dzeta, ewald_eta, hk, z, -x, -y)
                      ewald_rowPrep = grv_bhAccLong1P2I(Linv, ewald_dzeta, ewald_eta, hk, z, x, y)

                      ewald_row(0) = ewald_rowPrep(0)
                      ewald_row(1) = ewald_rowPrep(2)
                      ewald_row(2) = ewald_rowPrep(3)
                      ewald_row(3) = ewald_rowPrep(1)
                      ewald_row(4) = ewald_rowPrep(8)
                      ewald_row(5) = ewald_rowPrep(9)
                      ewald_row(6) = ewald_rowPrep(5)
                      ewald_row(8) = ewald_rowPrep(12)
                      ewald_row(9) = ewald_rowPrep(6)
                      ewald_row(12) = ewald_rowPrep(4)

                    case(3)
                      ewald_row = grv_bhAccLong2P1I(Linv, ewald_dzeta, ratio_pinv, hi, hj, x, y, z)
                    case(5)
                      ewald_rowPrep = grv_bhAccLong2P1I(Linv, ewald_dzeta, ratio_pinv, hi, hk, z, x, y)

                      ewald_row(0) = ewald_rowPrep(0)
                      ewald_row(1) = ewald_rowPrep(2)
                      ewald_row(2) = ewald_rowPrep(3)
                      ewald_row(3) = ewald_rowPrep(1)
                      ewald_row(4) = ewald_rowPrep(8)
                      ewald_row(5) = ewald_rowPrep(9)
                      ewald_row(6) = ewald_rowPrep(5)
                      ewald_row(8) = ewald_rowPrep(12)
                      ewald_row(9) = ewald_rowPrep(6)
                      ewald_row(12) = ewald_rowPrep(4)

                    case(6)
                      ewald_rowPrep = grv_bhAccLong2P1I(Linv, ewald_dzeta, ratio_pinv, hk, hj, y, z, x)

                      ewald_row(0) = ewald_rowPrep(0)
                      ewald_row(1) = ewald_rowPrep(3)
                      ewald_row(2) = ewald_rowPrep(1)
                      ewald_row(3) = ewald_rowPrep(2)
                      ewald_row(4) = ewald_rowPrep(12)
                      ewald_row(5) = ewald_rowPrep(6)
                      ewald_row(6) = ewald_rowPrep(9)
                      ewald_row(8) = ewald_rowPrep(4)
                      ewald_row(9) = ewald_rowPrep(5)
                      ewald_row(12) = ewald_rowPrep(8)

                    case (7)
                      ewald_row = grv_bhAccLong3P(Linv, ewald_dzeta, ratio_pinv1, ratio_pinv2, hi, hj, hk, x, y, z)

                  end select

                  field_EwaldRawLoc(:,i,j,k) = field_EwaldRawLoc(:,i,j,k) + ewald_row

                endif
              enddo
            enddo
          enddo

! subtract the 1/r or similar terms

           rni2 = x**2+y**2+z**2
           rni = sqrt(rni2)

           field_EwaldRawLoc(0,i,j,k) = field_EwaldRawLoc(0,i,j,k) - 1.0D0/rni

           field_EwaldRawLoc(1,i,j,k) = field_EwaldRawLoc(1,i,j,k) - x/(rni**3)
           field_EwaldRawLoc(2,i,j,k) = field_EwaldRawLoc(2,i,j,k) - y/(rni**3)
           field_EwaldRawLoc(3,i,j,k) = field_EwaldRawLoc(3,i,j,k) - z/(rni**3)

           field_EwaldRawLoc(4,i,j,k) = field_EwaldRawLoc(4,i,j,k) - (rni2-3*x*x)/(rni**5)
           field_EwaldRawLoc(5,i,j,k) = field_EwaldRawLoc(5,i,j,k) + 3*x*y/(rni**5)
           field_EwaldRawLoc(6,i,j,k) = field_EwaldRawLoc(6,i,j,k) + 3*x*z/(rni**5)
           field_EwaldRawLoc(8,i,j,k) = field_EwaldRawLoc(8,i,j,k) - (rni2-3*y*y)/(rni**5)
           field_EwaldRawLoc(9,i,j,k) = field_EwaldRawLoc(9,i,j,k) + 3*y*z/(rni**5)
           field_EwaldRawLoc(12,i,j,k) = field_EwaldRawLoc(12,i,j,k) - (rni2-3*z*z)/(rni**5)

! use symmetry (since second derivatives do not depend on the order of derivation)
          field_EwaldRawLoc(7,i,j,k) = field_EwaldRawLoc(5,i,j,k)
          field_EwaldRawLoc(10,i,j,k) = field_EwaldRawLoc(6,i,j,k)
          field_EwaldRawLoc(11,i,j,k) = field_EwaldRawLoc(9,i,j,k)

        endif ! end my chunk

      enddo
    enddo
  enddo

  call MPI_AllReduce(field_EwaldRawLoc,field_EwaldRaw, & 
  &    13*(nx+1)*(ny+1)*(nz+1),FLASH_REAL,FLASH_SUM,grv_meshComm,ierr)  

! set appropriate field at the origin (0,0,0) to zero
  field_EwaldRaw(:,0,0,0) = 0.0D0

! determine difference in gravitational potential between numerical and analytical
! value (it is needed for interpolating)
  select case (grv_bhEwald_periodicity)
    case(1)
      y = grv_bhDX*ny
      z = grv_bhDX*nz
      ew_rcmax = sqrt(y**2+z**2)
      grv_bhPotConst = 2.0D0*Linv*log(ew_rcmax) + &
      &                field_EwaldRaw(0,0,ny,nz) + 1.0D0/ew_rcmax

    case(2)
      x = grv_bhDX*nx
      z = grv_bhDX*nz
      ew_rcmax = sqrt(x**2+z**2)
      grv_bhPotConst = 2.0D0*Linv*log(ew_rcmax) + &
      &                field_EwaldRaw(0,nx,0,nz) + 1.0D0/ew_rcmax

    case(4)
      x = grv_bhDX*nx
      y = grv_bhDX*ny
      ew_rcmax = sqrt(x**2 + y**2)
      grv_bhPotConst = 2.0D0*Linv*log(ew_rcmax) + &
      &                field_EwaldRaw(0,nx,ny,0) + 1.0D0/ew_rcmax

    case(3)
      z = grv_bhDX*nz

      grv_bhPotConst = 2.0D0*pi*z/(grv_bhLx*grv_bhLy) + &
      &                field_EwaldRaw(0,0,0,nz) + 1.0D0/z

    case(5)
      y = grv_bhDX*ny
      grv_bhPotConst = 2.0D0*pi*y/(grv_bhLx*grv_bhLz) + &
      &                field_EwaldRaw(0,0,ny,0) + 1.0D0/y

    case(6)
      x = grv_bhDX*nx
      grv_bhPotConst = 2.0D0*pi*x/(grv_bhLy*grv_bhLz) + &
      &                field_EwaldRaw(0,nx,0,0) + 1.0D0/x

  end select


! reorganise data to get some speedup during interpolation

  do k = 0, nz
    do j = 0, ny
      do i = 0, nx

        x = i * grv_bhDX
        y = j * grv_bhDX
        z = k * grv_bhDX

! potential plus shift (this is compensation for offset, because we use x
! instead of dx in routines Gravity_bhEwaldAcc.F90 and
! Gravity_bhEwaldPot.F90)
        field_Ewald(0,i,j,k) = field_EwaldRaw(0,i,j,k) + x*field_EwaldRaw(1,i,j,k) + &
        &    y*field_EwaldRaw(2,i,j,k) + z*field_EwaldRaw(3,i,j,k)
! partial derivatives of potential
        field_Ewald(1,i,j,k) = -field_EwaldRaw(1,i,j,k)
        field_Ewald(2,i,j,k) = -field_EwaldRaw(2,i,j,k)
        field_Ewald(3,i,j,k) = -field_EwaldRaw(3,i,j,k)
! acceleration plus shift
! note that now the values of field for acceleration are very different from the 
! values of gradient of the potential
        field_Ewald(4,i,j,k) = field_EwaldRaw(1,i,j,k) - x*field_EwaldRaw(4,i,j,k) - &
        & y*field_EwaldRaw(5,i,j,k) - z*field_EwaldRaw(6,i,j,k)
        field_Ewald(5,i,j,k) = field_EwaldRaw(2,i,j,k) - x*field_EwaldRaw(5,i,j,k) - &
        & y*field_EwaldRaw(8,i,j,k) - z*field_EwaldRaw(9,i,j,k)
        field_Ewald(6,i,j,k) = field_EwaldRaw(3,i,j,k) - x*field_EwaldRaw(6,i,j,k) - &
        & y*field_EwaldRaw(9,i,j,k) - z*field_EwaldRaw(12,i,j,k)
! partial derivatives of potential - only 6 terms instead of 9 due to symmetry 
! of second partial derivatives
        field_Ewald(7,i,j,k) = field_EwaldRaw(4,i,j,k)
        field_Ewald(8,i,j,k) = field_EwaldRaw(5,i,j,k)
        field_Ewald(9,i,j,k) = field_EwaldRaw(6,i,j,k)
        field_Ewald(10,i,j,k) = field_EwaldRaw(8,i,j,k)
        field_Ewald(11,i,j,k) = field_EwaldRaw(9,i,j,k)
        field_Ewald(12,i,j,k) = field_EwaldRaw(12,i,j,k)

      enddo
    enddo
  enddo
  

  return
end subroutine grv_bhGenEwaldFieldLevel
