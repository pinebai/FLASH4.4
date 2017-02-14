!!****if* source/Grid/GridSolvers/Multigrid/gr_hgPoissonSolve3D
!!
!! NAME
!!  gr_hgPoissonSolve3D
!!
!! SYNOPSIS
!!
!!  call gr_hgPoissonSolve3D(real(INOUT) :: soln,
!!                   integer(IN)    :: nx,
!!                   integer(IN)    :: ny,
!!                   integer(IN)    :: nz,
!!                      real(IN)    :: dx,
!!                      real(IN)    :: dy,
!!                      real(IN)    :: dz,
!!                   integer(IN)    :: ibnd,
!!                      real(INOUT) :: Gx(0:nx),
!!                      real(INOUT) :: Gy(0:ny),
!!                      real(INOUT) :: Gz(0:nz),
!!                      real(INOUT) :: wk1(0:nwk1),
!!                      real(INOUT) :: wk2(0:nwk2),
!!                      real(INOUT) :: wk3(0:nwk3),
!!                   integer(INOUT) :: ip(0:nip)
!!                      real(IN)    :: norm,
!!                   integer(IN)    :: nwk1,
!!                   integer(IN)    :: nwk2,
!!                   integer(IN)    :: nwk3,
!!                   integer(IN)    :: nip,
!!                                          )
!!
!!
!! DESCRIPTION
!!
!! Solve Poisson's equation on a 3D rectangular domain using
!! Fourier, sine, or cosine transforms.  The type of transform
!! employed depends upon the requested boundary conditions.
!! Transforms are performed using Takuya Ooura's FFT package
!! (http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html).
!!
!! ARGUMENTS
!!
!!   soln(:)            -     On input, source array, sampled at zone centers;
!!                  on output, array to receive solution, sampled at zone centers
!!   nx,ny,nz           -     Dimensions of source and solution arrays
!!   dx,dy,dz           -     Zone widths in x, y, z
!!   ibnd               -   Type of boundary conditions to assume:
!!                           MG_BND_PERIODIC, MG_BND_DIRICHLET, MG_BND_NEUMANN
!!   Gx, Gy, Gz         -     The Green's functions in x, y, and z
!!   wk1,wk2,wk3        -     Probably space for a work array
!!   norm               -     A scaling factor for the solution
!!   nwk1,nwk2,nwk3     -     Size of work array
!!   ip                 -     work area for bit reversal.  
!!                            ip(0), ip(1) are pointers of the cos/sin table
!!   nip                -     Size of ip array
!!
!!
!!***


subroutine gr_hgPoissonSolve3D (soln, nx, ny, nz, dx, dy, dz, ibnd, &
     Gx, Gy, Gz, wk1, wk2, ip, norm, nwk1, nwk2, nip, G3D)

#include "Multigrid.h"

  !===============================================================================

  implicit none

  real, intent(in)       :: dx, dy, dz, norm
  integer, intent(in)    :: nx, ny, nz, ibnd, nwk1, nwk2, nip
  real, intent(inout)    :: soln(0:nx-1,0:ny-1,0:nz-1)
  real, intent(inout)    :: wk1(0:nwk1), wk2(0:nwk2), Gx(0:nx), Gy(0:ny), Gz(0:nz)
  integer, intent(inout) :: ip(0:nip)
  real, intent(inout)    :: G3D(0:nx,0:ny,0:nz)

  integer         :: i, j, k
  real            :: G

  !===============================================================================

  ! Solve the Poisson equation by transforming the source, applying the Green's
  ! function, and performing the inverse transform.

  select case (ibnd)

     !===============================================================================

     ! Periodic boundary conditions:  use real discrete Fourier transform

  case (MG_BND_PERIODIC)

     call rdft3d (nx, ny, nx, ny, nz, 1, soln, wk1, ip, wk2)

     do k = 0, nz-1
        do j = 0, ny-1
           do i = 1, nx/2-1
              !          G = GreenFctn3D(Gx(i), Gy(j), Gz(k), ibnd)
              G = G3D(i,j,k)
              soln(2*i,j,k)   = soln(2*i,j,k) * G
              soln(2*i+1,j,k) = soln(2*i+1,j,k) * G
           enddo
        enddo
     enddo

     do k = 0, nz-1
        do j = 1, ny/2-1
           i = 0
           !        G = GreenFctn3D(Gx(i), Gy(j), Gz(k), ibnd)
           G = G3D(i,j,k)
           soln(0,j,k) = soln(0,j,k) * G
           soln(1,j,k) = soln(1,j,k) * G
           i = nx/2
           !        G = GreenFctn3D(Gx(i), Gy(j), Gz(k), ibnd)
           G = G3D(i,j,k)
           soln(1,ny-j,k) = soln(1,ny-j,k) * G
           soln(0,ny-j,k) = soln(0,ny-j,k) * G
        enddo
     enddo

     do k = 1, nz/2-1
        i = 0
        j = 0
        !      G = GreenFctn3D(Gx(i), Gy(j), Gz(k), ibnd)
        G = G3D(i,j,k)
        soln(0,0,k) = soln(0,0,k) * G
        soln(1,0,k) = soln(1,0,k) * G
        j = ny/2
        !      G = GreenFctn3D(Gx(i), Gy(j), Gz(k), ibnd)
        G = G3D(i,j,k)
        soln(0,ny/2,k) = soln(0,ny/2,k) * G
        soln(1,ny/2,k) = soln(1,ny/2,k) * G
        i = nx/2
        j = 0
        !      G = GreenFctn3D(Gx(i), Gy(j), Gz(k), ibnd)
        G = G3D(i,j,k)
        soln(1,0,nz-k) = soln(1,0,nz-k) * G
        soln(0,0,nz-k) = soln(0,0,nz-k) * G
        j = ny/2
        !      G = GreenFctn3D(Gx(i), Gy(j), Gz(k), ibnd)
        G = G3D(i,j,k)
        soln(1,ny/2,nz-k) = soln(1,ny/2,nz-k) * G
        soln(0,ny/2,nz-k) = soln(0,ny/2,nz-k) * G
     enddo

     !    G = GreenFctn3D(Gx(0), Gy(0), Gz(0), ibnd)
     G = G3D(0,0,0)
     soln(0,0,0) = soln(0,0,0) * G
     !    G = GreenFctn3D(Gx(nx/2), Gy(0), Gz(0), ibnd)
     G = G3D(nx/2,0,0)
     soln(1,0,0) = soln(1,0,0) * G
     !    G = GreenFctn3D(Gx(0), Gy(0), Gz(nz/2), ibnd)
     G = G3D(0,0,nz/2)
     soln(0,0,nz/2) = soln(0,0,nz/2) * G
     !    G = GreenFctn3D(Gx(nx/2), Gy(0), Gz(nz/2), ibnd)
     G = G3D(nx/2,0,nz/2)
     soln(1,0,nz/2) = soln(1,0,nz/2) * G
     !    G = GreenFctn3D(Gx(0), Gy(ny/2), Gz(0), ibnd)
     G = G3D(0,ny/2,0)
     soln(0,ny/2,0) = soln(0,ny/2,0) * G
     !    G = GreenFctn3D(Gx(nx/2), Gy(ny/2), Gz(0), ibnd)
     G = G3D(nx/2,ny/2,0)
     soln(1,ny/2,0) = soln(1,ny/2,0) * G
     !    G = GreenFctn3D(Gx(0), Gy(ny/2), Gz(nz/2), ibnd)
     G = G3D(0,ny/2,nz/2)
     soln(0,ny/2,nz/2) = soln(0,ny/2,nz/2) * G
     !    G = GreenFctn3D(Gx(nx/2), Gy(ny/2), Gz(nz/2), ibnd)
     G = G3D(nx/2,ny/2,nz/2)
     soln(1,ny/2,nz/2) = soln(1,ny/2,nz/2) * G

     call rdft3d (nx, ny, nx, ny, nz, -1, soln, wk1, ip, wk2)

     do k = 0, nz-1
        do j = 0, ny-1
           do i = 0, nx-1
              soln(i,j,k) = soln(i,j,k) * norm
           enddo
        enddo
     enddo

     !===============================================================================

     ! Dirichlet boundary conditions:  use discrete sine transform

  case (MG_BND_DIRICHLET)

     call ddst3d (nx, ny, nx, ny, nz, -1, soln, wk1, ip, wk2)
     do k = 1, nz
        do j = 1, ny
           do i = 1, nx
              !          G = GreenFctn3D(Gx(i), Gy(j), Gz(k), ibnd)
              G = G3D(i,j,k)
              soln(mod(i,nx),mod(j,ny),mod(k,nz)) = &
                   soln(mod(i,nx),mod(j,ny),mod(k,nz)) * G
           enddo
        enddo
     enddo
     do j = 0, ny-1
        do i = 0, nx-1
           soln(i,j,0) = soln(i,j,0) * 0.5
        enddo
     enddo

     do k = 0, nz-1
        do i = 0, nx-1
           soln(i,0,k) = soln(i,0,k) * 0.5
        enddo
        do j = 0, ny-1
           soln(0,j,k) = soln(0,j,k) * 0.5
        enddo
     enddo

     call ddst3d (nx, ny, nx, ny, nz, 1, soln, wk1, ip, wk2)

     do k = 0, nz-1
        do j = 0, ny-1
           do i = 0, nx-1
              soln(i,j,k) = soln(i,j,k) * norm
           enddo
        enddo
     enddo

     !===============================================================================

     ! Neumann boundary conditions:  use discrete cosine transform

  case (MG_BND_NEUMANN)

     call ddct3d (nx, ny, nx, ny, nz, -1, soln, wk1, ip, wk2)

     do k = 0, nz-1
        do j = 0, ny-1
           do i = 0, nx-1
              !          G = GreenFctn3D(Gx(i), Gy(j), Gz(k), ibnd)
              G = G3D(i,j,k)
              soln(i,j,k) = soln(i,j,k) * G
           enddo
        enddo
     enddo

     do j = 0, ny-1
        do i = 0, nx-1
           soln(i,j,0) = soln(i,j,0) * 0.5
        enddo
     enddo

     do k = 0, nz-1
        do i = 0, nx-1
           soln(i,0,k) = soln(i,0,k) * 0.5
        enddo
        do j = 0, ny-1
           soln(0,j,k) = soln(0,j,k) * 0.5
        enddo
     enddo

     call ddct3d (nx, ny, nx, ny, nz, 1, soln, wk1, ip, wk2)

     do k = 0, nz-1
        do j = 0, ny-1
           do i = 0, nx-1
              soln(i,j,k) = soln(i,j,k) * norm
           enddo
        enddo
     enddo

     !===============================================================================

  end select

  !===============================================================================

  return
end subroutine gr_hgPoissonSolve3D
