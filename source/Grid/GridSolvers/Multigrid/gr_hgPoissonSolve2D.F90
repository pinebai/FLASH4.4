!!****if* source/Grid/GridSolvers/Multigrid/gr_hgPoissonSolve2D
!!
!! NAME
!!  gr_hgPoissonSolve2D
!!
!! SYNOPSIS
!!
!!  gr_hgPoissonSolve2D(real(INOUT) :: soln,
!!                   integer(IN)    :: nx,
!!                   integer(IN)    :: ny,
!!                      real(IN)    :: dx,
!!                      real(IN)    :: dy,
!!                   integer(IN)    :: ibnd,
!!                      real(INOUT) :: Gx(0:nx),
!!                      real(INOUT) :: Gy(0:ny),
!!                      real(INOUT) :: wk1(0:nwk1),
!!                      real(INOUT) :: wk2(0:nwk2),
!!                   integer(INOUT) :: ip(0:nip)
!!                      real(IN)    :: norm,
!!                   integer(IN)    :: nwk1,
!!                   integer(IN)    :: nwk2,
!!                   integer(IN)    :: nip,
!!                                          )
!!
!!
!! DESCRIPTION
!!
!! Solve Poisson's equation on a 2D rectangular domain using
!! Fourier, sine, or cosine transforms.  The type of transform
!! employed depends upon the requested boundary conditions.
!! Transforms are performed using Takuya Ooura's FFT package
!! (http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html).
!!
!! ARGUMENTS
!!
!!   soln(:)  -     On input, source array, sampled at zone centers;
!!                  on output, array to receive solution, sampled at zone centers
!!   nx,ny    -     Dimensions of source and solution arrays
!!   dx,dy    -     Zone widths in x, y
!!   ibnd     -     Type of boundary conditions to assume:
!!                    MG_BND_PERIODIC, MG_BND_DIRICHLET, MG_BND_NEUMANN
!!   Gx,Gy    -     The Green's functions in the x and y directions
!!   wk1,wk2  -     Probably space for a work array
!!   ip       -     Additional space used internally for the FFT solve
!!   norm     -     A scaling factor for the solution
!!   nwk1, nwk2     -     Size of work array
!!   nip      -     Size of ip array
!!
!!  NOTES
!!    This routine is called from gr_hgPoissonSolveBlock
!!
!!
!!***


subroutine gr_hgPoissonSolve2D (soln, nx, ny, dx, dy, ibnd, &
     Gx, Gy, wk1, wk2, ip, norm, nwk1, nwk2, nip, G2D)

  !===============================================================================
#include "Multigrid.h"

  implicit none

  real, intent(in)       :: dx, dy, norm
  integer, intent(in)    :: nx, ny, ibnd, nwk1, nwk2, nip
  real, intent(inout)    :: soln(0:nx-1,0:ny-1)
  real, intent(inout)    :: wk1(0:nwk1), wk2(0:nwk2), Gx(0:nx), Gy(0:ny)
  integer, intent(inout) :: ip(0:nip)
  real, intent(inout)    :: G2D(0:nx,0:ny)

  integer         :: i, j
  real            :: G

  !===============================================================================

  ! Solve the Poisson equation by transforming the source, applying the Green's
  ! function, and performing the inverse transform.

  select case (ibnd)

     !===============================================================================

     ! Periodic boundary conditions:  use real discrete Fourier transform

  case (MG_BND_PERIODIC)

     ! Forward transform
     call rdft2d (nx, nx, ny, 1, soln, wk1, ip, wk2)

     ! Apply Green's function
     do j = 1, ny-1
        do i = 1, nx/2-1
           !        G = GreenFctn2D(Gx(i), Gy(j), ibnd)
           G = G2D(i,j)
           soln(2*i,j)   = soln(2*i,j) * G
           soln(2*i+1,j) = soln(2*i+1,j) * G
        enddo
     enddo

     do i = 1, nx/2-1
        !      G = GreenFctn2D(Gx(i), Gy(0), ibnd)
        G = G2D(i,0)
        soln(2*i,0)   = soln(2*i,0) * G
        soln(2*i+1,0) = soln(2*i+1,0) * G
     enddo

     do j = 1, ny/2-1
        !      G = GreenFctn2D(Gx(0), Gy(j), ibnd)
        G = G2D(0,j)
        soln(0,j) = soln(0,j) * G
        soln(1,j) = soln(1,j) * G
        !      G = GreenFctn2D(Gx(nx/2), Gy(j), ibnd)
        G = G2D(nx/2,j)
        soln(1,ny-j) = soln(1,ny-j) * G
        soln(0,ny-j) = soln(0,ny-j) * G
     enddo

     !    G = GreenFctn2D(Gx(0), Gy(0), ibnd)
     G = G2D(0,0)
     soln(0,0) = soln(0,0) * G
     !    G = GreenFctn2D(Gx(nx/2), Gy(0), ibnd)
     G = G2D(nx/2,0)
     soln(1,0) = soln(1,0) * G
     !    G = GreenFctn2D(Gx(0), Gy(ny/2), ibnd)
     G = G2D(0,ny/2)
     soln(0,ny/2) = soln(0,ny/2) * G
     !    G = GreenFctn2D(Gx(nx/2), Gy(ny/2), ibnd)
     G = G2D(nx/2,ny/2)
     soln(1,ny/2) = soln(1,ny/2) * G

     ! Inverse transform
     call rdft2d (nx, nx, ny, -1, soln, wk1, ip, wk2)

     ! NOTE -- norm is calculated in gr_hgPoissonSolveBlock as
     !        norm(level) = 2./(float(nx)*float(ny)*float(nz))
     do j = 0, ny-1
        do i = 0, nx-1
           soln(i,j) = soln(i,j) * norm
        enddo
     enddo

     !===============================================================================

     ! Dirichlet boundary conditions:  use discrete sine transform

  case (MG_BND_DIRICHLET)

     call ddst2d (nx, nx, ny, -1, soln, wk1, ip, wk2)

     do j = 1, ny
        do i = 1, nx
           !        G = GreenFctn2D(Gx(i), Gy(j), ibnd)
           G = G2D(i,j)
           soln(mod(i,nx),mod(j,ny)) = soln(mod(i,nx),mod(j,ny)) * G
        enddo
     enddo

     do i = 0, nx-1
        soln(i,0) = soln(i,0) * 0.5
     enddo

     do j = 0, ny-1
        soln(0,j) = soln(0,j) * 0.5
     enddo

     call ddst2d (nx, nx, ny, 1, soln, wk1, ip, wk2)

     do j = 0, ny-1
        do i = 0, nx-1
           soln(i,j) = soln(i,j) * norm
        enddo
     enddo

     !===============================================================================

     ! Neumann boundary conditions:  use discrete cosine transform

  case (MG_BND_NEUMANN)

     call ddct2d (nx, nx, ny, -1, soln, wk1, ip, wk2)

     do j = 0, ny-1
        do i = 0, nx-1
           !        G = GreenFctn2D(Gx(i), Gy(j), ibnd)
           G = G2D(i,j)
           soln(i,j) = soln(i,j) * G
        enddo
     enddo

     do i = 0, nx-1
        soln(i,0) = soln(i,0) * 0.5
     enddo

     do j = 0, ny-1
        soln(0,j) = soln(0,j) * 0.5
     enddo

     call ddct2d (nx, nx, ny, 1, soln, wk1, ip, wk2)

     do j = 0, ny-1
        do i = 0, nx-1
           soln(i,j) = soln(i,j) * norm
        enddo
     enddo

     !===============================================================================

  end select

  !===============================================================================

  return
end subroutine gr_hgPoissonSolve2D
