!!****if* source/Grid/GridSolvers/Multigrid/gr_hgPoissonSolve1D
!!
!! NAME
!!  gr_hgPoissonSolve1D
!!
!! SYNOPSIS
!!
!!  gr_hgPoissonSolve1D(real(INOUT) :: soln,
!!                   integer(IN)    :: nx,
!!                      real(IN)    :: dx,
!!                   integer(IN)    :: ibnd,
!!                      real(INOUT) :: Gx(0:nx),
!!                      real(INOUT) :: wk1(0:nwk1),
!!                   integer(INOUT) :: ip(0:nip)
!!                      real(IN)    :: norm,
!!                   integer(IN)    :: nwk1,
!!                   integer(IN)    :: nip)
!!
!!
!! DESCRIPTION
!!
!!  Description: Solve Poisson's equation on a 1D rectangular domain using
!!               Fourier, sine, or cosine transforms.  The type of transform
!!               employed depends upon the requested boundary conditions.
!!               Transforms are performed using Takuya Ooura's FFT package
!!               (http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html).
!!
!! ARGUMENTS
!!
!!   soln(:)  -     On input, source array, sampled at zone centers;
!!                  on output, array to receive solution, sampled at zone centers
!!   nx       -     Dimensions of source and solution arrays
!!   dx       -     Zone widths in x
!!   ibnd     -     Type of boundary conditions to assume:
!!                    MG_BND_PERIODIC, MG_BND_DIRICHLET, MG_BND_NEUMANN
!!   Gx       -     The Green's function for the x direction
!!   wk1      -     Probably space for a work array
!!   ip       -     Additional space used internally for the FFT solve
!!   norm     -     A scaling factor for the solution
!!   nwk1     -     Size of work array
!!   nip      -     Size of ip array
!!
!!
!!***


subroutine gr_hgPoissonSolve1D (soln, nx, dx, ibnd, &
     Gx, wk1, ip, norm, nwk1, nip)

  !===============================================================================

#include "Multigrid.h"

  implicit none

  real, intent(in)       :: dx, norm
  integer, intent(in)    :: nx, ibnd, nwk1, nip
  real, intent(inout)    :: soln(0:nx-1)
  real, intent(inout)    :: wk1(0:nwk1), Gx(0:nx)
  integer, intent(inout) :: ip(0:nip)

  !internal functions
  real gr_hgPoissonGreenFctn1D

  integer         :: i
  real            :: G

  !===============================================================================

  ! Solve the Poisson equation by transforming the source, applying the Green's
  ! function, and performing the inverse transform.

  select case (ibnd)

     !===============================================================================

     ! Periodic boundary conditions:  use real discrete Fourier transform

  case (MG_BND_PERIODIC)

     call rdft (nx, 1, soln, ip, wk1)

     do i = 0, nx/2-1
        G = gr_hgPoissonGreenFctn1D(Gx(i), ibnd)
        soln(2*i) = soln(2*i) * G
     enddo

     do i = 1, nx/2-1
        G = gr_hgPoissonGreenFctn1D(Gx(i), ibnd)
        soln(2*i+1) = soln(2*i+1) * G
     enddo

     G = gr_hgPoissonGreenFctn1D(Gx(nx/2), ibnd)
     soln(1) = soln(1) * G


     call rdft (nx, -1, soln, ip, wk1)

     do i = 0, nx-1
        soln(i) = soln(i) * norm
     enddo

     !===============================================================================

     ! Dirichlet boundary conditions:  use discrete sine transform

  case (MG_BND_DIRICHLET)
     call ddst (nx, -1, soln, ip, wk1)

     do i = 1, nx
        G = gr_hgPoissonGreenFctn1D(Gx(i), ibnd)
        soln(mod(i,nx)) = soln(mod(i,nx)) * G
     enddo

     soln(0) = soln(0) * 0.5

     call ddst (nx, 1, soln, ip, wk1)

     do i = 0, nx-1
        soln(i) = soln(i) * norm
     enddo

     !===============================================================================

     ! Neumann boundary conditions:  use discrete cosine transform

  case (MG_BND_NEUMANN)
     call ddct (nx, -1, soln, ip, wk1)

     do i = 0, nx-1
        G = gr_hgPoissonGreenFctn1D(Gx(i), ibnd)
        soln(i) = soln(i) * G
     enddo

     soln(0) = soln(0) * 0.5

     call ddct (nx, 1, soln, ip, wk1)

     do i = 0, nx-1
        soln(i) = soln(i) * norm
     enddo

     !===============================================================================

  end select

  !===============================================================================

  return
end subroutine gr_hgPoissonSolve1D

!*******************************************************************************

!  Routine:     GreenFctn1D

!  Description: Evaluate the Green's function for the Poisson equation at a
!               specified k-value.  The denominator is accepted as an argument.


function gr_hgPoissonGreenFctn1D (Gx, ibnd)

  !===============================================================================

  implicit none

  real, intent(in)    :: Gx
  integer, intent(in) :: ibnd
  real                :: gr_hgPoissonGreenFctn1D

  real                :: Ginv

  !===============================================================================

  Ginv = Gx

  if (abs(Ginv) > 1.D-99) then
     gr_hgPoissonGreenFctn1D = 1./(Ginv)
  else
     gr_hgPoissonGreenFctn1D = 0.
  endif

  !===============================================================================

  return
end function gr_hgPoissonGreenFctn1D


