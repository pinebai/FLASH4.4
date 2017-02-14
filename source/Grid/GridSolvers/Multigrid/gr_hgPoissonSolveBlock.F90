!!****if* source/Grid/GridSolvers/Multigrid/gr_hgPoissonSolveBlock
!!
!! NAME
!!  gr_hgPoissonSolveBlock
!!
!! SYNOPSIS
!!  gr_hgPoissonSolveBlock(real, intent(inout) :: soln(:,:,:),
!!                         integer, intent(in) :: nx,
!!                         integer, intent(in) :: ny,
!!                         integer, intent(in) :: nz,
!!                         real, intent(in)    :: dx,
!!                         real, intent(in)    :: dy,
!!                         real, intent(in)    :: dz,
!!                         integer, intent(in) :: bnd_type,
!!                         integer, intent(in) :: level)
!!
!! DESCRIPTION
!!  
!!  Single block Poisson Solver using an FFT library.
!!
!! ARGUMENTS
!!
!!  soln       - the input and output arrays, sampled at zone centers
!!  nx, ny, nz - the dimensions of the array
!!  dx, dy, dz - the zone widths
!!  bnd_type   - MG_BND_PERIODIC, MG_BND_DIRICHLET, or MG_BND_NEUMANN
!!  level      - the current level to solve at (precalculation of dimensions)
!!
!! NOTES
!!
!! This routine is called as SolveBlock from gr_hgSolveLevel with arguments
!!    (soln,NXB,NYB,NZB,dx,dy,dz,bnd_type,level)
!!
!!***

subroutine gr_hgPoissonSolveBlock (soln, nx, ny, nz, dx, dy, dz, bnd_type,&
                                   level)

!=========================================================

!use runtime_parameters

#include "Multigrid.h"
#include "Flash.h"

#ifdef FLASH_GRID_PARAMESH2
use tree, ONLY: lrefine_max,ndim
#else
use paramesh_dimensions, ONLY: ndim
use tree, ONLY: lrefine_max
#endif

implicit none

integer, intent(in) :: nx, ny, nz, bnd_type, level
real, intent(in)    :: dx, dy, dz
real, intent(inout) :: soln(nx,ny,nz)

logical, save              :: first_call = .true.
integer, save              :: nwk1, nwk2, nip
integer, allocatable, save :: iinit(:), lastBndType(:), ip(:,:)
real, allocatable, save    :: wk1(:,:), wk2(:,:), Gx(:,:), Gy(:,:), Gz(:,:)
real, allocatable, save    :: norm(:), G2D(:,:,:), G3D(:,:,:,:)
real, parameter            :: pi = 3.1415926535898

real                       :: nxinv, nyinv, nzinv, dxinv, dyinv, dzinv
integer                    :: m, i, j, k

!===============================================================================

! No good way to remove this first call 
! DEV more importantly, how do you deallocate?
if (first_call) then
  allocate(iinit(lrefine_max))
  allocate(lastBndType(lrefine_max))
  iinit(:) = 1
  allocate(norm(lrefine_max))
  allocate(Gx(0:nx,lrefine_max))
  allocate(Gy(0:ny,lrefine_max))
  allocate(Gz(0:nz,lrefine_max))
  if (ndim == 2) then
    allocate(G2D(0:nx,0:ny,lrefine_max))
  elseif (ndim == 3) then
    allocate(G3D(0:nx,0:ny,0:nz,lrefine_max))
  endif

  if (ndim == 1) then
    nwk1 = nx*5/4-1                     ! nwk1 and nip are the problems with moving this to 
    nwk2 = 0                            ! initializations -- they are allocated below
    nip  = 1+int(sqrt(float(nx/2)))     ! but are routine arguments.
  else
!    select case (bnd_type)
!      case (0) ! periodic
        nwk1 = 8*max(ny, nz)-1
!        nwk2 = max(nx/4, ny/2, nz/2)+nx/4-1
        nip  = 1+int(sqrt(float(max(nx/2, ny, nz))))
!      case (1) ! Dirichlet
!        nwk1 = 4*max(ny, nz)-1
        nwk2 = max(nx*3/2, ny*3/2, nz*3/2)-1
!        nip  = 1+int(sqrt(float(max(nx/2, ny/2, nz/2))))
!      case (2) ! Neumann
!        nwk1 = 4*max(ny, nz)-1
!        nwk2 = max(nx*3/2, ny*3/2, nz*3/2)-1
!        nip  = 1+int(sqrt(float(max(nx/2, ny/2, nz/2))))
!    end select
  endif
  allocate (wk1(0:nwk1,lrefine_max))
  allocate (wk2(0:nwk2,lrefine_max))
  allocate (ip(0:nip,lrefine_max))
  first_call = .false.
endif

!  If level already initialized but for a different bc, force reinitialization here - KW
if (iinit(level) == 0) then
   if (bnd_type .NE. lastBndType(level)) iinit(level) = 1
end if

if (iinit(level) == 1) then
  dxinv = 2./dx**2
  dyinv = 0.0
  dzinv = 0.0
  if (NDIM >=2) dyinv = 2./dy**2
  if (NDIM ==3) dzinv = 2./dz**2
  nxinv = pi/nx
  nyinv = pi/ny
  nzinv = pi/nz
  Gx(:,level) = 0.
  Gy(:,level) = 0.
  Gz(:,level) = 0.
  select case (bnd_type)
    case (MG_BND_PERIODIC) ! periodic
      do i = 0, nx-1
        Gx(i,level) = dxinv * (cos(2.*i*nxinv) - 1.)
      enddo
      do j = 0, ny-1
        Gy(j,level) = dyinv * (cos(2.*j*nyinv) - 1.)
      enddo
      do k = 0, nz-1
        Gz(k,level) = dzinv * (cos(2.*k*nzinv) - 1.)
      enddo
      norm(level) = 2./(float(nx)*float(ny)*float(nz))
    case (MG_BND_DIRICHLET) ! Dirichlet
      do i = 1, nx
        Gx(i,level) = dxinv * (cos(i*nxinv) - 1.)
      enddo
      do j = 1, ny
        Gy(j,level) = dyinv * (cos(j*nyinv) - 1.)
      enddo
      do k = 1, nz
        Gz(k,level) = dzinv * (cos(k*nzinv) - 1.)
      enddo
      norm(level) = 2.**ndim / (float(nx)*float(ny)*float(nz))
    case (MG_BND_NEUMANN) ! Neumann
      do i = 0, nx-1
        Gx(i,level) = dxinv * (cos(i*nxinv) - 1.)
      enddo
      do j = 0, ny-1
        Gy(j,level) = dyinv * (cos(j*nyinv) - 1.)
      enddo
      do k = 0, nz-1
        Gz(k,level) = dzinv * (cos(k*nzinv) - 1.)
      enddo
      norm(level) = 2.**ndim / (float(nx)*float(ny)*float(nz))
  end select
  if (ndim == 2) then
    do j = 0, ny
      do i = 0, nx
        if (abs(Gx(i,level)+Gy(j,level)) > 1.D-99) then
          G2D(i,j,level) = 1./(Gx(i,level)+Gy(j,level))
        else
          G2D(i,j,level) = 0.
        endif
      enddo
    enddo
  elseif (ndim == 3) then
    do k = 0, nz
      do j = 0, ny
        do i = 0, nx
          if (abs(Gx(i,level)+Gy(j,level)+Gz(k,level)) > 1.D-99) then
            G3D(i,j,k,level) = 1./(Gx(i,level)+Gy(j,level)+Gz(k,level))
          else
            G3D(i,j,k,level) = 0.
          endif
        enddo
      enddo
    enddo
  endif
  ip(0,level) = 0
  lastBndType = bnd_type
  iinit(level) = 0
endif

if (ndim == 1) then
  call gr_hgPoissonSolve1D (soln, nx, dx, bnd_type, &
                      Gx(0,level), wk1(0,level), &
                      ip(0,level), norm(level), nwk1, nip)
elseif (ndim == 2) then
  call gr_hgPoissonSolve2D (soln, nx, ny, dx, dy, bnd_type, &
                      Gx(0,level), Gy(0,level), wk1(0,level), &
                      wk2(0,level), ip(0,level), norm(level), &
                      nwk1, nwk2, nip, G2D(0,0,level))
else
  call gr_hgPoissonSolve3D (soln, nx, ny, nz, dx, dy, dz, bnd_type, &
                      Gx(0,level), Gy(0,level), Gz(0,level), &
                      wk1(0,level), wk2(0,level), ip(0,level), norm(level), &
                      nwk1, nwk2, nip, G3D(0,0,0,level))
endif

!===============================================================================

return
end subroutine gr_hgPoissonSolveBlock
