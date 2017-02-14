!!****if* source/Simulation/SimulationMain/PoisTest/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Huang & Greengard
!!  Poisson solver test problem (their example 4.3).
!!
!!  Reference:   Huang, J., & Greengard, L. 2000, SIAM J. Sci. Comput.,
!!               21, 1551
!! 
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!
!! PARAMETERS
!!
!! sim_smlRho
!! sim_gam
!!
!!***

subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY: sim_smlRho, sim_gam
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putPointData

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockID
  

  real,allocatable, dimension(:) :: xl, xr, yl, yr, zl, zr
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX, sizeY, sizeZ
  integer, dimension(MDIM) :: axis

  !CD: Coordinates for guardcells are returned along with the interior cells.
  logical :: gcell = .true.
  integer :: istat 

  integer       :: i, j, k, n
  integer       :: Nint, ii, jj, kk, h
  real          :: xx, yy, zz
  real          :: xdist, ydist, zdist, dist
  real          :: Nintinv, sum_rho, exponent, Nintinv1, rho

  integer, parameter :: Npeak = 13
  real :: xctr(Npeak) = & 
       &          (/ 0., -1., -1., 0.28125, 0.5, 0.3046875, 0.3046875, & 
       &             0.375, 0.5625, -0.5, -0.125, 0.296875, 0.5234375 /)
  real :: yctr(Npeak) = & 
       &          (/ 0., 0.09375, 1., 0.53125, 0.53125, 0.1875, 0.125, & 
       &             0.15625, -0.125, -0.703125, -0.703125, -0.609375, & 
       &             -0.78125 /)
  real :: zctr(Npeak) = & 
       &          (/0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./)
  real :: sigma(Npeak) = & 
       &          (/ 0.01, 4000., 20000., 80000., 16., 360000., 400000., & 
       &             2000., 18200., 128., 49000., 37000., 18900. /)


  ! dump some output to stdout listing the paramters
!!$  if (sim_meshMe == MASTER_PE) then
!!$
!!$     !print *, "Initialising block:", blockID
!!$     !write (*,*) 'using ', Npeak, ' density peaks:'
!!$     !write (*,*)
!!$     !write (*,3) 'Peak', 'x', 'y', 'z', 'sigma'
!!$     !do i = 1, Npeak
!!$     !   write (*,4) i, xctr(i), yctr(i), zctr(i), sigma(i)
!!$     !enddo
!!$     !write (*,*)
!!$
!!$
!!$1    format (1X, 1P, 4(A13, E12.6, :, 1X))
!!$2    format (1X, 1P, A13, I12)
!!$3    format (1X, A5, 2X, 4(A13, :, 2X))
!!$4    format (1X, I5, 2X, 1P, 4(E13.6, :, 2X))
!!$
!!$  endif


  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  !CD: blkLimitsGC is the block limits from 1 to NXB+2*NGUARD (i.e. inc guardcells).
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xl(sizeX),stat=istat)
  allocate(xr(sizeX),stat=istat)
  allocate(yl(sizeY),stat=istat)
  allocate(yr(sizeY),stat=istat)
  allocate(zl(sizeZ),stat=istat)
  allocate(zr(sizeZ),stat=istat)

  xl = 0.0
  xr = 0.0
  yl = 0.0
  yr = 0.0
  zl = 0.0
  zr = 0.0

  if (NDIM == 3) then
     call Grid_getCellCoords(KAXIS, blockId, LEFT_EDGE, gcell, zl, sizeZ)
     call Grid_getCellCoords(KAXIS, blockId, RIGHT_EDGE, gcell, zr, sizeZ)
  end if

  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE, gcell, yl, sizeY)
     call Grid_getCellCoords(JAXIS, blockId, RIGHT_EDGE, gcell, yr, sizeY)
  end if

  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xl, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xr, sizeX)


  Nint    = 7
  Nintinv = 1./ real(Nint)
  Nintinv1= 1./ (real(Nint)-1.)

  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)    
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)       

           sum_rho = 0.

           do kk = 0, (Nint-1)*K3D
              zz    = zl(k) + kk*Nintinv1*(zr(k)-zl(k))
              do jj = 0, (Nint-1)*K2D
                 yy    = yl(j) + jj*Nintinv1*(yr(j)-yl(j))
                 do ii = 0, Nint-1
                    xx    = xl(i) + ii*Nintinv1*(xr(i)-xl(i))

                    do h = 1, Npeak
                       xdist = (xx - xctr(h))**2
                       ydist = (yy - yctr(h))**2 * K2D
                       zdist = (zz - zctr(h))**2 * K3D
                       dist  = xdist + ydist + zdist
                       exponent = sigma(h) * dist

                       !Avoid underflow by only computing exponential
                       !for "small" exponent magnitudes
                       if (exponent < 80.) then 
                          sum_rho = sum_rho + exp(-exponent)
                       end if

                    enddo

                 enddo
              enddo
           enddo

           rho = max(sum_rho * Nintinv**NDIM, sim_smlRho)

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           ! store the variables in the current zone via the database put methods.
           ! data is put stored one cell at a time with this call to Grid_putPointData.
           call Grid_putPointData(blockID, CENTER, GAME_VAR, EXTERIOR, axis, sim_gam)
           call Grid_putPointData(blockID, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gam)
           call Grid_putPointData(blockID, CENTER, DENS_VAR, EXTERIOR, axis, rho)

        enddo
     enddo
  enddo


  !! Cleanup!  Must deallocate arrays
  deallocate(xl)
  deallocate(xr)
  deallocate(yl)
  deallocate(yr)
  deallocate(zl)
  deallocate(zr)

!!$  if (sim_meshMe == MASTER_PE) then
!!$     print *, "Finished initialising block:", blockID
!!$  end if

  return
end subroutine Simulation_initBlock
