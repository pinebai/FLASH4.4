!!****if* source/Simulation/SimulationMain/MacLaurin/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Maclaurin spheroid problem.
!!
!!  References:  Maclaurin, C. 1742, 
!!               Chandrasekhar, S. 1987, Ellipsoidal Figures of Equilibrium
!!
!! ARGUMENTS
!!
!!  blockID     --   The number of the block to initialize
!!
!!***

subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY: sim_density, sim_gamma, sim_Omega2, sim_Pconst, &
       sim_nsubinv, sim_nsubzones, sim_xctr, sim_yctr, sim_zctr, &
       sim_initGeometry, sim_geom2DAxisymmetric, sim_geom3DCartesian, &
       sim_a3inv, sim_a1inv, &
       sim_smallRho, sim_smallP, sim_smallE
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putRowData

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in)  :: blockID

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(3) :: startingPos
  integer  :: sizeX, sizeY, sizeZ
  logical  :: gcell=.true.
  integer  :: i, j, k, ii, jj, kk
  real     :: xdist, ydist, zdist, dist2, rxy, rxyz2, rinv
  real     :: xx, yy, zz, dxx, dyy, dzz, vxfac, vyfac, vzfac
  real     :: sum_rho, sum_p, sum_vx, sum_vy, sum_vz, pres, vel

  real, dimension(:), allocatable :: x, y, z, xl, yl, zl, xr, yr, zr, dx, dy, dz
  real, dimension(:), allocatable :: vx, vy, vz, p, rho, e, ek, ei, gam


  ! Get the coordinate information for the current block

  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
  allocate(x(sizeX), xl(sizex), xr(sizeX), dx(sizeX))
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
  allocate(y(sizeY), yl(sizeY), yr(sizeY), dy(sizeY))
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1
  allocate(z(sizeZ), zl(sizeZ), zr(sizeZ), dz(sizeZ))
  if (NDIM == 3) then 
     call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, z, sizeZ)
     call Grid_getCellCoords(KAXIS, blockId, LEFT_EDGE, gcell, zl, sizeZ)
     call Grid_getCellCoords(KAXIS, blockId, RIGHT_EDGE, gcell, zr, sizeZ)
  endif
  if (NDIM >= 2) then 
     call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, y, sizeY)
     call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE, gcell, yl, sizeY)
     call Grid_getCellCoords(JAXIS, blockId, RIGHT_EDGE, gcell, yr, sizeY)
  endif
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, x, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xl, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xr, sizeX)

  dx(:) = xr(:) - xl(:)
  dy(:) = yr(:) - yl(:)
  dz(:) = zr(:) - zl(:)

  ! Set initial conditions in each zone

  allocate(vx(sizeX), vy(sizeX), vz(sizeX), p(sizeX), rho(sizeX), e(sizeX), & 
       ek(sizeX), ei(sizeX), gam(sizeX))
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     dzz = dz(k) * sim_nsubinv
     do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
        dyy = dy(j) * sim_nsubinv
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
           dxx = dx(i) * sim_nsubinv

           sum_rho = 0.0
           sum_p   = 0.0
           sum_vx  = 0.0
           sum_vy  = 0.0
           sum_vz  = 0.0

           ! Break the zone into nsubzones^ndim sub-zones and average the
           ! results to get the values for the zone.

           do kk = 0, sim_nsubzones-1
              zz    = zl(k) + (kk + 0.5)*dzz 
              zdist = (zz - sim_zctr) * K3D

              do jj = 0, sim_nsubzones-1
                 yy    = yl(j) + (jj + 0.5)*dyy
                 ydist = (yy - sim_yctr) * K2D

                 do ii = 0, sim_nsubzones-1
                    xx    = xl(i) + (ii + 0.5)*dxx
                    xdist = xx - sim_xctr

                    select case (sim_initGeometry)

                    case (sim_geom2DAxisymmetric)

                       dist2 = (xdist*sim_a1inv)**2 + (ydist*sim_a3inv)**2
                       rxy   = xdist
                       rxyz2 = dist2
                       rinv  = 1./sqrt(xdist**2 + ydist**2)
                       vxfac = 0.0
                       vyfac = 0.0
                       vzfac = 0.0

                    case (sim_geom3DCartesian)

                       dist2 = (xdist*sim_a1inv)**2 + (ydist*sim_a1inv)**2 + (zdist*sim_a3inv)**2
                       rxy   = sqrt(xdist**2 + ydist**2)
                       rxyz2 = (rxy*sim_a1inv)**2 + (zdist*sim_a3inv)**2
                       rinv  = 1./sqrt(xdist**2 + ydist**2 + zdist**2)
                       vxfac = -ydist*rinv
                       vyfac = xdist*rinv
                       vzfac = 0.0

                    end select

                    if (dist2 <= 1.) then    ! inside the spheroid
                       pres    = sim_Pconst * (1.0 - rxyz2)
                       vel     = rxy * sim_Omega2
                       sum_rho = sum_rho + sim_density
                       sum_p   = sum_p   + pres
                       sum_vx  = sum_vx  + vel*vxfac
                       sum_vy  = sum_vy  + vel*vyfac
                       sum_vz  = sum_vz  + vel*vzfac
                    else                     ! outside the spheroid
                       sum_rho = sum_rho + sim_smallRho*10
                       sum_p   = sum_p   + sim_smallP*10
                    endif

                 enddo
              enddo
           enddo

           rho(i) = max(sum_rho * sim_nsubinv**3, sim_smallRho)
           p(i)   = max(sum_p * sim_nsubinv**3, sim_smallP)
           vx(i)  = sum_vx*sim_nsubinv**3
           vy(i)  = sum_vy*sim_nsubinv**3
           vz(i)  = sum_vz*sim_nsubinv**3
           ek(i)  = 0.5*(vx(i)**2+vy(i)**2+vz(i)**2)
           gam(i) = sim_gamma
           ei(i)  = max(p(i)/(rho(i)*(gam(i)-1.0)), sim_smallE)
           e(i)   = ei(i) + ek(i)

        enddo

        startingPos(1) = 1
        startingPos(2) = j
        startingPos(3) = k
        call Grid_putRowData(blockID, CENTER, DENS_VAR, EXTERIOR, IAXIS, startingPos, rho, sizeX)
        call Grid_putRowData(blockID, CENTER, PRES_VAR, EXTERIOR, IAXIS, startingPos, p, sizeX)
        call Grid_putRowData(blockID, CENTER, ENER_VAR, EXTERIOR, IAXIS, startingPos, e, sizeX)
        call Grid_putRowData(blockID, CENTER, GAME_VAR, EXTERIOR, IAXIS, startingPos, gam, sizeX)
        call Grid_putRowData(blockID, CENTER, GAMC_VAR, EXTERIOR, IAXIS, startingPos, gam, sizeX)
        call Grid_putRowData(blockID, CENTER, VELX_VAR, EXTERIOR, IAXIS, startingPos, vx, sizeX)
        call Grid_putRowData(blockID, CENTER, VELY_VAR, EXTERIOR, IAXIS, startingPos, vy, sizeX)
        call Grid_putRowData(blockID, CENTER, VELZ_VAR, EXTERIOR, IAXIS, startingPos, vz, sizeX)
        call Grid_putRowData(blockID, CENTER, EINT_VAR, EXTERIOR, IAXIS, startingPos, ei, sizeX)

     enddo
  enddo
  deallocate(rho, p, vx, vy, vz, e, ei, ek, gam)
  deallocate(x, xl, xr, dx)
  deallocate(y, yl, yr, dy)
  deallocate(z, zl, zr, dz)

  return
end subroutine Simulation_initBlock

