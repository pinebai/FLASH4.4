!!****if* source/Simulation/SimulationMain/Pancake/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId, 
!!                       
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Zel'dovich pancake problem.
!!
!!  References:  Zel'dovich, Ya. B. 1970, AA, 5, 84
!!               Anninos, W. Y. and Norman, M. L. 1994, ApJ, 429, 434
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_xMin, sim_yMin, sim_zMin, &
       sim_xMax, sim_yMax, sim_zMax, sim_xcos, sim_ycos, sim_zcos, &
       sim_temp1, sim_temp2, sim_temp3, sim_kmag, sim_meanden, sim_bmeanden, &
       sim_zinitial, sim_Tfiducial, sim_zfiducial, sim_gascon, &
       sim_gamma, sim_smallE, sim_useParticles, sim_smallRho, sim_hubble, sim_zcaustic
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords, & 
                             Grid_putPointData
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) ::  blockId
  

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX, sizeY, sizeZ
  real, dimension(:), allocatable :: x, y, z
  logical :: gcell=.true.
  integer, dimension(MDIM) :: axis

  integer :: i, j, k, g, n
  real :: xe, xl, xlold, err, temp, vel, abarinv
  real :: rhodm, rho, p, t, e, ek, gam, vx, vy, vz
  real, dimension(SPECIES_END) :: xn


! Get the coordinate information for the current block

call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
allocate(x(sizeX))
sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
allocate(y(sizeY))
sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1
allocate(z(sizeZ))
call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, z, sizeZ)
call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, y, sizeY)
call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, x, sizeX)

! Initialize the density, pressure, and velocity fields.

do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
  do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
    do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

! Calculate Lagrangian coordinate

      xe = (x(i) - 0.5*(sim_xMin+sim_xMax)) * sim_xcos + &
           (y(j) - 0.5*(sim_yMin+sim_yMax)) * sim_ycos + &
           (z(k) - 0.5*(sim_zMin+sim_zMax)) * sim_zcos   
      xlold = xe
      g = 0
      err = huge(1.0)
      do while ((err > 1.E-6*abs(xlold)) .and. (g <= 500))  
        xl = xe + sim_temp2 * sin(sim_kmag*xlold)
        err = abs(xl-xlold)
        xlold = xl
        g = g + 1
      enddo
      if (g > 500) then
        print *, i, j, k
        print *, xe, xl, err
        call Driver_abortFlash("Simulation_initBlock:  Lagrange iteration failed to converge!")
      endif

! Calculate analytical solution
      
      rhodm = (sim_meanden - sim_bmeanden) / (1. - sim_temp1*cos(sim_kmag*xl))
      
#ifdef PRES_VAR
      rho   = sim_bmeanden / (1. - sim_temp1*cos(sim_kmag*xl))
      vel   = -sim_hubble * (1.+sim_zcaustic) * sqrt(1.+sim_zinitial) * sin(sim_kmag*xl) / &
                          sim_kmag
      vx    = vel * sim_xcos
      vy    = vel * sim_ycos
      vz    = vel * sim_zcos

      temp  = sim_Tfiducial * (1.+sim_zinitial)**2 * &
              ( ((1.+sim_zfiducial)/(1.+sim_zinitial))**3 * &  
                (1.-sim_temp1*cos(sim_kmag*xl))/(1.-sim_temp3*cos(sim_kmag*xl)) ) ** (1.-sim_gamma)

!      call query_mfluid_suminv (mf_prop_A, xn(i,:), abarinv)
      abarinv = 0.75+0.25/4.0
      p     = rho * sim_gascon * temp * abarinv

! If OmegaBaryon is zero, render the gas harmless

      if (rho < sim_smallRho) then
        rho  = sim_smallRho
        temp = sim_Tfiducial * (1.+sim_zinitial)**2 * &
               ( ((1.+sim_zfiducial)/(1.+sim_zinitial))**3 ) ** (1.-sim_gamma)
        p    = rho * sim_gascon * temp * abarinv
        vx   = 0.
        vy   = 0.
        vz   = 0.
      endif

! Compute total energy

      ek = 0.5 * (vx**2 + vy**2 + vz**2)
      e  = p/(rho*(sim_gamma-1.0))
      e  = e + ek
      e  = max(e, sim_smallE)

      gam = sim_gamma

#endif
      
! Put initial data into database for this zone

      axis(IAXIS) = i
      axis(JAXIS) = j
      axis(KAXIS) = k

#ifdef PRES_VAR
      call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
      call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
      call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
      call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
      call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
      call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, gam)
      call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, gam)
      call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)
#endif

      if (NSPECIES == 1) then 
        xn(1) = 1.0
        call Grid_putPointData(blockId, CENTER, SPECIES_BEGIN, EXTERIOR, axis, xn(1))
      else if (NSPECIES > 1) then
        xn(:) = 1.0E-10
        xn(1) = 0.75
        xn(2) = 0.25
        do n = SPECIES_BEGIN, SPECIES_END
          call Grid_putPointData(blockId, CENTER, n, EXTERIOR, axis, xn(n))
        enddo
      endif

! If particles are included, set up their mesh density (in case this is
! used in the initial refinement iteration)

      if (sim_useParticles) then
        rhodm = max(rhodm, sim_smallRho)
        call Grid_putPointData(blockId, CENTER, PDEN_VAR, EXTERIOR, axis, rhodm)
      endif

    enddo
  enddo
enddo

deallocate(x)
deallocate(y)
deallocate(z)
return
end subroutine Simulation_initBlock
