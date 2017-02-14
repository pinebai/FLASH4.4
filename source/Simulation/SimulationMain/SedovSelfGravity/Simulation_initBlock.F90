!!****if* source/Simulation/SimulationMain/SedovSelfGravity/Simulation_initBlock
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
!!  a specified block.  This version sets up the Sedov spherical
!!  explosion problem with self gravity
!!
!!  References:  Sedov, L. I., 1959, Similarity and Dimensional Methods
!!                 in Mechanics (New York:  Academic)
!!
!!               Landau, L. D., & Lifshitz, E. M., 1987, Fluid Mechanics,
!!                 2d ed. (Oxford:  Pergamon)
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  
!!
!!***
subroutine Simulation_initBlock(blockID)
  use Simulation_data, ONLY : sim_nsubzones, sim_inSubzm1, sim_inSubZones,&
                              sim_rProf,sim_vProf, sim_pProf, sim_rhoProf, &
                              sim_gamma, SIM_NPROFILE,&
                              sim_smlrho, sim_smallp, sim_xn
  use Grid_interface, ONLY : Grid_getBlkIndexLimits,Grid_putPointData,&
                             Grid_getCellCoords
  
  implicit none
#include "constants.h"
#include "Flash.h"

  integer,intent(IN)  ::  blockID
  

  integer  ::  i, jlo, jhi
  integer  ::  ii,n

  real     ::  sum_rho, sum_p, sum_vx, vel
  real     ::  vx, pres, rho, eint, ek,gam


  real,allocatable, DIMENSION(:) :: xCenter, xLeft, xRight,blkDensity

  real :: xx, frac
  integer , dimension(MDIM) :: pos
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: iSize

  pos = 1

  call Grid_getBlkIndexLimits(blockID,blkLimits, blkLimitsGC)

  iSize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) +1
  allocate(xLeft(iSize))
  allocate(xRight(iSize))
  allocate(xCenter(iSize))
  allocate(blkDensity(iSize))

  call Grid_getCellCoords(IAXIS,blockID,LEFT_EDGE,.true.,xLeft,iSize)
  call Grid_getCellCoords(IAXIS,blockID,RIGHT_EDGE,.true.,xRight,iSize)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,xCenter,iSize)



!
! For each cell, loop over phi (3-d only), theta (2- and 3-d), and r.
! We will radially subzone in r.  It doesn't really matter in phi or 
! theta, since this is a purely radial problem.   
!  

  do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
     pos(IAXIS)=i
     sum_rho = 0.
     sum_p   = 0.
     sum_vx  = 0.
     
     !       Break the cell into nsubzones sub-zones, and look up the
     !       appropriate quantities along the 1d profile for that subzone.  
     !
     !       Have the final values for the zone be equal to the average of
     !       the subzone values.
     ! 
     do ii = 0, (sim_nsubzones-1)
        xx    = xCenter(i) + (ii*sim_inSubzm1-.5)*(xRight(i)-xLeft(i))
        
        call find (sim_rProf, SIM_NPROFILE, xx, jlo)
        
        !
        !  a point at `dist' is frac-way between jlo and jhi.   We do a
        !  linear interpolation of the quantities at jlo and jhi and sum those.
        ! 
        
        if (jlo .eq. 0) then
           jlo = 1
           jhi = 1
           frac = 0.
           
        else if (jlo .eq. SIM_NPROFILE) then
           jlo = SIM_NPROFILE
           jhi = SIM_NPROFILE
           frac = 0.
           
        else
           jhi = jlo + 1
           frac = (xx - sim_rProf(jlo)) / (sim_rProf(jhi)-sim_rProf(jlo))
        endif
        
        ! 
        !   Now total these quantities.   Note that  v is a radial velocity; 
        !   we multiply by the tangents of the appropriate angles to get
        !   the projections in the x, y and z directions.
        !
        sum_p = sum_p + sim_pProf(jlo) + &
             frac*(sim_pProf(jhi)  - sim_pProf(jlo))
        sum_rho = sum_rho + sim_rhoProf(jlo) + &
             frac*(sim_rhoProf(jhi)- sim_rhoProf(jlo))
        vel = sim_vProf(jlo) + frac*(sim_vProf(jhi)  - sim_vProf(jlo))
        
        sum_vx  = sum_vx  + vel
        
     enddo
     
     pres = max(sum_p   * sim_inSubzones, sim_smallp)
     rho  = max(sum_rho * sim_inSubzones, sim_smlrho)
     vx   = sum_vx  * sim_inSubzones
     
     ek  = 0.5*(vx*vx)

!
!  assume gamma-law equation of state
!
     gam = sim_gamma
     eint = pres/(gam-1.)
     eint = eint/rho + ek
     eint = max (eint, sim_smallp)
           
     do n=SPECIES_BEGIN,SPECIES_END
        call Grid_putPointData(blockID,CENTER,n,EXTERIOR,pos,&
            sim_xn(n-SPECIES_BEGIN+1))
     enddo
     
     call Grid_putPointData(blockID,CENTER,DENS_VAR,EXTERIOR,pos,rho)
     call Grid_putPointData(blockID,CENTER,PRES_VAR,EXTERIOR,pos,pres)
     call Grid_putPointData(blockID,CENTER,VELX_VAR,EXTERIOR,pos,vx)
     call Grid_putPointData(blockID,CENTER,GAME_VAR,EXTERIOR,pos,gam)
     call Grid_putPointData(blockID,CENTER,GAMC_VAR,EXTERIOR,pos,gam)
     call Grid_putPointData(blockID,CENTER,ENER_VAR,EXTERIOR,pos,eint)
  end do
  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(blkDensity)  
  return
end subroutine Simulation_initBlock
