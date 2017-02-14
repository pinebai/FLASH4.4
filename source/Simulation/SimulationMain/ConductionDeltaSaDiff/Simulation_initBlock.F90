!!****if* source/Simulation/SimulationMain/ConductionDeltaSaDiff/Simulation_initBlock
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
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up a linear conduction problem
!!  to test conduction in a medium with constant isochoric conduction coefficient.
!!  The exact three-dimensional solution for the initial condition
!!  
!!                 T (r, t = 0) = Q delta (0)
!!
!!  is 
!!
!!       T (r, t) = Q / (4 pi \chi t)^(3/2) exp [-r^2 / 4 \chi t]
!!
!!  Here we set up the initial condition with the exact solution
!!  slightly offset from t = 0.
!!
!!  Reference: 
!!
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data
  use Conductivity_interface, ONLY : Conductivity
  use Eos_interface, ONLY : Eos
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData
  use Driver_interface, ONLY: Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"   

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  real :: xx, yy, zz, x0, y0, z0, r2, toff, pi
  real :: nexpo, xi, fOfXi
  integer :: ivelx, ively, ivelz

  real,allocatable,dimension(:) :: xCenter,yCenter,zCenter

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ,vecLen
  integer, dimension(MDIM) :: axis

  real, dimension(:), allocatable :: eosData
  
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       enerZone, ekinZone, eintZone, tempZone
  real :: cond_zone, diff_coeff, cond_chi
  real, parameter :: bogusUnusedTemperature = 3.e-5
  real, save :: bogusUnusedMassFrac(NSPECIES)
  
  logical :: gcell = .true.
  
  
  
  pi = 4. * atan (1.d0)
  
  ! WE COULD dump some output to stdout listing the parameters.
  ! But currently we don't.
  if (sim_meshMe == MASTER_PE) then
     
1    format (1X, 1P, 4(A7, E13.7, :, 1X))
2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
     
  endif
  
  
  !get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  
  allocate(xCenter(sizeX))
      allocate(yCenter(sizeY)) 
  allocate(zCenter(sizeZ))
  
  vecLen = 1
  
  allocate(eosData(vecLen * EOS_NUM))
  
  if (NDIM == 3) call Grid_getCellCoords&
       (KAXIS, blockId, CENTER,gcell, zCenter, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
       (JAXIS, blockId, CENTER,gcell, yCenter, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  !------------------------------------------------------------------------------
  
  !Coordinates of initial point source
  x0 = sim_xCenter
  y0 = sim_yCenter
  z0 = sim_zCenter
  
  ! Loop over cells in the block.  For each, compute the physical position of 
  !  the state using an exact solution.
  
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)     
     ! get the coordinates of the cell center in the z-direction
     zz = zCenter(k)     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)        
        ! get the coordinates of the cell center in the y-direction
        yy = yCenter(j)        
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           xx  = xCenter(i)

           ! In planar geometry, set coordinate equal to distance along
           !  a cardinal axis
           if (sim_orientation .eq. 1) then
              r2 = (xx - x0) * (xx - x0)              
           else if (sim_orientation .eq. 2) then
              r2 = (yy - y0) * (yy - y0)              
           else if (sim_orientation .eq. 3) then
              r2 = (zz - z0) * (zz - z0)
              ! Otherwise, in 3-D spherical geometry, use the 3d radial distance,              
           else if (sim_orientation .eq. 0) then
              if (NDIM .EQ. 1) then
!!$                 call Driver_abortFlash("Conduction point source in 3D and 2D only.")
                 r2 =  (xx - x0) * (xx - x0)
              else if (NDIM .EQ. 2) then
                 r2 = (xx - x0) * (xx - x0) + (yy - y0) * (yy - y0)
              else
                 r2 = (xx - x0) * (xx - x0) + (yy - y0) * (yy - y0) + &
                      (zz - z0) * (zz - z0)
              endif
           end if

           call Conductivity(bogusUnusedTemperature, sim_rhoInit, bogusUnusedMassFrac, cond_zone, diff_coeff, 2)
           if (diff_coeff .NE. 0.0) then
              cond_chi = diff_coeff
           else
              cond_chi = 2.0E9 / (3.5 * 8.3145119843E7)
           end if         

           ! Now, initialize temperature solution according to the exact solution
           !       T (r, t) = Q / (4 pi \chi t)^(3/2) exp [-r^2 / 4 \chi t]
           ! for a point source in 3D, and
           !       T (x, t) = Q / (4 pi \chi t)^(1/2) exp [-x^2 / 4 \chi t]
           ! for a planar source in any dimension. 
           !
           ! Here \chi is the thermal diffusivity obtained in the Conduction subroutine.
           ! (In the case of constant conductivity, \chi is computed from the value stored
           ! in cond_constantIsochoric.)
           ! The exact solution at t = 0 is a Dirac delta function, so we initialize
           !   with some offset in time toff. The gas is initially stationary.
           
           toff = sim_toffset
           rhoZone = sim_rhoInit

           if (sim_orientation .eq. 0) then              
              ! Initialize point source 
              if (NDIM == 1) then 
                 tempZone = sim_Q / (4. * pi * cond_chi * toff)**(3./2.) * &
                      exp (-r2 / (4. * cond_chi * toff) )
              else if (NDIM == 2) then 
                 tempZone = sim_Q / (4. * pi * cond_chi * toff)* &
                      exp (-r2 / (4. * cond_chi * toff) )
              else
                 tempZone = sim_Q / (4. * pi * cond_chi * toff)**(3./2.) * &
                      exp (-r2 / (4. * cond_chi * toff) )
              endif
              
           else
              ! Initialize planar source              
              if (sim_iniCondTemperatureExponent==0) then
                 tempZone = sim_Q / (4. * pi * cond_chi * toff)**(1./2.) * &
                      exp (-r2 / (4. * cond_chi * toff) )
              else
                 nexpo = sim_iniCondTemperatureExponent
                 call sim_xToXi(sqrt(r2),0.0,nexpo,xi)
                 if (xi .GE. sim_xi0) then
                    tempZone = 0.0
                 else
999                 format('sqrt(r2),xi,xi/x0,f=',4(1PG15.8))
                 fOfXi = ( nexpo*sim_xi0**2/(2*(nexpo+2)) * (1-(xi/sim_xi0)**2) )**(1/nexpo)
                 print 999,sqrt(r2),xi,xi/sim_xi0,fOfXi
                 tempZone = (sim_Q**2 / (sim_alpha * sim_toffset) )**(1.0/(nexpo+2)) * fOfXi
              end if
           end if
        endif
        
        ! Now add this Guassian shape to the temperature background
        tempZone = sim_tempBackground + tempZone
        
        velxZone = 0.
        velyZone = 0.
        velzZone = 0.
        
#define DO_CALL_EOS
        
#ifdef DO_CALL_EOS
        ! Get thermodynamic quantities and load them on the grid using
        ! the EOS routine        
        vecLen = 1
        
        eosData(EOS_TEMP) = tempZone
        eosData(EOS_DENS) = rhoZone
        
        CALL Eos(MODE_DENS_TEMP, vecLen, eosData)
        
        if (tempZone .NE. eosData(EOS_TEMP)) then
           print*,'Something went wrong in Simulation_initBlock, tempZone=',tempZone,' but eosData(EOS_TEMP)=',eosData(EOS_TEMP)
           tempZone = eosData(EOS_TEMP)
        end if
        if (rhoZone .NE. eosData(EOS_DENS)) then
           print*,'Something went wrong in Simulation_initBlock, rhoZone=',rhoZone,' but eosData(EOS_DENS)=',eosData(EOS_DENS)
           rhoZone = eosData(EOS_DENS)
        end if
        presZone = eosData(EOS_PRES)
        
        if (presZone .lt. sim_smallP) then
!!$          print *, 'Error : pressure < smallP.'
!!$           call Driver_abortFlash('Error : pressure < smallP.')
        endif
        
        eintZone = eosData(EOS_EINT)
        
        !       block_data(EINT_VAR, loop_i, loop_j, loop_k) = eosData(EOS_EINT)
        !       block_data(GAMC_VAR, loop_i, loop_j, loop_k) = eosData(EOS_GAMC)
        !       block_data(GAME_VAR, loop_i, loop_j, loop_k) = eosData(EOS_PRES)/(eosData(EOS_EINT)*eosData(EOS_DENS)) + 1.0E 
#else
        
        presZone = rhoZone * tempZone
        
        eintZone = presZone / (sim_gamma-1.)
        eintZone = eintZone / rhoZone
        
#endif
        
        ! Compute the gas energy and set the gamma-values needed for the equation of 
        ! state.
        ekinZone = 0.5 * (velxZone**2 + & 
             velyZone**2 + & 
             velzZone**2)
        
        enerZone = eintZone + ekinZone
        enerZone = max(enerZone, sim_smallP)
        
        axis(IAXIS) = i
        axis(JAXIS) = j
        axis(KAXIS) = k
        
#if NSPECIES > 0
        call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
             axis, 1.0e0-(NSPECIES-12)*sim_smallX)
        
        !if there is only 1 species, this loop will not execute
        do n = SPECIES_BEGIN+1,SPECIES_END
           call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                axis, sim_smallX)
        enddo
#endif
        
        ! store the variables in the current zone via the database put methods
        ! data is put stored one cell at a time with this call to Grid_putData           
        
        call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
        call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
        
        call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
        call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone) 
        call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)

#ifdef EINT_VAR
        call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eintZone)   
#endif
#ifdef ENER_VAR
        call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)   
#endif
#ifdef GAME_VAR          
        call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
#endif
#ifdef GAMC_VAR
        call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
#endif
        
        
     enddo
  enddo
enddo
deallocate(eosData)
deallocate(xCenter)
deallocate(yCenter)
deallocate(zCenter)

return

end subroutine Simulation_initBlock
