!!****if* source/Simulation/SimulationMain/GrayDiffRadShock/Simulation_initBlock
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
!!  a specified block.
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  
!!
!!
!!***



subroutine Simulation_initBlock(blockId)
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putPointData, Grid_getPointData
  use Driver_interface, ONLY: Driver_abortFlash, Driver_getMype
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT
  use Eos_interface, ONLY : Eos_wrapped

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"    

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer, intent(in) :: blockId
  
  integer :: i, j, k
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)
  real, allocatable :: xcent(:), ycent(:), zcent(:)

  real :: rho   ! Density
  real :: tion  ! Ion temperature
  real :: tele  ! Electron temperature 
  real :: trad  ! Radiation temperature
  real :: tradActual

  
  real :: gasConst = 8.31447E+07
  real :: gamma    = 5.0/3.0
  real :: dens_ref = 1.0
  real :: clight   = 2.99792E+10  
  real :: L_ref  = 1.0
  real :: radConst = 7.56577E-15

  real :: abs_opac    = 1.0E6
  real :: kappa       = 1.0  

  real :: cs_ref, trans_opac_actual, abs_opac_actual, velx

  real :: rho1, t1

  integer :: meshMe
  
  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  
  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., &
       xcent, blkLimitsGC(HIGH, IAXIS))
  
  allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., &
       ycent, blkLimitsGC(HIGH, JAXIS))
  
  allocate(zcent(blkLimitsGC(HIGH, KAXIS)))
  call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., &
       zcent, blkLimitsGC(HIGH, KAXIS)) 
  
  
  cs_ref = sqrt(gamma*gasConst*sim_temp)    
  
  dens_ref = radConst*(sim_temp**4)/(sim_P0*cs_ref**2)
  
  ! sim_rho =  dens_ref
  
  trans_opac_actual = (1.0/kappa)*(clight/(3.0*cs_ref*L_ref))
  
  abs_opac_actual   = abs_opac*(cs_ref/(L_ref*clight))      
  

  if (sim_M0 == 1.05) then

     t1   = 1.04945451752
     rho1 = 1.0749588726     

  else if (sim_M0 == 1.2) then

     t1   = 1.19475152105
     rho1 = 1.29732134522

  else if (sim_M0 == 1.4) then

     t1   = 1.39173649786
     rho1 = 1.5807111889
     
  else if (sim_M0 == 2.0) then
     
     if (sim_P0 == 1.0e-4) then
        t1   = 2.07756999533
        rho1 = 2.28607489893    
        print *, "HERE"
     else if (sim_P0 == 1.0e-1) then
        t1   = 1.80622128471
        rho1 = 2.42750848291
     endif
     
  else if (sim_M0 == 3.0) then
     
     t1   = 3.66191266581
     rho1 = 3.00216769711
     
  else if (sim_M0 == 5.0) then
     
     t1   = 8.55719921848
     rho1 = 3.59791065306
     
  end if


!!$  call Driver_getMype(MESH_COMM,meshMe)
!!$
!!$  if (meshMe == 0) then   
!!$     
!!$     write (*,*) "UPSTREAM CONDITIONS"
!!$     write (*,*) "DENS:", dens_ref
!!$     write (*,*) "TELE:", sim_temp
!!$     write (*,*) "TION:", sim_temp
!!$     write (*,*) "TRAD:", sim_temp
!!$     write (*,*) "VELX:", sim_M0*cs_ref  
!!$     
!!$     write (*,*) "DOWNSTREAM CONDITIONS"
!!$     write (*,*) "DENS:", rho1*dens_ref
!!$     write (*,*) "TELE:", t1*sim_temp
!!$     write (*,*) "TION:", t1*sim_temp
!!$     write (*,*) "TRAD:", t1*sim_temp
!!$     write (*,*) "VELX:", (sim_M0*cs_ref)/rho1
!!$     
!!$     write (*,*) "OPACITIES"
!!$     write (*,*) "Transport opacity :", trans_opac_actual
!!$     write (*,*) "Absorption opacity:", abs_opac_actual 
!!$     write (*,*) "Emission opacity  :", abs_opac_actual   
!!$
!!$  end if

  
  
  ! Loop over cells and set the initial state
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)          
           
           if (xcent(i) > 0.0) then !! downstream
              ! rho  = rho1*dens_ref
              rho = rho1
              tele = t1*sim_temp
              tion = t1*sim_temp
              trad = t1*sim_temp
              velx = (1.0/rho1)*sim_M0*cs_ref
           else !! upstream
              ! rho  = 1.0*dens_ref
              rho = sim_rho
              tele = 1.0*sim_temp
              tion = 1.0*sim_temp
              trad = 1.0*sim_temp
              velx = 1.0*sim_M0*cs_ref
           end if
           
           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k
           
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, tion)
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, tele)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velx)           
           
           call RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
           
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, tradActual)
           
           
        enddo
     enddo
  enddo
    
  call Eos_wrapped(MODE_DENS_TEMP_GATHER, blkLimits, blockId)

  
  !call Simulation_computeAnalytical(blockID, 0.0)
  
  
  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)
  
  return

end subroutine Simulation_initBlock
