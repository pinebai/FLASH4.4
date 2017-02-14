!!****if* source/Simulation/SimulationMain/Blast2/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!
!! ARGUMENTS
!!  
!!   none
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data 
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  implicit none
#include "constants.h"
#include "Flash.h"

  call Driver_getMype(MESH_COMM, sim_meshMe)
  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smalle', sim_smallE)
  call RuntimeParameters_get('smlrho', sim_smallRho)
  
  call RuntimeParameters_get('gamma', sim_gamma)
  
  call PhysicalConstants_get("electron mass",sim_eMassInUAmu,unitMass="amu")

  call RuntimeParameters_get('sim_rhoLeft', sim_rhoLeft)
  call RuntimeParameters_get('sim_rhoRight', sim_rhoRight)
  call RuntimeParameters_get('sim_rhoMid', sim_rhoMid)
  
  call RuntimeParameters_get('sim_pLeft', sim_pLeft)
  call RuntimeParameters_get('sim_pRight', sim_pRight)
  call RuntimeParameters_get('sim_pMid', sim_pMid)
  
  call RuntimeParameters_get('sim_uLeft', sim_uLeft)
  call RuntimeParameters_get('sim_uRight', sim_uRight)
  call RuntimeParameters_get('sim_uMid', sim_uMid)
  
  call RuntimeParameters_get('sim_xangle', sim_xAngle)
  call RuntimeParameters_get('sim_yangle', sim_yAngle)
  
  call RuntimeParameters_get('sim_posnL', sim_posnL)
  call RuntimeParameters_get('sim_posnR', sim_posnR)


  ! convert the shock angle paramters
  sim_xAngle = sim_xAngle * 0.0174532925 ! Convert to radians.
  sim_yAngle = sim_yAngle * 0.0174532925

  sim_xCos = cos(sim_xAngle)
  
  if (NDIM == 1) then
     sim_xCos = 1.
     sim_yCos = 0.
     sim_zCos = 0.
     
  elseif (NDIM == 2) then
     sim_yCos = sqrt(1. - sim_xCos*sim_xCos)
     sim_zCos = 0.
     
  elseif (NDIM == 3) then
     sim_yCos = cos(sim_yAngle)
     sim_zCos = sqrt( max(0., 1. - sim_xCos*sim_xCos - sim_yCos*sim_yCos) )
  endif

 if (sim_meshMe == MASTER_PE) then

    write (*,*)
    write (*,*) 'flash:  initializing two blast wave problem.'
    write (*,*)
    write (*,1) 'sim_xCos = ', sim_xCos, 'sim_yCos = ', sim_yCos, 'sim_zCos = ', sim_zCos
    write (*,1) 'sim_rhoLeft = ', sim_rhoLeft,  'sim_pLeft = ', sim_pLeft, 'sim_uLeft = ', sim_uLeft
    write (*,1) 'sim_rhoMid = ', sim_rhoMid,'sim_pMid = ', sim_pMid,  'sim_uMid = ', sim_uMid
    write (*,1) 'sim_rhoRight =',sim_rhoRight,'sim_pRight  = ',sim_pRight,'sim_uRight = ',sim_uRight
    write (*,1) 'sim_posnL = ', sim_posnL,'sim_posnR = ', sim_posnR,'sim_gamma = ', sim_gamma
    write (*,2) 'ndim = ', NDIM
    write (*,*)
    
1   format (1X, 1P, 4(A14, E13.7, :, 1X))
2   format (1X, 1P, A7, I13)
    
 endif
         
end subroutine Simulation_init







