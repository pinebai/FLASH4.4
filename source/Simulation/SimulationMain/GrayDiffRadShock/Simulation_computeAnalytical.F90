!!****if* source/Simulation/SimulationMain/GrayDiffRadShock/Simulation_computeAnalytical
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_computeAnalytical(integer(IN) :: blockID, 
!!                                    real(IN)    :: tcurr)
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
!!  tcurr   -        current time
!!
!!
!!***

subroutine Simulation_computeAnalytical(blockId, tcurr)

  use Simulation_data, ONLY : sim_rho, sim_temp
  use Conductivity_interface, ONLY : Conductivity
  use Eos_interface, ONLY : Eos
  use Driver_interface, ONLY: Driver_abortFlash
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_advanceDiffusion, Grid_getBlkIndexLimits, Grid_fillGuardCells, &
      Grid_getDeltas, Grid_getCellCoords

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"   

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  real,    intent(in) :: tcurr

  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec

  real, allocatable :: xcent(:)
  
  call Grid_getBlkPtr(blockID,solnVec)     
  
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., &
       xcent, blkLimitsGC(HIGH, IAXIS))
  
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)         
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)     
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)                    
           
           solnVec(TION_VAR,i,j,k) = solnVec(TION_VAR,i,j,k) / sim_temp
           solnVec(TELE_VAR,i,j,k) = solnVec(TELE_VAR,i,j,k) / sim_temp
           solnVec(TRAD_VAR,i,j,k) = solnVec(TRAD_VAR,i,j,k) / sim_temp
           
           solnVec(DENS_VAR,i,j,k) = solnVec(DENS_VAR,i,j,k) / sim_rho                  
           
           write (*,'(1p8e14.6)')  xcent(i), solnVec(DENS_VAR,i,j,k) / sim_rho, 1.1604505E6*solnVec(TELE_VAR,i,j,k) / sim_temp, 1.1604505E6*solnVec(TRAD_VAR,i,j,k) / sim_temp
           
        enddo
     enddo
  enddo
  
  deallocate (xcent)
  
  call Grid_releaseBlkPtr(blockID,solnVec) 
    
  return

end subroutine Simulation_computeAnalytical
