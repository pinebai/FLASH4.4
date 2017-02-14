!!****if* source/Simulation/SimulationMain/ReinickeMeyer/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId)
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for a
!!  specified block. Initial profiles are set according to the
!!  analytic solution for a conducting blast wave, derived in:
!!
!! Reinicke, P., Meyer-ter-Vehn, J., The point explosion with heat conduction, 
!! Phys. Fluids A, 1807, 3, 1991
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!
!! PARAMETERS
!!
!!***
!!REORDER(4): U
subroutine Simulation_initBlock(blockId)

  use Grid_interface, ONLY: Grid_getBlkIndexLimits
  use Grid_interface, ONLY: Grid_getBlkPtr 
  use Grid_interface, ONLY: Grid_releaseBlkPtr

  use Simulation_interface, ONLY: Simulation_computeAnalytical

  use Simulation_data, ONLY: sim_tinitial

  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer, intent(IN) ::  blockId
    
  integer :: i
  integer :: j
  integer :: k
  integer :: blkLimits(2,MDIM)
  integer :: blkLimitsGC(2,MDIM)

  real, pointer :: U(:,:,:,:)
  
  ! Compute the analytic solution. The results will be stored in the
  ! UNK variables:
  ! 
  ! ARHO_VAR - density
  ! ATMP_VAR - temperature
  ! AVLX_VAR - X velocity
  ! AVLY_VAR - Y velocity
  ! AVLZ_VAR - Z velocity
  call Simulation_computeAnalytical(blockID, sim_tinitial)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkPtr(blockId,U)
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)    
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH, IAXIS)
                      
           ! Set initial conditions using the analytic solution:
           U(DENS_VAR,i,j,k) = U(ARHO_VAR,i,j,k)
           U(TEMP_VAR,i,j,k) = U(ATMP_VAR,i,j,k)
           U(VELX_VAR,i,j,k) = U(AVLX_VAR,i,j,k)
           U(VELY_VAR,i,j,k) = U(AVLY_VAR,i,j,k)
           U(VELZ_VAR,i,j,k) = U(AVLZ_VAR,i,j,k)
        enddo
     enddo
  enddo
  call Grid_releaseBlkPtr(blockId,U)

  return
  
end subroutine Simulation_initBlock

