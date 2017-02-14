!!****if* source/Simulation/SimulationMain/Plasma/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockId)
!!
!!
!!
!! DESCRIPTION
!!  This routine applies initial conditions of an
!!  ion beam into plasma test problem for hybrid-pic
!!
!! 
!! ARGUMENTS
!!
!!  blockId -         the number of the block to update
!!
!!***
subroutine Simulation_initBlock(blockId)
  use Simulation_data, ONLY : sim_bx,sim_by, sim_bz
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) ::  blockId
  
  real, dimension(:,:,:,:), pointer :: U

  call Grid_getBlkPtr(blockId, U)
  U(GRBX_VAR,:,:,:) = sim_bx
  U(GRBY_VAR,:,:,:) = sim_by
  U(GRBZ_VAR,:,:,:) = sim_bz

  ! External field
  U(GBX1_VAR,:,:,:) = 0.0
  U(GBY1_VAR,:,:,:) = 0.0
  U(GBZ1_VAR,:,:,:) = 0.0

  call Grid_releaseBlkPtr(blockId, U)

end subroutine Simulation_initBlock
