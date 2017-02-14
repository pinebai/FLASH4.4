!!****if* source/Simulation/SimulationMain/radflaHD/BondiAccretion/Simulation_computeAnalytical
!!
!! NAME
!!
!!  Simulation_computeAnalytical
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
!!  Initializes analuytical solution for
!!  a specified block.
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  tcurr   -        current time (unused here)
!!
!!
!!***

subroutine Simulation_computeAnalytical(blockId, tcurr)

  use Simulation_data, ONLY : sim_rho_vac
  use Driver_interface, ONLY: Driver_abortFlash
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_advanceDiffusion, Grid_getBlkIndexLimits, Grid_fillGuardCells, &
      Grid_getDeltas, Grid_getCellCoords

  use sim_interface, ONLY: sim_computeAnaBondi

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

           call sim_computeAnaBondi(xcent(i),solnVec(VELA_VAR,i,j,k),solnVec(DENA_VAR,i,j,k))
           
!           write (*,'(1p8e14.6)')  xcent(i), solnVec(DENA_VAR,i,j,k) / sim_rho_vac, &
!                solnVec(VELA_VAR,i,j,k), solnVec(VELA_VAR,i,j,k)*solnVec(DENA_VAR,i,j,k)
           
           
           solnVec(FLAD_VAR,i,j,k) = solnVec(FLLM_VAR,i,j,k) - solnVec(FLLA_VAR,i,j,k)
           solnVec(FLAQ_VAR,i,j,k) = solnVec(FLLM_VAR,i,j,k) / solnVec(FLLA_VAR,i,j,k)
           solnVec(ERDD_VAR,i,j,k) = solnVec(ERAD_VAR,i,j,k) - solnVec(ERDA_VAR,i,j,k)
           solnVec(ERDQ_VAR,i,j,k) = solnVec(ERAD_VAR,i,j,k) / solnVec(ERDA_VAR,i,j,k)
           solnVec(URDD_VAR,i,j,k) = solnVec(ERAD_VAR,i,j,k)*solnVec(DENS_VAR,i,j,k) - solnVec(URDA_VAR,i,j,k)
           solnVec(URDQ_VAR,i,j,k) = solnVec(ERAD_VAR,i,j,k)*solnVec(DENS_VAR,i,j,k) / solnVec(URDA_VAR,i,j,k)

        enddo
     enddo
  enddo
  
  deallocate (xcent)
  
  call Grid_releaseBlkPtr(blockID,solnVec) 
    
  return

end subroutine Simulation_computeAnalytical
