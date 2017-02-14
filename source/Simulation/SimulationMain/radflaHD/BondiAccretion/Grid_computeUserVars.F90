!!****if* source/Simulation/SimulationMain/radflaHD/BondiAccretion/Grid_computeUserVars
!!
!! NAME
!!  Grid_computeUserVars
!!
!!
!! SYNOPSIS
!!
!!  call Grid_computeUserVars() 
!!
!!  
!! DESCRIPTION 
!!  
!!  Prepare variables for output.
!!
!!  This version just calls Eos_wrapped because we want flux-limiter-scaled
!!  radiation pressures to show in plot and checkpoint files.
!!
!!  Nothing is done unless the variable sim_plotScaledPressures is TRUE.
!!
!! ARGUMENTS 
!!  none
!!
!!
!! NOTES
!!
!!  The present version calls Eos_wrapped on all blocks, not just leaf blocks.
!!***

subroutine Grid_computeUserVars()
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getCellCoords,  &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped
  use Simulation_data, ONLY : sim_plotScaledPressures

  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer :: blockID
  integer :: localNumBlocks
  integer :: j,k

  integer,dimension(2,MDIM) :: bl,blkLimitsGC
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  real,allocatable :: xc(:)

  call Grid_getLocalNumBlks(localNumBlocks)


  do blockID=1,localNumBlocks

     call Grid_getBlkIndexLimits(blockID,bl,blkLimitsGC)
     allocate(xc(bl(1,IAXIS):bl(2,IAXIS)))
     call Grid_getCellCoords(IAXIS, blockID, CENTER, .FALSE., xc, size(xc))

     call Grid_getBlkPtr(blockID,solnVec)     
     solnVec(FLAD_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) = &
          solnVec(FLLM_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) - &
          solnVec(FLLA_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS))
     solnVec(FLAQ_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) = &
          solnVec(FLLM_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) / &
          solnVec(FLLA_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS))
     solnVec(ERDD_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) = &
          solnVec(ERAD_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) - &
          solnVec(ERDA_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS))
     solnVec(ERDQ_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) = &
          solnVec(ERAD_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) / &
          solnVec(ERDA_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS))
     solnVec(URDD_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) = &
          solnVec(ERAD_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS))&
          *solnVec(DENS_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) - &
          solnVec(URDA_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS))
     solnVec(URDQ_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) = &
          solnVec(ERAD_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS))&
          *solnVec(DENS_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS)) / &
          solnVec(URDA_VAR,bl(1,IAXIS):bl(2,IAXIS),bl(1,JAXIS):bl(2,JAXIS),bl(1,KAXIS):bl(2,KAXIS))
     do k=bl(1,KAXIS),bl(2,KAXIS)
        do j=bl(1,JAXIS),bl(2,JAXIS)
           solnVec(MAR0_VAR,bl(1,IAXIS):bl(2,IAXIS),j,k) = &
                4.0 * PI * &
                solnVec(VELX_VAR,bl(1,IAXIS):bl(2,IAXIS),j,k) &
                *solnVec(DENS_VAR,bl(1,IAXIS):bl(2,IAXIS),j,k) &
                *xc(:)*xc(:)
        end do
     end do
     call Grid_releaseBlkPtr(blockID,solnVec)     
     deallocate(xc)

     if (sim_plotScaledPressures) then
        call Eos_wrapped(MODE_DENS_EI_MAT_GATHER_PRADSCALE, bl, blockID)
     end if
  end do


end subroutine Grid_computeUserVars
