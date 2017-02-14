!!****ih* source/physics/Hydro/HydroMain/split/PPM/chomboCompatible/hy_ppm_chombo_interface
!!
!! NAME
!!
!!  hy_ppm_chombo_interface
!!
!! SYNOPSIS
!!
!!  use hy_ppm_chombo_interface
!!
!! DESCRIPTION
!!
!!  Interface module for internal use within the split PPM Hydro implementation.
!!
!! NOTES
!!  CD: I added the interface module to work around a bug in the Intel compiler:
!!  http://software.intel.com/en-us/articles/error-8000-interface-block-function-returning-expression-sized-array/
!!  Previously I had an interface for hy_ppm_chomboLikeUpdateSoln inside
!!  hy_ppm_updateSoln subroutine, but this was not accepted by Intel compiler.
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData, tempFlx

Module hy_ppm_chombo_interface 
  implicit none

#include "constants.h"
#include "Flash.h"
    
  interface
     subroutine hy_ppm_chomboLikeUpdateSoln(rangeSwitch, &
          xyzswp, dt,                          &
          blkLimits,blkLimitsGC,numCells,      &
          tempArea, tempGrav1d_o, tempGrav1d,  &
          tempDtDx, tempFict,                  &
          tempFlx,  solnData )
       
       implicit none
       integer, intent(IN) :: rangeSwitch
       integer, intent(IN) :: xyzswp
       real,    intent(IN) :: dt
       integer, intent(IN) :: numCells
       integer, intent(IN),dimension(2,MDIM)::blkLimitsGC,blkLimits
#ifdef FIXEDBLOCKSIZE
       real, intent(IN), DIMENSION(GRID_ILO_GC:GRID_IHI_GC, &
            GRID_JLO_GC:GRID_JHI_GC, &
            GRID_KLO_GC:GRID_KHI_GC  ) :: &
            tempArea, tempGrav1d_o,  &
            tempGrav1d, &
            tempDtDx, tempFict
       real, intent(IN), DIMENSION(NFLUXES,GRID_ILO_GC:GRID_IHI_GC, &
            GRID_JLO_GC:GRID_JHI_GC, &
            GRID_KLO_GC:GRID_KHI_GC) :: tempFlx
#else
       real, intent(IN), DIMENSION(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)  ) :: &
            tempArea, tempGrav1d_o, &
            tempGrav1d, &
            tempDtDx, tempFict
       real, intent(IN), DIMENSION(NFLUXES,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: tempFlx
       
#endif
       real, pointer :: solnData(:,:,:,:) 
       
     end subroutine hy_ppm_chomboLikeUpdateSoln
  end interface


  interface
     subroutine hy_ppm_flux_conserve(operation, xyzswp, blkCount, blkList)
       implicit none
       integer, intent(in) :: operation, xyzswp, blkCount
       integer, dimension(blkCount), intent(in) :: blkList
     end subroutine hy_ppm_flux_conserve
  end interface

end Module hy_ppm_chombo_interface
