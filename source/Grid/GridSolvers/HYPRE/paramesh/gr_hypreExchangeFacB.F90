!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreExchangeFacB
!!
!!  NAME 
!!
!!  gr_hypreExchangeFacB
!!
!!  SYNOPSIS
!!
!!  call gr_hypreExchangeFacB (integer,intent(IN) :: iFactorB,
!!                             integer,intent(IN) :: blockCount,
!!                             integer, dimension(blockCount),intent(IN):: blockList)
!!
!!  DESCRIPTION 
!!   This routine helps exchange iFactorB at fine-coarse boundaries. If called
!!   in UG mode this rountine returns without any action. In AMR mode it uses (tricks) 
!!   Grid_conserveFluxes to succesfully communicate the values across processors. Doing
!!   this saves us from performing more complicated point-point communication.
!!
!!
!! ARGUMENTS
!!   iFactorB   : Factor from the diffusion equation.
!!   blockCount : The number of blocks in the list.   
!!   blockList  : The list of blocks on which the solution must be updated.  
!!
!!
!! SIDE EFFECTS
!!
!!   Upon return, data will have been stored (and corrected) in the global storage
!!   area that he Grid unit uses for storage of flux data.  The data can be retrieved
!!   by calling Grid_getFluxData.
!!
!! NOTES
!!
!!***


subroutine gr_hypreExchangeFacB (iFactorB, blockCount, blockList)
  
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_putFluxData, Grid_getBlkData, &
    Grid_conserveFluxes
  use gr_interface, ONLY   : gr_hypreGetFaceB
  use Timers_interface, ONLY : Timers_start, Timers_stop  
  use gr_hypreData,     ONLY : gr_hypreNParts

  use physicaldata, ONLY: consv_flux_densities
  
  implicit none
  
#include "Flash.h"  
#include "constants.h"
  
  integer,intent(IN) :: iFactorB
  integer,intent(IN) :: blockCount
  integer, dimension(blockCount),intent(IN):: blockList
  
  
  !! LOCAL VARIABLES
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits  
  integer :: blockID
  integer :: lb
  integer, parameter :: nFluxVars = 2**(NDIM-1)
  integer, dimension(2**(MDIM-1)),parameter :: pSlots = (/1,2,3,4/)

  real, allocatable, dimension(:,:,:,:) :: anyFlux

  real, allocatable :: areaLeft  (:,:,:)  
  
  integer :: datasizeGC(MDIM)

  real :: scaleFlux
  
  if (gr_hypreNParts == 1) then !! paramesh in UG mode.     
     return
  end if
     
  call Timers_start("gr_hypreExchangeFacB")       

  
  
  if (blockCount > 0) then 

     scaleFlux = 1.0
     
     if (consv_flux_densities) scaleFlux = 2**(NDIM-1)
     
     
     do lb = 1, blockCount
        
        blockID = blockList(lb)    
     
        call Grid_getBlkPtr(blockID, solnVec)          
        
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     
        
        datasizeGC  (1:MDIM)= blkLimitsGC  (HIGH,1:MDIM)-blkLimitsGC  (LOW,1:MDIM)+1          

        allocate(areaLeft(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),   &
             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),   &
             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
        allocate(anyFlux(NFLUXES,blkLimitsGC(HIGH,IAXIS), &
             blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
     
        anyFlux = 0.0
        call gr_hypreGetFaceB (IAXIS, iFactorB, blkLimits, blkLimitsGC, solnVec, anyFlux, numVars=1)
        anyFlux = anyFlux*scaleFlux
        call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR,& 
             blkLimitsGC(LOW,:), areaLeft(:,:,:), datasizeGC) 
        call Grid_putFluxData(blockID,IAXIS,anyFlux, blkLimitsGC(HIGH,:), &
             pressureSlots=pSlots(1:nFluxVars),areaLeft=areaLeft)

#if NDIM >= 2    
        anyFlux = 0.0
        call gr_hypreGetFaceB (JAXIS, iFactorB, blkLimits, blkLimitsGC, solnVec, anyFlux, numVars=1)
        anyFlux = anyFlux*scaleFlux
        call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR,& 
             blkLimitsGC(LOW,:), areaLeft(:,:,:), datasizeGC)
        call Grid_putFluxData(blockID,JAXIS,anyFlux,blkLimitsGC(HIGH,:), &
             pressureSlots=pSlots(1:nFluxVars),areaLeft=areaLeft)
     
#if NDIM == 3
        anyFlux = 0.0
        call gr_hypreGetFaceB (KAXIS, iFactorB, blkLimits, blkLimitsGC, solnVec, anyFlux, numVars=1)
        anyFlux = anyFlux*scaleFlux
        call Grid_getBlkData(blockID, CELL_FACEAREA, KLO_FACE, EXTERIOR,& 
             blkLimitsGC(LOW,:), areaLeft(:,:,:), datasizeGC)
        call Grid_putFluxData(blockID,KAXIS,anyFlux,blkLimitsGC(HIGH,:), &
             pressureSlots=pSlots,areaLeft=areaLeft)
#endif     
#endif  
        deallocate (anyFlux)
        
        call Grid_releaseBlkPtr(blockID, solnVec)            

        deallocate(areaLeft)

     end do
     
  end if
  
  
!!$  call Grid_conserveFluxes(IAXIS, 0)  
!!$#if NDIM >= 2  
!!$  call Grid_conserveFluxes(JAXIS, 0)  
!!$#if NDIM == 3
!!$  call Grid_conserveFluxes(KAXIS, 0)
!!$#endif
!!$#endif

  !! ALLDIR: 0
  call Grid_conserveFluxes(ALLDIR, 0)
  
  
  call Timers_stop("gr_hypreExchangeFacB") 
  
  return
  
end subroutine gr_hypreExchangeFacB
