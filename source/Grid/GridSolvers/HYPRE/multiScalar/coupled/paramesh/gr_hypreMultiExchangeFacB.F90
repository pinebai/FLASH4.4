!!****if* source/Grid/GridSolvers/HYPRE/multiScalar/coupled/paramesh/gr_hypreMultiExchangeFacB
!!
!!  NAME 
!!
!!  gr_hypreMultiExchangeFacB
!!
!!  SYNOPSIS
!!
!!  call gr_hypreMultiExchangeFacB (integer,intent(IN) :: iFactorB,
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


subroutine gr_hypreMultiExchangeFacB (unkVarsDesc, firstHypreVar,diffCoeffDesc, blockCount, blockList)
  
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_putFluxData, Grid_getBlkData, &
    Grid_conserveFluxes, &
    Grid_ascGetBlkPtr, Grid_ascReleaseBlkPtr
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_interface,           ONLY: gr_hypreGetFaceB
  use gr_hypreLocalInterface, ONLY: gr_hypreGetFaceBFcB
  use gr_hypreData,     ONLY : gr_hypreNParts

  use physicaldata, ONLY: consv_flux_densities
  
  implicit none
  
#include "Flash.h"  
#include "constants.h"
#include "FortranLangFeatures.fh"
  
  integer,intent(IN) :: unkVarsDesc(VARDESC_SIZE)
  integer,intent(IN) :: firstHypreVar
  integer,intent(IN) :: diffCoeffDesc(VARDESC_SIZE)
  integer,intent(IN) :: blockCount
  integer, dimension(blockCount),intent(IN):: blockList
  
  
  !! LOCAL VARIABLES
  integer :: iFactorB
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  real, POINTER, DIMENSION(:,:,:,:) :: facBptr
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits
  integer :: i, numVars, iv
  integer :: blockID
  integer :: lb
  integer, parameter :: fpv = 2**(NDIM-1) ! fluxes per variable
  integer :: nFluxVars
  integer, dimension(:),allocatable :: pSlots
  logical :: doFcB

  real, allocatable, dimension(:,:,:,:) :: anyFlux
  real, POINTER, DIMENSION(:,:,:) :: facBptr1

  real, allocatable :: areaLeft  (:,:,:)  
  
  integer :: datasizeGC(MDIM)

  real :: scaleFlux
  
  if (gr_hypreNParts == 1) then !! paramesh in UG mode.     
     return
  end if
     
  call Timers_start("gr_hypreMultiExchangeFacB")       

  numVars = unkVarsDesc(VARDESC_NUM)
  iFactorB = diffCoeffDesc(VARDESC_VAR)
  nFluxVars = fpv * numVars
  doFcB = (diffCoeffDesc(VARDESC_GDS) == FACES .AND. diffCoeffDesc(VARDESC_DURATION) == VD_DUR_GASC)
  
  if (blockCount > 0) then 

     allocate(pSlots(numVars))
     do i=1,nFluxVars
        pSlots(i) = i
     end do

     scaleFlux = 1.0
     
     if (consv_flux_densities) scaleFlux = real(2**(NDIM-1))
     
     
     do lb = 1, blockCount
        
        blockID = blockList(lb)    
     
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     
        
        datasizeGC  (1:MDIM)= blkLimitsGC  (HIGH,1:MDIM)-blkLimitsGC  (LOW,1:MDIM)+1          

        allocate(areaLeft(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),   &
             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),   &
             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
        allocate(anyFlux(NFLUXES,blkLimitsGC(HIGH,IAXIS), &
             blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)))
     
        anyFlux = 0.0
        if (doFcB) then
           call Grid_ascGetBlkPtr(blockID,facBptr,FACEX)
           do iv = 0, numVars-1
              call AssoMed(facBptr1,facBptr,iFactorB+iv) ! facBptr1(:,:,:) => facBptr(iFactorB+iv,:,:,:)
              call gr_hypreGetFaceBFcB (IAXIS, blkLimits, blkLimitsGC, facBptr1, anyFlux, iv)
           end do
        else
           call Grid_getBlkPtr(blockID, solnVec)          
           call gr_hypreGetFaceB (IAXIS, iFactorB, blkLimits, blkLimitsGC, solnVec, anyFlux, numVars=numVars)
        end if
        anyFlux = anyFlux*scaleFlux
        call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR,& 
             blkLimitsGC(LOW,:), areaLeft(:,:,:), datasizeGC) 
        call Grid_putFluxData(blockID,IAXIS,anyFlux, blkLimitsGC(HIGH,:), &
             pressureSlots=pSlots(1:nFluxVars),areaLeft=areaLeft)
        if (doFcB) call Grid_ascReleaseBlkPtr(blockID,facBptr,FACEX)

#if NDIM >= 2    
        anyFlux = 0.0
        if (doFcB) then
           call Grid_ascGetBlkPtr(blockID,facBptr,FACEY)
           do iv = 0, numVars-1
              call AssoMed(facBptr1,facBptr,iFactorB+iv) ! facBptr1(:,:,:) => facBptr(iFactorB+iv,:,:,:)
              call gr_hypreGetFaceBFcB (JAXIS, blkLimits, blkLimitsGC, facBptr1, anyFlux, iv)
           end do
        else
           call gr_hypreGetFaceB (JAXIS, iFactorB, blkLimits, blkLimitsGC, solnVec, anyFlux, numVars=numVars)
        end if
        anyFlux = anyFlux*scaleFlux
        call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR,& 
             blkLimitsGC(LOW,:), areaLeft(:,:,:), datasizeGC)
        call Grid_putFluxData(blockID,JAXIS,anyFlux,blkLimitsGC(HIGH,:), &
             pressureSlots=pSlots(1:nFluxVars),areaLeft=areaLeft)
        if (doFcB) call Grid_ascReleaseBlkPtr(blockID,facBptr,FACEY)
     
#if NDIM == 3
        anyFlux = 0.0
        if (doFcB) then
           call Grid_ascGetBlkPtr(blockID,facBptr,FACEZ)
           do iv = 0, numVars-1
              call AssoMed(facBptr1,facBptr,iFactorB+iv) ! facBptr1(:,:,:) => facBptr(iFactorB+iv,:,:,:)
              call gr_hypreGetFaceBFcB (KAXIS, blkLimits, blkLimitsGC, facBptr1, anyFlux, iv)
           end do
        else
           call gr_hypreGetFaceB (KAXIS, iFactorB, blkLimits, blkLimitsGC, solnVec, anyFlux, numVars=numVars)
        end if
        anyFlux = anyFlux*scaleFlux
        call Grid_getBlkData(blockID, CELL_FACEAREA, KLO_FACE, EXTERIOR,& 
             blkLimitsGC(LOW,:), areaLeft(:,:,:), datasizeGC)
        call Grid_putFluxData(blockID,KAXIS,anyFlux,blkLimitsGC(HIGH,:), &
             pressureSlots=pSlots,areaLeft=areaLeft)
        if (doFcB) call Grid_ascReleaseBlkPtr(blockID,facBptr,FACEZ)
#endif     
#endif  
        deallocate (anyFlux)
        
        if (.NOT. doFcB) call Grid_releaseBlkPtr(blockID, solnVec)            

        deallocate(areaLeft)

     end do
     
     deallocate(pSlots)
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
  
  
  call Timers_stop("gr_hypreMultiExchangeFacB") 
  
  return
  
contains
  subroutine AssoMed(pp, mm, varNo)
    real,POINTER_INTENT_OUT :: pp(:,:,:)
    real,POINTER_INTENT_IN  :: mm(:,:,:,:)
    integer,intent(in) :: varNo
    call AssoFin(pp,mm(varNo,:,:,:),lbound(mm,1),lbound(mm,2),lbound(mm,3),lbound(mm,4))
  end subroutine AssoMed

  subroutine AssoFin(pp, dd, lb1,lb2,lb3,lb4)
    real,POINTER_INTENT_OUT :: pp(:,:,:)
    integer, intent(in) :: lb1,lb2,lb3,lb4
    real,   intent(in),target :: dd(lb2:,lb3:,lb4:)
    pp => dd
  end subroutine AssoFin
end subroutine gr_hypreMultiExchangeFacB
