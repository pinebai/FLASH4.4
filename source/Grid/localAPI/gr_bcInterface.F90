!!****ih* source/Grid/localAPI/gr_bcInterface
!!
!! NAME 
!!
!!   gr_bcInterface
!!
!! SYNOPSIS
!!   
!!   use gr_bcInterface
!!
!! DESCRIPTION
!! 
!! This is the header file for the GridBoundaryConditions
!! subunit that defines its interfaces.
!!
!!***

module gr_bcInterface
  implicit none
#include "constants.h"
#include "Flash.h"

  interface
     subroutine gr_bcApplyToAllBlks(axis,isWork)
       integer, intent(in) :: axis
       logical, intent(in) :: isWork
     end subroutine gr_bcApplyToAllBlks
  end interface

  interface
     subroutine gr_bcApplyToOneFace(axis,bcType,gridDataStruct,varCount,&
          regionType,blkLimits,blkLimitsGC,blockHandle,idest)
       integer, intent(in) :: axis,bcType,gridDataStruct,varCount,blockHandle,idest
       integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC,blkLimits
       integer,dimension(MDIM),intent(IN) :: regionType
     end subroutine gr_bcApplyToOneFace
  end interface

  interface 
     subroutine gr_bcGetRegion(gridDataStruct,axis,endPoints,regionSize,mask,&
          region,blockID,idest)
       integer, intent(in) :: gridDataStruct,axis
       integer,dimension(LOW:HIGH,MDIM),intent(IN) :: endPoints
       integer,intent(IN) :: regionSize(REGION_DIM)
       logical,dimension(regionSize(STRUCTSIZE)),intent(OUT) :: mask
       real,dimension(regionSize(BC_DIR),regionSize(SECOND_DIR),&
            regionSize(THIRD_DIR),regionSize(STRUCTSIZE)),&
            intent(OUT) :: region
       integer, intent(in) :: blockID
       integer,intent(IN) :: idest
     end subroutine gr_bcGetRegion
     subroutine gr_bcGetRegionsMixedGds(gridDataStruct,axis,secondDir,thirdDir,endPoints,&
          regionSize,&
          regionDataC,regionDataFN,regionDataFT1,regionDataFT2,&
          blockID,idest)
       implicit none
       integer, intent(in) :: gridDataStruct,axis, secondDir,thirdDir
       integer,dimension(LOW:HIGH,MDIM),intent(IN) :: endPoints
       integer,intent(IN) :: regionSize(REGION_DIM)
       real,pointer,dimension(:,:,:,:) :: regionDataFN, regionDataFT1, regionDataFT2, regionDataC
       integer, intent(in) :: blockID
       integer,intent(IN) :: idest
     end subroutine gr_bcGetRegionsMixedGds
  end interface

  interface 
     subroutine gr_bcPutRegion(gridDataStruct,axis,endPoints,regionSize,mask,&
          region,blockID,idest)
       integer, intent(in) :: gridDataStruct,axis
       integer,dimension(LOW:HIGH,MDIM),intent(IN) :: endPoints
       integer,intent(IN) :: regionSize(REGION_DIM)
       logical,dimension(regionSize(STRUCTSIZE)),intent(IN) :: mask
       real,dimension(regionSize(BC_DIR),regionSize(SECOND_DIR),&
            regionSize(THIRD_DIR),regionSize(STRUCTSIZE)),&
            intent(IN) :: region
       integer, intent(in) :: blockID
       integer,intent(IN) :: idest
     end subroutine gr_bcPutRegion
     subroutine gr_bcPutRegionsMixedGds(gridDataStruct,axis,secondDir,thirdDir,endPoints,&
          regionSize,&
          regionDataC,regionDataFN,regionDataFT1,regionDataFT2,&
          blockID,idest)
       implicit none
       integer, intent(in) :: gridDataStruct,axis, secondDir,thirdDir
       integer,dimension(LOW:HIGH,MDIM),intent(IN) :: endPoints
       integer,intent(IN) :: regionSize(REGION_DIM)
       real,pointer,dimension(:,:,:,:) :: regionDataFN, regionDataFT1, regionDataFT2, regionDataC
       integer, intent(in) :: blockID
       integer,intent(IN) :: idest
     end subroutine gr_bcPutRegionsMixedGds
  end interface

  interface
     subroutine gr_bcMapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct, &
          axis,face,idest)
       integer, intent(OUT) :: bcTypeToApply
       integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct,axis,face
       integer,intent(IN),OPTIONAL:: idest
     end subroutine gr_bcMapBcType
  end interface

  interface
     subroutine gr_bcInit()
     end subroutine gr_bcInit
  end interface

  interface
     subroutine gr_bcHseInit()
     end subroutine gr_bcHseInit
  end interface

  interface
     subroutine gr_hseStep(dens, temp, ye, sumy, n, inputg, delta, direction, order, mode, massFrac)
       implicit none
       real, dimension(:), intent(IN)    :: ye, sumy
       real, dimension(:), intent(INOUT) :: dens, temp
       integer, intent(IN)               :: n               ! index to update
       real, intent(IN)                  :: inputg, delta
       integer, intent(IN)               :: direction, order, mode
       real, optional, intent(IN)        :: massFrac(:)
     end subroutine gr_hseStep
  end interface

  interface
     subroutine gr_bcFinalize()
     end subroutine gr_bcFinalize
  end interface
end module gr_bcInterface






