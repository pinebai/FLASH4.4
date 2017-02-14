!!****if* source/Grid/GridBoundaryConditions/gr_bcApplyToOneFaceAllGds
!!
!! NAME
!!  gr_bcApplyToOneFaceAllGds
!!
!! SYNOPSIS
!!
!!  call gr_bcApplyToOneFaceAllGds(integer(IN) :: axis,
!!                      integer(IN) :: bcType,
!!                      integer(IN) :: gridDataStruct,
!!                      integer(IN) :: varCount,
!!                      integer(IN) :: regionType(MDIM)
!!                      integer(IN) :: blkLimits(LOW:HIGH,MDIM)
!!                      integer(IN) :: blkLimitsGC(LOW:HIGH,MDIM)
!!                      integer(IN) :: blockHandle,
!!                      integer(IN) :: idest)
!!  
!! DESCRIPTION 
!!
!!  
!! 
!! ARGUMENTS
!!  
!!    axis           - the direction for applying BC, one of IAXIS, JAXIS, or KAXIS
!!    bcType         - the type of boundary condition
!!    gridDataStruct - In PM3 and PM4 it can have values (CENTER,FACEX,FACEY,FACEZ,WORK),
!!                     in UG (CENTER,FACEX,FACEY,FACEZ), and in PM2 (CENTER,WORK).
!!    varCount       - the number of variable in the data structure specified in gridDataStruct. IGNORED??
!!    regionType     - The part of the block that is involved in the boundary condition. This integer
!!                     array can have values (LEFT_EDGE, CENTER, RIGHT_EDGE, WHOLE_VEC and NO_VEC)
!!                     for each of the three dimensions of the physical data.
!!                     LEFT_EDGE implies guard cells along lower face. If this value is specified
!!                     for the axis along which the BC are being applies, that means that we 
!!                     are applying BC to the lowerface, if it is one of the other axes, then
!!                     the BC are being applied to one of the corners. Same is true of RIGHT_EDGE
!!                     except that implies the upperface. CENTER, WHOLE_VEC and NO_VEC values
!!                     are valid only for dimensions other than the one specified in "axis". 
!!                     CENTER implies only the interior cells, WHOLE_VEC implies all cells in 
!!                     the block and NO_VEC implies that the correspoding dimension is not
!!                     a part of the region. Normally this value is most likely to be used
!!                     along KAXIS in a 2D problems, and JAXIS and KAXIS in a 1D problem
!!    blkLimits      - the endpoints of the block cell (or face) indices without the guardcells
!!    blkLimitsGC    - the endpoints of the block cell (or face) indices including the guardcells
!!    blockHandle    - local block number; with Paramesh 4, this may be a handle for a remote block
!!                     for which Paramesh has cached information locally.
!!    idest         - this is used in Paramesh 4, where boundary condition handling is applied
!!                    not to a block's solution data in their permanent location (named UNK, etc.),
!!                    but to a buffer that contain temporary copies of a few (normally, two) blocks'
!!                    solution data (named UNK1, etc.). The dummy argument idest selects a buffer
!!                    slot in that case. It is ignored in other Grid implementations.
!!
!! NOTES
!!  A specific direction is required in axis - no ALLDIR at this point.
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_bcApplyToOneFaceAllGds(axis,bcType,gridDataStruct,varCount,&
     regionType,blkLimits,blkLimitsGC,blockHandle,idest)

  use Grid_interface, ONLY : Grid_bcApplyToRegion, &
       Grid_bcApplyToRegionSpecialized, &
       Grid_bcApplyToRegionMixedGds
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_bcInterface, ONLY : gr_bcGetRegion, gr_bcPutRegion, &
                             gr_bcGetRegionsMixedGds, gr_bcPutRegionsMixedGds
  implicit none
  
  integer, intent(in) :: axis,bcType,gridDataStruct,varCount,blockHandle,idest
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC,blkLimits
  integer,dimension(MDIM),intent(IN) :: regionType

  real,pointer,dimension(:,:,:,:) :: regionDataFN, regionDataFT1, regionDataFT2, regionDataC
  integer :: face,guard
  integer,dimension(LOW:HIGH,MDIM) :: endpointsCtr
  integer,dimension(REGION_DIM) :: regionSizeCtr
  integer,dimension(MDIM-1) :: nextDir
  integer :: i,j,k,n,var,sb,se,sbf,sef
  integer :: numAllVars

  logical :: isFaceVarNormalDir,applied
  logical :: rightHanded

  if(regionType(axis)==LEFT_EDGE) then
     face=LOW
  else
     face=HIGH
  end if


  guard=blkLimits(LOW,axis)-blkLimitsGC(LOW,axis)

  isFaceVarNormalDir = (gridDataStruct==FACEX).and.(axis==IAXIS)
  isFaceVarNormalDir = isFaceVarNormalDir.or.((gridDataStruct==FACEY).and.(axis==JAXIS))
  isFaceVarNormalDir = isFaceVarNormalDir.or.((gridDataStruct==FACEZ).and.(axis==KAXIS))
!!$  isFaceVarNormalDir = isFaceVarNormalDir.or.(gridDataStruct==CENTER_FACES)
  isFaceVarNormalDir = .TRUE.

  if(regionType(axis)==LEFT_EDGE) then
     sb = blkLimits(LOW,axis)-guard
     se = blkLimits(LOW,axis)+guard-1
  else
     sb = blkLimits(HIGH,axis)-guard+1
     se = blkLimits(HIGH,axis)+guard
  end if
  sbf = sb; sef = se+1

  rightHanded = .TRUE.
  if (axis==IAXIS) then
     nextDir(1) = JAXIS
     nextDir(2) = KAXIS
  else if (axis==JAXIS) then
     nextDir(1) = IAXIS         !!
     nextDir(2) = KAXIS         !!
     rightHanded = .FALSE.
  else
     nextDir(1) = IAXIS
     nextDir(2) = JAXIS
  end if

  do i = 1,2
     j=nextDir(i)
     select case(regionType(j))
     case(LEFT_EDGE)
        endpointsCtr(LOW,j)=blkLimitsGC(LOW,j)
        endpointsCtr(HIGH,j)=blkLimits(LOW,j)-1
     case(CENTER)
        endpointsCtr(:,j)=blkLimits(:,j)
     case(RIGHT_EDGE)
        endpointsCtr(LOW,j)=blkLimits(HIGH,j)+1
        endpointsCtr(HIGH,j)=blkLimitsGC(HIGH,j)
     case(WHOLE_VECTOR)
        endpointsCtr(:,j)=blkLimitsGC(:,j)
     case(NO_VEC)
        endpointsCtr(:,j)=1
     end select
     regionSizeCtr(i+1)=endpointsCtr(HIGH,j)-endpointsCtr(LOW,j)+1
  end do

  endpointsCtr(LOW,axis)=sb; endpointsCtr(HIGH,axis)=se
  regionSizeCtr(BC_DIR)=se-sb+1
  regionSizeCtr(STRUCTSIZE)=varCount

  call gr_bcGetRegionsMixedGds(gridDataStruct,axis,nextDir(1),nextDir(2),endpointsCtr,regionSizeCtr,&
       regionDataC,regionDataFN,regionDataFT1,regionDataFT2,blockHandle,idest)
!!$  print*,'Done with gr_bcGetRegionsMixedGds'
#ifdef __INTEL_COMPILER
#define SUBASSERT(asser) call abo(#asser)
#else
#define SUBASSERT(asser) call abo('asser')
#endif
#define ASSERT(assertion) if (.NOT.(assertion)) SUBASSERT(assertion)

  ASSERT(size(regionDataC,1) == regionSizeCtr(BC_DIR))
  ASSERT( size(regionDataC,2) == regionSizeCtr(SECOND_DIR) )
  ASSERT( size(regionDataC,3) == regionSizeCtr(THIRD_DIR) )
  numAllVars = 0
!!$  if (associated(regionDataC  )) print*,'numAllVars = size(regionDataC,4)=',size(regionDataC,4)
!!$  if (associated(regionDataFN )) print*,'numAllVars = numAllVars + size(regionDataFN,4)=',numAllVars,' +',size(regionDataFN,4)
!!$  if (associated(regionDataFT1)) print*,'numAllVars = numAllVars + size(regionDataFT1,4)=',numAllVars,' +', size(regionDataFT1,4)
!!$  if (associated(regionDataFT2)) print*,'numAllVars = numAllVars + size(regionDataFT2,4)=',numAllVars,' +', size(regionDataFT2,4)
  if (associated(regionDataC  )) numAllVars = size(regionDataC,4)
  if (associated(regionDataFN )) numAllVars = numAllVars + size(regionDataFN,4)
  if (associated(regionDataFT1)) numAllVars = numAllVars + size(regionDataFT1,4)
  if (associated(regionDataFT2)) numAllVars = numAllVars + size(regionDataFT2,4)
  ASSERT( numAllVars == regionSizeCtr(STRUCTSIZE) )
!!$  call Grid_bcApplyToRegionSpecialized(bcType,gridDataStruct,&
!!$       guard,axis,face,regionData,regionSize,mask,applied,&
!!$       blockHandle,nextDir(1),nextDir(2),endpointsCtr,blkLimitsGC, idest)
!!$  if(.not.applied) then
!!$     call Grid_bcApplyToRegion(bcType,gridDataStruct,&
!!$          guard,axis,face,regionData,regionSize,mask,applied,&
!!$       blockHandle,nextDir(1),nextDir(2),endpointsCtr,blkLimitsGC, idest)
!!$  end if
  applied = .TRUE.
  call Grid_bcApplyToRegionMixedGds(bcType,gridDataStruct,&
          guard,axis,face,&
          regionDataC,regionDataFN,regionDataFT1,regionDataFT2,&
          regionSizeCtr,&
          applied,&
       blockHandle,nextDir(1),nextDir(2),endpointsCtr,blkLimitsGC, rightHanded,idest)
  if(.not.applied) then
     print*,'gr_bcApplyToOneFace: Unhandled boundary type',bcType, 'axis,regionType=',axis,regionType(axis)
     if (regionType(axis)==LEFT_EDGE) then
        call Driver_abortFlash("unsupported boundary condition on Lower Face")
     else if (regionType(axis)==RIGHT_EDGE) then
         call Driver_abortFlash("unsupported boundary condition on Upper Face")
     else
        call Driver_abortFlash("unexpected regionType!")
     end if
  end if
  call gr_bcPutRegionsMixedGds(gridDataStruct,axis,nextDir(1),nextDir(2),endpointsCtr,regionSizeCtr,&
       regionDataC,regionDataFN,regionDataFT1,regionDataFT2,blockHandle,idest)

  return

contains
  subroutine abo(msg)
    character(len=*) msg
    call Driver_abortFlash("Failed assertion in gr_bcApplyToOneFaceAllGds:"//msg)
  end subroutine abo
end subroutine gr_bcApplyToOneFaceAllGds
