!!****if* source/Grid/GridBoundaryConditions/Chombo/gr_bcPutRegion
!!
!! NAME
!!  gr_bcPutRegion
!!
!! SYNOPSIS
!!
!!  call gr_bcPutRegion(integer(IN) :: gridDataStruct,
!!                 integer(IN)   :: axis,
!!                 integer(IN)   :: endPoints(LOW:HIGH,MDIM),
!!                 integer(IN)   :: regionSize(REGION_DIM),
!!                 integer(IN)   :: mask(regionSize(STRUCTSIZE)),
!!                 real(IN)      :: region(regionSize(1),regionSize(2),regionSize(3),regionSize(4)),
!!                 integer(IN)   :: blockID,
!!                 integer(IN)   :: idest)
!!  
!! DESCRIPTION 
!!  This routine puts the data with boundary conditions applied back into
!!  the Grid-maintained storage for the specified Grid data structure for 
!!  all supported data structures, described for argument gridDataStruct below.
!!  The data are returned in a four-dimensional array, region, where the fourth
!!  dimension represents the individual variables of the data structure.
!!  The first three dimensions store a set of rows with relevant sections of
!!  the block on which the boundary conditions have been applied. Each row
!!  contains complete set of data points to correctly apply the boundary
!!  conditions along the specified axis. The endPoints argument specifies
!!  the bounding box of the regions being selected. For more details, see
!!  the example in gr_bcGetRegion:
!!  
!!
!! ARGUMENTS 
!!
!!
!!  gridDataStruct : optional integer value specifying data structure. 
!!                   The options are defined in constants.h and they are :
!!                   CENTER cell centered variables (default)
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!                   WORK   work array specific to paramesh
!!  axis           : The axis on which boundary condition is being applied
!!  endPoints      : the boundaries of the region to be extracted from the 
!!                   Grid block
!!  regionSize     : regionSize(BC_DIR) contains the size of each row,
!!                   regionSize(SECOND_DIR) contains the number of rows along the
!!                   second direction, and regionSize(THIRD_DIR) has the number of rows
!!                   along the third direction. regionSize(STRUCTSIZE) contains the
!!                   number of variables in the data structure.
!!  mask           : Mask to apply if selected variables are getting boundary filled.
!!                   Currently this only has meaning for PARAMESH4.
!!                   Currently ignored here.
!!  region         : the extracted region
!!  blockID        : the local block ID.
!!  idest          : has meaning only for PARAMESH4, where it distinguishes 
!!                   between leaf and parent nodes; see NOTES below.
!!
!!
!! NOTES
!!  Beginning with PARAMESH3: The updated solution data in the
!!  region array are not actually copied directly to "permanent" storage
!!  (UNK,WORK,etc.), but to the one-block arrays (UNK1,WORK1,etc.)
!!  that are used by PARAMESH while it is processing a block's data.
!!  Calls to gr_bcGetRegion are therefore only meaningful in certain
!!  contexts.  The idest argument is then taken as an index into these
!!  one-block arrays. It distinguishes between the slots available and
!!  must match the slot that is actually being used by PARAMESH.
!!
!! SEE ALSO
!!   gr_bcGetRegion
!!***
#include "Flash.h"

!!REORDER(4): solnData

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_bcPutRegion(gridDataStruct,axis,endPoints,regionSize,mask,&
     region,blockID,idest)
  
#include "constants.h"
  
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  implicit none
  
  integer, intent(in) :: gridDataStruct,axis
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: endPoints
  integer,intent(IN) :: regionSize(REGION_DIM)
  logical,dimension(regionSize(STRUCTSIZE)),intent(IN) :: mask
  real,dimension(regionSize(BC_DIR),regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),regionSize(STRUCTSIZE)),intent(IN) :: region
  integer, intent(in) :: blockID
  integer,intent(IN) :: idest

  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: var,i,j,k,n,m,strt,fin
  logical :: validGridDataStruct
  integer :: varCount, bcVecEnd


  strt = endPoints(LOW,axis)
  fin  = endPoints(HIGH,axis)
  varCount = regionSize(STRUCTSIZE)
  bcVecEnd = regionSize(BC_DIR)

  !An invalid gridDataStruct is caught by Grid_getBlkPtr so
  !there is no extra checking in gr_bcPutRegion
  call Grid_getBlkPtr(blockID, solnData, gridDataStruct)
  
  if(axis==IAXIS) then
     do k=endPoints(LOW,KAXIS),endPoints(HIGH,KAXIS)
        m=k-endPoints(LOW,KAXIS)+1
        do j=endPoints(LOW,JAXIS),endPoints(HIGH,JAXIS)
           n=j-endPoints(LOW,JAXIS)+1
           do var=1,varCount
              solnData(var,strt:fin,j,k)=region(1:bcVecEnd,n,m,var)
           end do
        end do
     end do

  elseif(axis==JAXIS) then
     do k=endPoints(LOW,KAXIS),endPoints(HIGH,KAXIS)
        m=k+1-endPoints(LOW,KAXIS)
        do i=endPoints(LOW,IAXIS),endPoints(HIGH,IAXIS)
           n=i+1-endPoints(LOW,IAXIS)
           do var=1,varCount
              solnData(var,i,strt:fin,k)=region(1:bcVecEnd,n,m,var)
           end do
        end do
     end do
  elseif(axis==KAXIS) then
     do j=endPoints(LOW,JAXIS),endPoints(HIGH,JAXIS)
        m=j+1-endPoints(LOW,JAXIS)
        do i=endPoints(LOW,IAXIS),endPoints(HIGH,IAXIS)
           n=i+1-endPoints(LOW,IAXIS)
           do var=1,varCount
              solnData(var,i,j,strt:fin)=region(1:bcVecEnd,n,m,var)
           end do
        end do
     end do
  end if

  call Grid_releaseBlkPtr(blockID, solnData, gridDataStruct)

end subroutine gr_bcPutRegion

!!****if* source/Grid/GridBoundaryConditions/Chombo/gr_bcPutRegionsMixedGds
!!
!! NAME
!!  gr_bcPutRegionsMixedGds
!!
!! SYNOPSIS
!!
!!  call gr_bcPutRegionsMixedGds(integer(IN) :: gridDataStruct,
!!                               integer(IN)   :: axis,
!!                               integer(IN)   :: secondDir,
!!                               integer(IN)   :: thirdDir,
!!                               integer(IN)   :: endPoints(LOW:HIGH,MDIM),
!!                               integer(IN)   :: regionSize(REGION_DIM),
!!                               real(in),POINTER :: regionC(:,:,:,:),
!!                               real(in),POINTER :: regionFN(:,:,:,:),
!!                               real(in),POINTER :: regionFT1(:,:,:,:),
!!                               real(in),POINTER :: regionFT2(:,:,:,:),
!!                               integer(IN)   :: blockID,
!!                               integer(IN)   :: idest)
!!  
!! DESCRIPTION 
!!  This routine takes pointers for regions for the application of boundary conditions
!!  to all supported data structures (GDSs), which may access to several GDSs in the same
!!  invocation, and makes sure that the data is stored in the appropriate place.
!!  The pointers should have been gotten earlier with a corrresponding
!!  gr_bcGetRegionsMixedGds call.
!!
!! ARGUMENTS 
!!
!!  regionC     - Pointer to a region of cell-centered data.
!!                May point to data in UNK directly, or in a temporary buffer,
!!                depending on direction of axis and implementation.
!!                May be returned nonassociated, i.e. NULL(), depending on the
!!                combination og GDSs requested with the gridDataStruct argument,
!!  regionFN    - Pointer to a region of face-centered data in the normal direction,
!!                i.e., the direction given by axis.
!!                May point to data in global permanent solution storage directly,
!!                or to data in a temporary buffer,
!!                depending on direction of axis and implementation.
!!                May be returned nonassociated, i.e. NULL(), depending on the
!!                combination og GDSs requested with the gridDataStruct argument,
!!  regionFT1   - Pointer to a region of face-centered data in the first transvers,
!!                direction, i.e., the direction given by secondDir.
!!                May point to data in global permanent solution storage directly,
!!                or to data in a temporary buffer,
!!                depending on direction of axis and implementation.
!!                May be returned nonassociated, i.e. NULL(), depending on the
!!                combination og GDSs requested with the gridDataStruct argument,
!!  regionFT2   - Pointer to a region of face-centered data in the first transvers,
!!                direction, i.e., the direction given by thirdDir.
!!                May point to data in global permanent solution storage directly,
!!                or to data in a temporary buffer,
!!                depending on direction of axis and implementation.
!!                May be returned nonassociated, i.e. NULL(), depending on the
!!                combination og GDSs requested with the gridDataStruct argument,
!!
!! NOTES
!!
!!  Not really implemented for Chombo, since our current use of Chombo as a Grid implementation
!!  does not support face variables.
!!
!! SEE ALSO
!!
!!  gr_bcGetRegionsMixedGds
!!***
subroutine gr_bcPutRegionsMixedGds(gridDataStruct,axis,secondDir,thirdDir,endPoints,&
     regionSize,&
     regionC,regionFN,regionFT1,regionFT2,&
     blockID,idest)
  
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  
  integer, intent(in) :: gridDataStruct,axis, secondDir,thirdDir
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: endPoints
  integer,intent(IN) :: regionSize(REGION_DIM)
  real,pointer,dimension(:,:,:,:) :: regionFN, regionFT1, regionFT2, regionC
  integer, intent(in) :: blockID
  integer,intent(IN) :: idest

  logical :: doCenter, doFaces


  doCenter = (gridDataStruct==CENTER_FACES .OR. gridDataStruct==CENTER)
#if NFACE_VARS>0
  doFaces = (gridDataStruct==CENTER_FACES .OR. gridDataStruct==FACES)
#else
  doFaces = .FALSE.
#endif

  if (doCenter .OR. doFaces) then
     call Driver_abortFlash('gr_gcGetRegionsMixedGds: invalid call when using Chombo as Grid.')
  end if

  return
end subroutine gr_bcPutRegionsMixedGds
