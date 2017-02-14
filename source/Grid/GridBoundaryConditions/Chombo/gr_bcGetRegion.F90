!!****if* source/Grid/GridBoundaryConditions/Chombo/gr_bcGetRegion
!!
!! NAME
!!  gr_bcGetRegion
!!
!! SYNOPSIS
!!
!!  call gr_bcGetRegion(integer(IN) :: gridDataStruct,
!!                 integer(IN)   :: axis,
!!                 integer(IN)   :: endPoints(LOW:HIGH,MDIM),
!!                 integer(IN)   :: regionSize(REGION_DIM),
!!                 integer(OUT)  :: mask(regionSize(4)),
!!                 real(out)     :: region(regionSize(1),regionSize(2),regionSize(3),regionSize(4)),
!!                 integer(IN)   :: blockID,
!!                 integer(IN)   :: idest)
!!  
!! DESCRIPTION 
!!  This routine creates a region for the application of boundary condition
!!  to all supported data structures, described for argument 
!!  gridDataStruct below.
!!  The region is stored in a four-dimensional array, where the fourth
!!  dimension represents the individual variables of the data structure.
!!  The other dimensions store a set of rows containing relevant sections of
!!  the block on which the boundary conditions are being applied. Each row
!!  contains a complete set of data points to correctly apply the boundary
!!  conditions along the specified axis. The endPoints argument specifies
!!  the bounding box of the regions being selected. For more details, see
!!  the example below:
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
!!  mask           : Mask to be used if selected variables are getting boundary
!!                   filled. Currently this only has meaning for PARAMESH4.
!!  region         : the extracted region
!!  blockID        : the local block ID.
!!                   (With Paramesh4 f. this may actually be a blockHandle that refers
!!                   to a remote block, but this implementation does not actually use
!!                   the blockID at all if the grid is Paramesh4 f. - see idest instead).
!!  idest          : has meaning only for PARAMESH4, where it distinguishes between 
!!                   leaf and parent nodes, should be 1 or 2; see NOTES below.
!!
!!
!! EXAMPLE 
!!   In this example with 2D data on a LEAF block, 
!!   we want to apply boundary conditions on the right face of IAXIS, 
!!   and we wish to fetch columns 7 and 8 of the interior
!!   data and all the columns of the guardcell data. Since this example has
!!   4 guardcells on each side, the actual column numbers are 11,12 for the
!!   interior and 13-16 for the guardcells to be filled.
!!
!!       ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!     8 ----|-|-|-|-|-|-|*|*|----
!!     7 ----|-|-|-|-|-|-|*|*|----
!!     6 ----|-|-|-|-|-|-|*|*|----
!!     5 ----|-|-|-|-|-|-|*|*|----
!!     4 ----|-|-|-|-|-|-|*|*|----
!!     3 ----|-|-|-|-|-|-|*|*|----
!!     2 ----|-|-|-|-|-|-|*|*|----
!!     1 ----|-|-|-|-|-|-|*|*|----
!!     j ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!       ---- - - - - - - - - ----
!!         i  1-2-3-4 5-6-7-8 
!!     
!!     Then the values in the arguement endpoints should be
!!
!!          endPoint(LOW,IAXIS) = 11; endPoint(HIGH,IAXIS)=16
!!          endPoint(LOW,JAXIS) = 1 ; endPoint(HIGH,JAXIS)=16
!!          endPoint(LOW:HIGH,KAXIS) = 1
!!
!!     RegionSize argument should have 
!!         RegionSize(:)=endPoint(HIGH,:)-endPoint(LOW,:)+1
!!
!!     The argument Region will contain the data being fetched, so
!!     should be allocated as 
!!     Region(RegionSize(IAXIS),RegionSize(JAXIS),RegionSize(KAXIS),vars)
!!     where vars is the number of variables in the data structure;
!!     NUNK_VARS for cell centered, NFACE_VARS for face centered along IAXIS etc.
!!
!!     Please Note that if we were interested in rows (in other words
!!     the top regions) then the allocation would have been
!!     Region(RegionSize(JAXIS),RegionSize(IAXIS),RegionSize(KAXIS),vars)
!!
!!     The call will have the following syntax:
!!
!!     call gr_bcGetRegion(CENTER,IAXIS,endPoint,regionSize,mask,Region,blockID,
!!     idest)
!!
!! NOTES
!!  Beginning with PARAMESH3: The solution data used to fill the
!!  region array are not copied directly from "permanent" storage
!!  (UNK,WORK,etc.), but from the one-block arrays (UNK1,WORK1,etc.)
!!  that are filled by PARAMESH while it is processing a block's data.
!!  Calls to gr_bcGetRegion are therefore only valid in certain
!!  contexts.  The idest argument is then taken as an index into these
!!  one-block arrays. It distinguishes between the slots available and
!!  must match the slot that has actually been filled by PARAMESH.
!!  
!!
!!***
#include "Flash.h"

!!REORDER(4): solnData

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_bcGetRegion(gridDataStruct,axis,endPoints,regionSize,mask,&
     region,blockID,idest)
  
#include "constants.h"
  
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  implicit none
  
  integer, intent(in) :: gridDataStruct,axis
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: endPoints
  integer,intent(IN) :: regionSize(REGION_DIM)
  logical,dimension(regionSize(STRUCTSIZE)),intent(OUT) :: mask
  real,dimension(regionSize(BC_DIR),regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),regionSize(STRUCTSIZE)),intent(OUT) :: region
  integer, intent(in) :: blockID
  integer,intent(IN) :: idest

  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: var,i,j,k,n,m,strt,fin, varCount,bcVecEnd
  logical :: validGridDataStruct


  strt = endPoints(LOW,axis)
  fin  = endPoints(HIGH,axis)
  varCount=regionSize(STRUCTSIZE)
  bcVecEnd=regionSize(BC_DIR)

  !An invalid gridDataStruct is caught by Grid_getBlkPtr so
  !there is no extra checking in gr_bcGetRegion.
  call Grid_getBlkPtr(blockID, solnData, gridDataStruct)

  mask=.true.
  if(axis==IAXIS) then
     do k=endPoints(LOW,KAXIS),endPoints(HIGH,KAXIS)
        m=k-endPoints(LOW,KAXIS)+1
        do j=endPoints(LOW,JAXIS),endPoints(HIGH,JAXIS)
           n=j-endPoints(LOW,JAXIS)+1
           do var=1,varCount
              region(1:bcVecEnd,n,m,var)=solnData(var,strt:fin,j,k)
           end do
        end do
     end do

  elseif(axis==JAXIS) then
     do k=endPoints(LOW,KAXIS),endPoints(HIGH,KAXIS)
        m=k-endPoints(LOW,KAXIS)+1
        do i=endPoints(LOW,IAXIS),endPoints(HIGH,IAXIS)
           n=i-endPoints(LOW,IAXIS)+1
           do var=1,varCount
              region(1:bcVecEnd,n,m,var)=solnData(var,i,strt:fin,k)
           end do
        end do
     end do
  elseif(axis==KAXIS) then
     do j=endPoints(LOW,JAXIS),endPoints(HIGH,JAXIS)
        m=j-endPoints(LOW,JAXIS)+1
        do i=endPoints(LOW,IAXIS),endPoints(HIGH,IAXIS)
           n=i-endPoints(LOW,IAXIS)+1
           do var=1,varCount
              region(1:bcVecEnd,n,m,var)=solnData(var,i,j,strt:fin)
           end do
        end do
     end do
  end if

  call Grid_releaseBlkPtr(blockID, solnData, gridDataStruct)

end subroutine gr_bcGetRegion

!!****if* source/Grid/GridBoundaryConditions/Chombo/gr_bcGetRegionsMixedGds
!!
!! NAME
!!  gr_bcGetRegionsMixedGds
!!
!! SYNOPSIS
!!
!!  call gr_bcGetRegionsMixedGds(integer(IN) :: gridDataStruct,
!!                               integer(IN)   :: axis,
!!                               integer(IN)   :: secondDir,
!!                               integer(IN)   :: thirdDir,
!!                               integer(IN)   :: endPoints(LOW:HIGH,MDIM),
!!                               integer(IN)   :: regionSize(REGION_DIM),
!!                               real(out),POINTER :: regionC(:,:,:,:),
!!                               real(out),POINTER :: regionFN(:,:,:,:),
!!                               real(out),POINTER :: regionFT1(:,:,:,:),
!!                               real(out),POINTER :: regionFT2(:,:,:,:),
!!                               integer(IN)   :: blockID,
!!                               integer(IN)   :: idest)
!!  
!! DESCRIPTION 
!!  This routine returns pointers for regions for the application of boundary conditions
!!  to all supported data structures (GDSs), giving access to several GDSs in the same
!!  invocation.
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
!!***
subroutine gr_bcGetRegionsMixedGds(gridDataStruct,axis,secondDir,thirdDir,endPoints,&
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

  ! Should not be reached - may pacify compiler warnings.
  regionC   => NULL()
  regionFN  => NULL()
  regionFT1 => NULL()
  regionFT2 => NULL()
  return
end subroutine gr_bcGetRegionsMixedGds
