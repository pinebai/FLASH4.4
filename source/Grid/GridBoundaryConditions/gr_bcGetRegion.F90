!!****if* source/Grid/GridBoundaryConditions/gr_bcGetRegion
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
!!  regionSize     : regionSize(BC_DIR) contains the size of the each row
!!                   regionSize(SECOND_DIR) contains the number of rows along the
!!                   second direction, and regionSize(THIRD_DIR) has the number of rows
!!                   along the third direction. regionSize(STRUCTSIZE) contains the
!!                   number of variables in the data structure
!!  mask           : Mask to be used if selected variables are getting boundary
!!                   filled. Currently this has meaning for only PM3 and PM4.
!!  region         : the extracted region
!!  blockID        : the local block ID.
!!                   (With Paramesh3 f. this may actually be a blockHandle that refers
!!                   to a remote block, but this implementation does not actually use
!!                   the blockID at all if the grid is Paramesh3 f. - see idest instead).
!!  idest          : has meaning only for PM3 and PM4, where it distinguishes between 
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

#ifdef FLASH_GRID_PARAMESH3OR4
!!REORDER(5): unk1, facevar1[xyz]
#endif
#ifdef FLASH_GRID_UG
!!REORDER(5): unk, facevar[xyz]
#endif
#ifdef FLASH_GRID_PARAMESH2
!!REORDER(5) : unk
#endif

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_bcGetRegion(gridDataStruct,axis,endPoints,regionSize,mask,&
     region,blockID,idest)
  
#include "constants.h"
  
  use Driver_interface, ONLY : Driver_abortFlash
  
#ifdef FLASH_GRID_UG
  use physicaldata, ONLY: unk,facevarx,facevary,facevarz
#endif
#ifdef FLASH_GRID_PARAMESH3OR4
  use physicaldata, ONLY: unk1,facevarx1, facevary1, facevarz1, gcell_on_cc, gcell_on_fc
  use workspace, ONLY : work1
#endif

#ifdef FLASH_GRID_PARAMESH2
  use workspace, ONLY : work
  use physicaldata, ONLY : unk
#endif

  implicit none
  
  integer, intent(in) :: gridDataStruct,axis
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: endPoints
  integer,intent(IN) :: regionSize(REGION_DIM)
  logical,dimension(regionSize(STRUCTSIZE)),intent(OUT) :: mask
  real,dimension(regionSize(BC_DIR),regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),regionSize(STRUCTSIZE)),intent(OUT) :: region
  integer, intent(in) :: blockID
  integer,intent(IN) :: idest

  integer :: var,i,j,k,n,m,strt,fin, varCount,bcVecEnd
  logical :: validGridDataStruct


  strt = endPoints(LOW,axis)
  fin  = endPoints(HIGH,axis)
  varCount=regionSize(STRUCTSIZE)
  bcVecEnd=regionSize(BC_DIR)

#ifdef DEBUG_GRID

  validGridDataStruct = .false.
  validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
  
  if(.not.validGridDataStruct) then
     print *, "gr_bcGetRegion: gridDataStruct set to improper value"
     print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ,WORK " // &
          " (defined in constants.h)"
     call Driver_abortFlash("gr_bcGetRegion gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,WORK(see constants.h)")
  end if

  if((gridDataStruct==WORK).and.(varCount/=1)) &
       call Driver_abortFlash("gr_bcGetRegion: varCount must be 1 for work array")

  if((fin-strt+1)/=bcVecEnd)&
       call Driver_abortFlash("gr_bcGetRegion: mismatch between rowSize and the region size")
       

#endif
  
  mask=.true.
#ifdef FLASH_GRID_PARAMESH3OR4
  if (gridDataStruct==CENTER) then
     mask(1:varCount)=gcell_on_cc(1:varCount)
  elseif(gridDataStruct==FACEX) then
     mask(1:varCount)=gcell_on_fc(1,1:varCount)
  elseif(gridDataStruct==FACEY) then
     mask(1:varCount)=gcell_on_fc(2,1:varCount)
  elseif(gridDataStruct==FACEZ) then
     mask(1:varCount)=gcell_on_fc(3,1:varCount)
  end if
#endif
  if(axis==IAXIS) then
     do k=endPoints(LOW,KAXIS),endPoints(HIGH,KAXIS)
        m=k-endPoints(LOW,KAXIS)+1
        do j=endPoints(LOW,JAXIS),endPoints(HIGH,JAXIS)
           n=j-endPoints(LOW,JAXIS)+1
           do var=1,varCount
              select case(gridDataStruct)
#ifdef FLASH_GRID_PARAMESH3OR4
                 !! since PM3 f. insists on using unk1 etc
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk1(var,strt:fin,j,k,idest)
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx1(var,strt:fin,j,k,idest)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary1(var,strt:fin,j,k,idest)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz1(var,strt:fin,j,k,idest)
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work1(strt:fin,j,k,idest)
#else
                 !! this section in play if the grid is UG or PM2
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk(var,strt:fin,j,k,blockID)
#if NFACE_VARS>0
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx(var,strt:fin,j,k,blockID)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary(var,strt:fin,j,k,blockID)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz(var,strt:fin,j,k,blockID)
#endif
#ifdef FLASH_GRID_PARAMESH2
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work(strt:fin,j,k,blockID,1)
#endif
#endif
              end select
           end do
        end do
     end do

  elseif(axis==JAXIS) then
     do k=endPoints(LOW,KAXIS),endPoints(HIGH,KAXIS)
        m=k-endPoints(LOW,KAXIS)+1
        do i=endPoints(LOW,IAXIS),endPoints(HIGH,IAXIS)
           n=i-endPoints(LOW,IAXIS)+1
           do var=1,varCount
              select case(gridDataStruct)
#ifdef FLASH_GRID_PARAMESH3OR4
                 !! since PM3 f. insists on using unk1 etc
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk1(var,i,strt:fin,k,idest)
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx1(var,i,strt:fin,k,idest)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary1(var,i,strt:fin,k,idest)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz1(var,i,strt:fin,k,idest)
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work1(i,strt:fin,k,idest)
#else
                 !! this section in play if the grid is UG or PM2
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk(var,i,strt:fin,k,blockID)
#if NFACE_VARS>0
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx(var,i,strt:fin,k,blockID)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary(var,i,strt:fin,k,blockID)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz(var,i,strt:fin,k,blockID)
#endif
#ifdef FLASH_GRID_PARAMESH2
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work(i,strt:fin,k,blockID,1)
#endif
#endif
              end select
           end do
        end do
     end do
  elseif(axis==KAXIS) then
     do j=endPoints(LOW,JAXIS),endPoints(HIGH,JAXIS)
        m=j-endPoints(LOW,JAXIS)+1
        do i=endPoints(LOW,IAXIS),endPoints(HIGH,IAXIS)
           n=i-endPoints(LOW,IAXIS)+1
           do var=1,varCount
              select case(gridDataStruct)
#ifdef FLASH_GRID_PARAMESH3OR4
                 !! since PM3 f. insists on using unk1 etc
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk1(var,i,j,strt:fin,idest)
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx1(var,i,j,strt:fin,idest)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary1(var,i,j,strt:fin,idest)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz1(var,i,j,strt:fin,idest)
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work1(i,j,strt:fin,idest)
#else
                 !! this section in play if the grid is UG or PM2
              case(CENTER)
                 region(1:bcVecEnd,n,m,var)=unk(var,i,j,strt:fin,blockID)
#if NFACE_VARS>0
              case(FACEX)
                 region(1:bcVecEnd,n,m,var)=facevarx(var,i,j,strt:fin,blockID)
              case(FACEY)
                 region(1:bcVecEnd,n,m,var)=facevary(var,i,j,strt:fin,blockID)
              case(FACEZ)
                 region(1:bcVecEnd,n,m,var)=facevarz(var,i,j,strt:fin,blockID)
#endif
#ifdef FLASH_GRID_PARAMESH2
              case(WORK)
                 region(1:bcVecEnd,n,m,varCount)=work(i,j,strt:fin,blockID,1)
#endif
#endif
              end select
           end do
        end do
     end do
  end if
  return
end subroutine gr_bcGetRegion


!!****if* source/Grid/GridBoundaryConditions/gr_bcGetRegionsMixedGds
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
!!                combination of GDSs requested with the gridDataStruct argument,
!!  regionFN    - Pointer to a region of face-centered data in the normal direction,
!!                i.e., the direction given by axis.
!!                May point to data in global permanent solution storage directly,
!!                or to data in a temporary buffer,
!!                depending on direction of axis and implementation.
!!                May be returned nonassociated, i.e. NULL(), depending on the
!!                combination of GDSs requested with the gridDataStruct argument,
!!  regionFT1   - Pointer to a region of face-centered data in the first transverse
!!                direction, i.e., the direction given by secondDir.
!!                May point to data in global permanent solution storage directly,
!!                or to data in a temporary buffer,
!!                depending on direction of axis and implementation.
!!                May be returned nonassociated, i.e. NULL(), depending on the
!!                combination of GDSs requested with the gridDataStruct argument,
!!  regionFT2   - Pointer to a region of face-centered data in the first transverse
!!                direction, i.e., the direction given by thirdDir.
!!                May point to data in global permanent solution storage directly,
!!                or to data in a temporary buffer,
!!                depending on direction of axis and implementation.
!!                May be returned nonassociated, i.e. NULL(), depending on the
!!                combination of GDSs requested with the gridDataStruct argument,
!!***
subroutine gr_bcGetRegionsMixedGds(gridDataStruct,axis,secondDir,thirdDir,endPoints,&
     regionSize,&
     regionC,regionFN,regionFT1,regionFT2,&
     blockID,idest)
  
  use Driver_interface, ONLY : Driver_abortFlash
  
#ifdef FLASH_GRID_UG
  use physicaldata, ONLY: unk,facevarx,facevary,facevarz
#endif
#ifdef FLASH_GRID_PARAMESH3OR4
  use physicaldata, ONLY: unk1,facevarx1, facevary1, facevarz1, gcell_on_cc, gcell_on_fc
  use workspace, ONLY : work1
#endif

#ifdef FLASH_GRID_PARAMESH2
  use workspace, ONLY : work
  use physicaldata, ONLY : unk
#endif

  implicit none
  
  integer, intent(in) :: gridDataStruct,axis, secondDir,thirdDir
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: endPoints
  integer,intent(IN) :: regionSize(REGION_DIM)
  real,pointer,dimension(:,:,:,:) :: regionFN, regionFT1, regionFT2, regionC
  integer, intent(in) :: blockID
  integer,intent(IN) :: idest

  integer,parameter :: ndim=NDIM
  integer :: var,i,j,k,n,m,strt,fin, varCount,bcVecEnd
  integer :: i1,imax,j1,jmax,k1,kmax, nmax,mmax
  logical :: validGridDataStruct
  logical :: doCenter, doFaces

  real,pointer,dimension(:,:,:,:) :: regFN, regFT1, regFT2, regC
  real,pointer,dimension(:,:,:,:) :: pUnk, pFaceVarX,pFaceVarY,pFaceVarZ

  strt = endPoints(LOW,axis)
  fin  = endPoints(HIGH,axis)
  varCount=regionSize(STRUCTSIZE) !This is pretty useless here...
  bcVecEnd=regionSize(BC_DIR)

  nullify(regC);nullify(regFN);nullify(regFT1);nullify(regFT2)

#ifdef DEBUG_GRID

  validGridDataStruct = .false.
  validGridDataStruct= (gridDataStruct == CENTER_FACES).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACES).or.validGridDataStruct
  
  if(.not.validGridDataStruct) then
     print *, "gr_bcGetRegionsMixedGds: gridDataStruct set to improper value"
     print *, "gridDataStruct must be CENTER_FACES or FACES " // &
          " (defined in constants.h)"
     call Driver_abortFlash("gr_bcGetRegionsMixedGds gridDataStruct must be one of CENTER_FACES or FACES (see constants.h)")
  end if

  if((fin-strt+1)/=bcVecEnd)&
       call Driver_abortFlash("gr_bcGetRegionsMixedGds: mismatch between rowSize and the region size")
       

#endif

  doCenter = (gridDataStruct==CENTER_FACES .OR. gridDataStruct==CENTER)
#if NFACE_VARS>0
  doFaces = (gridDataStruct==CENTER_FACES .OR. gridDataStruct==FACES)
#else
  doFaces = .FALSE.
#endif

#ifdef FLASH_GRID_PARAMESH3OR4
  pUnk => unk1(:,:,:,:,idest)
  pFaceVarX => facevarx1(:,:,:,:,idest)
  pFaceVarY => facevary1(:,:,:,:,idest)
  pFaceVarZ => facevarz1(:,:,:,:,idest)
#else
  pUnk => unk(:,:,:,:,blockID)
#if NFACE_VARS>0
  pFaceVarX => facevarx(:,:,:,:,blockID)
  pFaceVarY => facevary(:,:,:,:,blockID)
  pFaceVarZ => facevarz(:,:,:,:,blockID)
#endif
#endif

  if(axis==IAXIS) then
     k1  = endPoints(LOW,KAXIS)
     kmax= endPoints(HIGH,KAXIS)
!     m1  = 1
     mmax=kmax-k1+1
     j1  = endPoints(LOW,JAXIS)
     jmax= endPoints(HIGH,JAXIS)
!     n1  = 1
     nmax=jmax-j1+1
     
#ifdef INDEXREORDER
     if (doCenter) regC  =>     pUnk(1:NUNK_VARS,strt:fin,  j1:jmax,k1:kmax)
     if (doFaces) then
        regFN =>pFaceVarX(1:NFACE_VARS,strt:fin+1,j1:jmax,k1:kmax)
        if(ndim>1) regFT1=>pFaceVarY(1:NFACE_VARS,strt:fin,j1:jmax+1,k1:kmax)
        if(ndim>2) regFT2=>pFaceVarZ(1:NFACE_VARS,strt:fin,j1:jmax,k1:kmax+1)
     end if
#else        
!!$     print*,'Now will allo regC'
     if (doCenter) then
        allocate( regC  (bcVecEnd,  nmax,mmax,NUNK_VARS) )
!!$     print*,'shape(regC)=',shape(regC)
        do var=1,NUNK_VARS
!!$        print*,'Now will copy var=',var
           regC  (1:bcVecEnd,1:nmax,1:mmax,var)=pUnk(var,strt:fin,  j1:jmax,k1:kmax)
        end do
     end if
     if (doFaces) then
        allocate( regFN (bcVecEnd+1,nmax,mmax,NFACE_VARS) )
        if(ndim>1) allocate( regFT1(bcVecEnd,nmax+1,mmax,NFACE_VARS) )
        if(ndim>2) allocate( regFT2(bcVecEnd,nmax,mmax+1,NFACE_VARS) )
        do var=1,NFACE_VARS
           regFN (1:bcVecEnd+1,1:nmax,1:mmax,var)=pFaceVarX(var,strt:fin+1,j1:jmax,k1:kmax)
           if(ndim>1) regFT1(1:bcVecEnd,1:nmax+1,1:mmax,var)=pFaceVarY(var,strt:fin,j1:jmax+1,k1:kmax)
           if(ndim>2) regFT2(1:bcVecEnd,1:nmax,1:mmax+1,var)=pFaceVarZ(var,strt:fin,j1:jmax,k1:kmax+1)
        end do
     end if
#endif

  elseif(axis==JAXIS) then
     k1  = endPoints(LOW,thirdDir)
     kmax= endPoints(HIGH,thirdDir)
!     m1  = 1
     mmax=kmax-k1+1
     i1  = endPoints(LOW,secondDir)
     imax= endPoints(HIGH,secondDir)
!     n1  = 1
     nmax=imax-i1+1
     
     if (doCenter) then
        allocate( regC  (bcVecEnd,  nmax,mmax,NUNK_VARS) )
        if (secondDir==IAXIS) then
           do var=1,NUNK_VARS
              do k=k1,kmax
                 m=k-k1+1
                 regC(1:bcVecEnd,1:nmax,m,var)=TRANSPOSE(pUnk(var,i1:imax,strt:fin,k))
              end do
           end do
        else
           do var=1,NUNK_VARS
              do k=k1,kmax
                 m=k-k1+1
                 regC  (1:bcVecEnd,  1:nmax, m ,var)=     pUnk(var, k ,strt:fin,  i1:imax)
              end do
           end do
        end if
     end if
     if (doFaces) then
        allocate( regFN (bcVecEnd+1,nmax,mmax,NFACE_VARS) )
        if(ndim>1) allocate( regFT1(bcVecEnd,nmax+1,mmax,NFACE_VARS) )
        if(ndim>2) allocate( regFT2(bcVecEnd,nmax,mmax+1,NFACE_VARS) )
        if (secondDir==IAXIS) then
           do var=1,NFACE_VARS
              do k=k1,kmax
                 m=k-k1+1
                 regFN (1:bcVecEnd+1,1:nmax, m ,var)=TRANSPOSE(pFaceVarY(var,i1:imax,strt:fin+1, k ))
                 if(ndim>1) regFT1(1:bcVecEnd,1:nmax+1, m ,var)=TRANSPOSE(pFaceVarX(var,i1:imax+1,strt:fin, k ))
                 if(ndim>2) regFT2(1:bcVecEnd,1:nmax,   m ,var)=TRANSPOSE(pFaceVarZ(var,i1:imax,strt:fin,   k ))
              end do
              if(ndim>2) regFT2  (1:bcVecEnd,1:nmax,mmax+1,var)=TRANSPOSE(pFaceVarZ(var,i1:imax,strt:fin,kmax+1))
           end do
        else
           do var=1,NFACE_VARS
              do k=k1,kmax
                 m=k-k1+1
                 regFN (1:bcVecEnd+1,1:nmax, m ,var)=pFaceVarY(var, k ,strt:fin+1,i1:imax)
                 if(ndim>1) regFT1(1:bcVecEnd,1:nmax+1, m ,var)=pFaceVarZ(var, k ,strt:fin,i1:imax+1)
                 if(ndim>2) regFT2(1:bcVecEnd,1:nmax,   m ,var)=pFaceVarX(var,   k ,strt:fin,i1:imax)
              end do
              if(ndim>2) regFT2(1:bcVecEnd,1:nmax,  mmax+1,var)=pFaceVarX(var,   kmax+1,strt:fin,i1:imax)
           end do
        end if
     end if

  elseif(axis==KAXIS) then
     j1  = endPoints(LOW,thirdDir)
     jmax= endPoints(HIGH,thirdDir)
!     m1  = 1
     mmax=jmax-j1+1
     i1  = endPoints(LOW,secondDir)
     imax= endPoints(HIGH,secondDir)
!     n1  = 1
     nmax=imax-i1+1

     if (doCenter) then
        allocate( regC  (bcVecEnd,  nmax,mmax,NUNK_VARS) )
        do var=1,NUNK_VARS
           do i=strt,fin
              regC  ( i-strt+1 ,  1:nmax,1:mmax,var)=     pUnk(var,i1:imax,j1:jmax, i   )
           end do
        end do
     end if
     if (doFaces) then
        allocate( regFN (bcVecEnd+1,nmax,mmax,NFACE_VARS) )
        if(ndim>1) allocate( regFT1(bcVecEnd,nmax+1,mmax,NFACE_VARS) )
        if(ndim>2) allocate( regFT2(bcVecEnd,nmax,mmax+1,NFACE_VARS) )
        do var=1,NFACE_VARS
           do i=strt,fin
              regFN ( i-strt+1   ,1:nmax,1:mmax,var)=pFaceVarZ(var,i1:imax,j1:jmax,   i )
              if(ndim>1) regFT1( i-strt+1 ,1:nmax+1,1:mmax,var)=pFaceVarX(var,i1:imax+1,j1:jmax, i )
              if(ndim>2) regFT2( i-strt+1 ,1:nmax,1:mmax+1,var)=pFaceVarY(var,i1:imax,j1:jmax+1, i )
           end do
           regFN (  bcVecEnd+1,1:nmax,1:mmax,var)=pFaceVarZ(var,i1:imax,j1:jmax,     fin+1)
        end do
     end if

  end if

!!$  print*,'Now will ptr assign regionC   => regC'
  regionC   => regC
!!$  print*,'Done wth ptr assign regionC   => regC'
  regionFN  => regFN
  regionFT1 => regFT1
  regionFT2 => regFT2
  return
end subroutine gr_bcGetRegionsMixedGds
