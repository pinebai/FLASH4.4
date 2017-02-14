!!****if* source/Grid/GridBoundaryConditions/OneRow/Grid_bcApplyToRegionSpecialized
!!
!!
!! NAME
!!  Grid_bcApplyToRegionSpecialized
!!
!! SYNOPSIS
!!
!!  call Grid_bcApplyToRegionSpecialized(integer(IN)  :: bcType,
!!                                       integer(IN)  :: gridDataStruct,
!!                                       integer(IN)  :: guard,
!!                                       integer(IN)  :: axis,
!!                                       integer(IN)  :: face,
!!                                       real(INOUT)  :: regionData(:,:,:,:),
!!                                       integer(IN)  :: regionSize(:),
!!                                       logical(IN)  :: mask(:),
!!                                       logical(OUT) :: applied,
!!                                       integer(IN)  :: blockHandle,
!!                                       integer(IN)  :: secondDir,
!!                                       integer(IN)  :: thirdDir,
!!                                       integer(IN)  :: endPoints(LOW:HIGH,MDIM),
!!                                       integer(IN)  :: blkLimitsGC(LOW:HIGH,MDIM),
!!                              OPTIONAL,integer(IN)  :: idest )
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Applies the boundary conditions to the specified data structure.
!!  The routine is handed a region that has been extracted from the
!!  data structure, on which it should apply the boundary conditions. 
!!  The direction along which the BC are to be applied is always the first
!!  dimension in the given region, and the last dimension contains the
!!  the variables in the data structure. The middle two dimension contain
!!  the size of the region along the two dimensions of the physical grid
!!  that are not having the BC applied.
!!
!!  This routine applies the boundary conditions on a given face (lowerface
!!  or upperface) along a given axis, by using and setting values
!!  for all variables in the gridDataStruct that are not masked out. The 
!!  argument "mask" has the information about the masked variables.
!! 
!!     If (face=LOW)  
!!       regionData(1:guard,:,:,masked(variables) =  boundary values
!!     If (face=HIGH) 
!!       regionData(regionSize(BC_DIR)-guard+1:regionSize(BC_DIR),:,:,masked(variables) =  boundary values
!!
!!
!!   This interface serves two purposes. One it provides a mechanism by which users can
!!   supply their own custom boundary conditions, and two it allows FLASH to implement
!!   more complex boundary conditions such as those with hydrostatic equilibrium in
!!   isolation from the machinery of setting up the region etc. This interface has
!!   several extra arguments over Grid_bcApplyToRegion. One is the blockHandle, which allows
!!   access to the coordinates information.
!!   This routine is always called first when handling boundary conditions, and if it
!!   finds a match with one of bcTypes that it implements, it sets "applied" to .true., otherwise
!!   "applied" is set to .false. If applied is false, then Grid_bcApplyToRegion is called.
!!   If that routine does not handle the given bcType either, Driver_abortFlash is called.
!!
!!
!! ARGUMENTS 
!!
!! 1. ARGUMENTS SHARED WITH Grid_bcApplyToRegion
!!
!!    bcType - the type of boundary condition being applied.
!!    gridDataStruct - the Grid dataStructure, should be given as
!!                     one of the contants CENTER, FACEX, FACEY, FACEZ.
!!    guard -    number of guardcells
!!    axis  - the dimension along which to apply boundary conditions,
!!            can take values of IAXIS, JAXIS and KAXIS
!!    face    -  can take values LOW and HIGH, defined in constants.h,
!!               to indicate whether to apply boundary on lowerface or 
!!               upperface
!!    regionData     : the extracted region from a block of permanent storage of the 
!!                     specified data structure. Its size is given by regionSize.
!!                     NOTE that the first three dimensions of this array do not necessarily
!!                     correspond to the (IAXIS, JAXIS, KAXIS) directions in this order;
!!                     rather, the axes are permuted such that the first index
!!                     of regionData always corresponds to the direction given by axis.
!!                     See regionSize for more information.
!!    regionSize     : regionSize(BC_DIR) contains the size of each row;
!!                     regionSize(SECOND_DIR) contains the number of rows along the
!!                     second direction, and regionSize(THIRD_DIR) has the number of rows
!!                     along the third direction. (See also below under secondDir,thirdDir
!!                     for the meaning of second and third direction; and see also NOTE (1)
!!                     below.)
!!                     Finally, regionSize(GRID_DATASTRUCT) contains the
!!                     number of variables in the data structure.
!!    mask - if present applies boundary conditions to only selected variables.
!!           However, an implementation of this interface may ignore a mask argument;
!!           a mask should be understood as a possible opportunity for optimization which
!!           an implementation may ignore.
!!    applied - is set true if this routine has handled the given bcType, otherwise it is 
!!              set to false.
!!
!!  idest - Only meaningful with PARAMESH 3 or later.  The argument indicates which slot
!!          in its one-block storage space buffers ("data_1blk.fh") PARAMESH is in the
!!          process of filling.
!!          The following applies when guard cells are filled as part of regular
!!          Grid_fillGuardCells processing (or, in NO_PERMANENT_GUARDCELLS mode,
!!          in order to satisfy a Grid_getBlkPtr request): The value is 1 if guard cells
!!          are being filled in the buffer slot in UNK1 and/or FACEVAR{X,Y,Z}1 or WORK1
!!          that will end up being copied to permanent block data storage (UNK and/or
!!          FACEVAR{X,Y,Z} or WORK, respectively) and/or returned to the user.
!!          The value is 2 if guard cells are being filled in the alternate slot in
!!          the course of assembling data to serve as input for coarse-to-fine
!!          interpolation.
!!          When guard cells are being filled in order to provide input data for
!!          coarse-to-fine interpolation as part of amr_prolong processing (which
!!          is what happens when Grid_updateRefinement is called for an AMR Grid),
!!          the value is always 1.
!!
!!          In other words, an implementation can nearly always ignore this optional
!!          argument.  As of FLASH 3.0, it is only used internally within the
!!          Grid unit and is handled by the GridBoundaryConditions/Grid_bcApplyToRegion
!!          implementation. It is used within the Grid unit by a Multigrid GridSolver
!!          implementation which requires some special handling, but this is only
!!          applied to the WORK data structure.  The argument has been added to the
!!          Grid_bcApplyToRegionSpecialized interface for consistency with
!!          Grid_bcApplyToRegion.
!!
!!
!! 2. ADDITIONAL ARGUMENTS
!!
!!  blockHandle - Handle for the block for which guardcells are to be filled.
!!              In grid implementations other than Paramesh 4, this is always
!!              a local blockID.
!!
!!              With Paramesh 4:
!!              This may be a block actually residing on the local processor,
!!              or the handle may refer to a block that belong to a remote processor
!!              but for which cached information is currently available locally.
!!              The two cases can be distinguished by checking whether 
!!              (blockHandle .LE. lnblocks): this is true only for blocks that
!!              reside on the executing processor.
!!              The block ID is available for passing on to some handlers for 
!!              boundary conditions that may need it, ignored in the default 
!!              implementation.
!!
!!  secondDir,thirdDir -   Second and third coordinate directions.
!!                         These are the transverse directions perpendicular to
!!                         the sweep direction.
!!                         This is not needed for simple boundary condition types
!!                         such as REFLECTIVE or OUTFLOW, It is provided for
!!                         convenience so that more complex boundary condition
!!                         can make use of it.
!!                         The values are currently fully determined by the sweep
!!                         direction bcDir as follows:
!!                          bcDir   |    secondDir       thirdDir
!!                          ------------------------------------------
!!                          IAXIS   |    JAXIS             KAXIS
!!                          JAXIS   |    IAXIS             KAXIS
!!                          KAXIS   |    IAXIS             JAXIS
!!
!!  endPoints - starting and endpoints of the region of interest.
!!              See also NOTE (1) below.
!!
!!  blkLimitsGC - the starting and endpoint of the whole block including
!!                the guard cells, as returned by Grid_getBlkIndexLimits.
!!              See also NOTE (1) below.
!!
!! NOTES
!!
!! (1)        NOTE that the second index of the endPoints and
!!            blkLimitsGC arrays count the (IAXIS, JAXIS, KAXIS)
!!            directions in the usual order, not permuted as in
!!            regionSize.
!!
!! (2)        The preprocessor symbols appearing in this description
!!            as well as in the dummy argument declarations (i.e.,
!!            all the all-caps token (other than IN and OUT)) are
!!            defined in constants.h.
!!
!! (3)        This routine is common to all the mesh packages supported.
!!            The mesh packages extract the small vectors relevant to
!!            boundary conditions calculations from their Grid data 
!!            structures. 
!!
!! (4)   This is an implementation of Grid_bcApplyToRegionSpecialized that
!!       dispatches any boundary conditions other than the simple ones
!!       (REFLECTING,OUTFLOW,DIODE) to the old-style handlers
!!       Grid_applyBCEdgeAllUnkVars and Grid_applyBCEdge, providing
!!       coordinate information as maybe needed by those handlers.
!!
!! (5)   The mask is currently not handled here.
!!
!! SEE ALSO
!!
!!   Grid_bcApplyToRegion            
!!
!!***


subroutine Grid_bcApplyToRegionSpecialized(bcType,gridDataStruct,&
     guard,axis,face,regionData,regionSize,mask,&
     applied,blockHandle,secondDir,thirdDir,endPoints,blkLimitsGC, idest)

#include "constants.h"
#include "Flash.h"

  use Grid_interface, ONLY : Grid_applyBCEdge, Grid_applyBCEdgeAllUnkVars
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(IN) :: bcType,axis,face,guard,gridDataStruct
  integer,intent(IN) :: secondDir,thirdDir
  integer,dimension(REGION_DIM),intent(IN) :: regionSize
  real,dimension(regionSize(BC_DIR),&
       regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),&
       regionSize(STRUCTSIZE)),intent(INOUT)::regionData
  logical,intent(IN),dimension(regionSize(STRUCTSIZE)):: mask
  integer,intent(IN) :: blockHandle
  integer,intent(IN),dimension(LOW:HIGH,MDIM) :: endPoints, blkLimitsGC
  logical, intent(OUT) :: applied
  integer,intent(IN),OPTIONAL:: idest

  integer :: j, k,ivar,je,ke,varCount
  real,allocatable,dimension(:) :: cellCenterSweepCoord, secondCoord,thirdCoord
  integer :: istrt, iend, jOffset, kOffset
  logical :: isFaceVarNormalDir


  !! This is an implementation of Grid_bcApplyToRegionSpecialized that
  !! dispatches any boundary conditions other than the simple ones
  !! (REFLECTING,OUTFLOW,DIODE) to the old-style handlers
  !! Grid_applyBCEdgeAllUnkVars and Grid_applyBCEdge, providing
  !! coordinate information as maybe needed by those handlers. - KW


  je=regionSize(SECOND_DIR)
  ke=regionSize(THIRD_DIR)
  varCount=regionSize(STRUCTSIZE)
  isFaceVarNormalDir = (gridDataStruct==FACEX).and.(axis==IAXIS)
  isFaceVarNormalDir = isFaceVarNormalDir.or.((gridDataStruct==FACEY).and.(axis==JAXIS))
  isFaceVarNormalDir = isFaceVarNormalDir.or.((gridDataStruct==FACEZ).and.(axis==KAXIS))


  select case (bcType)
  case (REFLECTING,OUTFLOW,DIODE,GRIDBC_ZERO)
     applied = .FALSE.          !Leave it to Grid_bcApplyToRegion default handling - KW
  case default
     applied = .TRUE.           !For now this just says we will try to apply - KW
  end select

  if (applied) then


     call gr_bcGetCoords_internal


     jOffset = endPoints(LOW,secondDir) - 1
     kOffset = endPoints(LOW,thirdDir) - 1

     if (gridDataStruct==CENTER) then
        istrt = endPoints(LOW,axis)
        iend = endPoints(HIGH,axis)
        do k = 1,ke
           do j = 1,je
              call Grid_applyBCEdgeAllUnkVars(bcType,axis,guard,regionData(:,j,k,:),face, &
                   cellCenterSweepCoord(istrt:iend), &
                   secondCoord(j+jOffset),thirdCoord(k+kOffset), blockHandle)
           end do
        end do
     else
        do ivar=1,varCount
           do k = 1,ke
              do j = 1,je
                 call Grid_applyBCEdge(bcType,axis,guard,ivar,regionData(:,j,k,ivar),face, &
                      gridDataStruct, blockHandle, secondCoord(j+jOffset),thirdCoord(k+kOffset))
              end do
           end do
        end do
     end if

     call deallocateMem_internal

  end if
     

  return
  
contains
  subroutine gr_bcGetCoords_internal
    use Grid_data, ONLY : gr_meshMe
    integer :: sizeGC

    sizeGC = blkLimitsGC(HIGH,axis)
    allocate(cellCenterSweepCoord(sizeGC))
    if (isFaceVarNormalDir) then
       call gr_extendedGetCellCoords(axis, blockHandle, gr_meshMe, FACES, .true., cellCenterSweepCoord, sizeGC)
    else
       call gr_extendedGetCellCoords(axis, blockHandle, gr_meshMe, CENTER, .true., cellCenterSweepCoord, sizeGC)
    end if
    sizeGC = blkLimitsGC(HIGH,secondDir)
    allocate(secondCoord(sizeGC))
#if NDIM > 1
    if ((gridDataStruct==FACEX .AND. secondDir==IAXIS) .OR. &
         (gridDataStruct==FACEY .AND. secondDir==JAXIS) .OR. &
         (gridDataStruct==FACEZ .AND. secondDir==KAXIS)) then
       call gr_extendedGetCellCoords(secondDir, blockHandle, gr_meshMe, FACES, .true., secondCoord, sizeGC)
    else
       call gr_extendedGetCellCoords(secondDir, blockHandle, gr_meshMe, CENTER, .true., secondCoord, sizeGC)
    end if
#endif
    sizeGC = blkLimitsGC(HIGH,thirdDir)
    allocate(thirdCoord(sizeGC))
#if NDIM > 2
    if ((gridDataStruct==FACEX .AND. thirdDir==IAXIS) .OR. &
         (gridDataStruct==FACEY .AND. thirdDir==JAXIS) .OR. &
         (gridDataStruct==FACEZ .AND. thirdDir==KAXIS)) then
       call gr_extendedGetCellCoords(thirdDir, blockHandle, gr_meshMe, FACES, .true., thirdCoord, sizeGC)
    else
       call gr_extendedGetCellCoords(thirdDir, blockHandle, gr_meshMe, CENTER, .true., thirdCoord, sizeGC)
    end if
#endif

  end subroutine gr_bcGetCoords_internal

  subroutine deallocateMem_internal
    deallocate(cellCenterSweepCoord)
    deallocate(secondCoord)
    deallocate(thirdCoord)
  end subroutine deallocateMem_internal

end subroutine Grid_bcApplyToRegionSpecialized
