!!****if* source/Simulation/SimulationMain/HydroStatic/Grid_bcApplyToRegionSpecialized
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

  use Simulation_data, ONLY : sim_gamma, sim_smallX, sim_xyzRef, sim_presRef, &
       sim_gravVector, sim_gravDirec, sim_gravConst, &
       sim_xyzRef, sim_densRef, sim_tempRef, &
       sim_molarMass, sim_gasconstant
  use Grid_data, ONLY : gr_myPE

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

  integer :: je,ke,varCount
  logical :: isFaceVarNormalDir
  integer :: sign

  integer :: offsets(MDIM)
  integer :: i1, i2, i3, ia, ib, ja, jb, ka, kb

  integer :: i, j, k, ivar, n
  real,allocatable,dimension(:) :: cellCenterSweepCoord, secondCoord,thirdCoord
  integer :: istrt, iend, jOffset, kOffset
  real :: xx, yy, zz, gh
  real, allocatable, dimension(:) :: xCoord
  real, allocatable, dimension(:) :: yCoord
  real, allocatable, dimension(:) :: zCoord
  real :: rhoZone, presZone, enerZone, ekinZone, tempZone

  integer :: sizeGC

! ////////////////////////////////////////////////////////////

  je=regionSize(SECOND_DIR)
  ke=regionSize(THIRD_DIR)
  varCount=regionSize(STRUCTSIZE)
  ! the following may be needed by facevar logic
  isFaceVarNormalDir = (gridDataStruct==FACEX).and.(axis==IAXIS)
  isFaceVarNormalDir = isFaceVarNormalDir.or.((gridDataStruct==FACEY).and.(axis==JAXIS))
  isFaceVarNormalDir = isFaceVarNormalDir.or.((gridDataStruct==FACEZ).and.(axis==KAXIS))


  select case (bcType)
  case (REFLECTING,OUTFLOW,DIODE)
     applied = .FALSE.          !Leave it to Grid_bcApplyToRegion default handling - KW
!!! The following part is copied from Grid/GridBoundaryConditions/Flash2HSE/Grid_bcApplyToRegionSpecialized
!!! so that this implementation could handle BOTH User-defined and several variants of hydrostatic BCs!

!!! BEGIN COMMENTED-OUT CASES - KW 2009-12-11 !!!
#if(0)
  case (HYDROSTATIC_F2,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVREFL)
     if (gridDataStruct==CENTER) then

        sizeGC = blkLimitsGC(HIGH,axis)
        allocate(cellCenterSweepCoord(sizeGC))
        call gr_extendedGetCellCoords(axis, blockHandle, gr_myPE, CENTER, .true., cellCenterSweepCoord, sizeGC)
        sizeGC = blkLimitsGC(HIGH,secondDir)
        allocate(secondCoord(sizeGC))
        if (NDIM > 1) &
           call gr_extendedGetCellCoords(secondDir, blockHandle, gr_myPE, CENTER, .true., secondCoord, sizeGC)
        sizeGC = blkLimitsGC(HIGH,thirdDir)
        allocate(thirdCoord(sizeGC))
        if (NDIM > 2) &
           call gr_extendedGetCellCoords(thirdDir, blockHandle, gr_myPE, CENTER, .true., thirdCoord, sizeGC)

        jOffset = endPoints(LOW,secondDir) - 1
        kOffset = endPoints(LOW,thirdDir) - 1
        istrt = endPoints(LOW,axis)
        iend = endPoints(HIGH,axis)
        do k = 1,ke
           do j = 1,je
              call gr_applyFlash2HSEBC(bcType,axis,guard,regionData(:,j,k,:),face,&
                   cellCenterSweepCoord(istrt:iend), &
                   secondCoord(j+jOffset),thirdCoord(k+kOffset))
           end do
        end do
        
        deallocate(cellCenterSweepCoord)
        deallocate(secondCoord)
        deallocate(thirdCoord)

        applied = .TRUE.           !For now this just says we will try to apply - KW
        return
     else
        applied = .FALSE.           ! Face variables, do not do anything here.
     end if
  case (HYDROSTATIC,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVOUT,HYDROSTATIC_NVREFL)
     if (gridDataStruct==CENTER) then

        sizeGC = blkLimitsGC(HIGH,axis)
        allocate(cellCenterSweepCoord(sizeGC))
        call gr_extendedGetCellCoords(axis, blockHandle, gr_myPE, CENTER, .true., cellCenterSweepCoord, sizeGC)
        sizeGC = blkLimitsGC(HIGH,secondDir)
        allocate(secondCoord(sizeGC))
        if (NDIM > 1) &
           call gr_extendedGetCellCoords(secondDir, blockHandle, gr_myPE, CENTER, .true., secondCoord, sizeGC)
        sizeGC = blkLimitsGC(HIGH,thirdDir)
        allocate(thirdCoord(sizeGC))
        if (NDIM > 2) &
           call gr_extendedGetCellCoords(thirdDir, blockHandle, gr_myPE, CENTER, .true., thirdCoord, sizeGC)

        jOffset = endPoints(LOW,secondDir) - 1
        kOffset = endPoints(LOW,thirdDir) - 1
        istrt = endPoints(LOW,axis)
        iend = endPoints(HIGH,axis)
        do k = 1,ke
           do j = 1,je
              call gr_applyFlash3HSEBC(bcType,axis,guard,regionData(:,j,k,:),face,&
                   cellCenterSweepCoord(istrt:iend), &
                   secondCoord(j+jOffset),thirdCoord(k+kOffset))
           end do
        end do
        
        deallocate(cellCenterSweepCoord)
        deallocate(secondCoord)
        deallocate(thirdCoord)

        applied = .TRUE.           !For now this just says we will try to apply - KW
        return
     else
        applied = .FALSE.           ! Face variables, do not do anything here.
     end if
#endif
!!! END COMMENTED-OUT CASES - KW 2009-12-11 !!!
  case (USER_DEFINED)
     applied = .TRUE.           !For now this just says we will try to apply - KW
  case default
     applied = .FALSE.          !Leave it to Grid_bcApplyToRegion default handling - KW
  end select

  if (.NOT. applied) return


! END OF PRELIMINARY STUFF...





! There's a bit of nasty index logic here to make the following possible.
! This could probably be done in a much nicer way.  - KW


  offsets(1:MDIM) = endPoints(LOW,1:MDIM) - blkLimitsGC(LOW,1:MDIM)

! Note: in the following,
! i,ia,ib - refer to the IAXIS direction
! j,ja,jb - refer to the JAXIS direction
! k,ka,kb - refer to the KAXIS direction
! i1  -   refers to the normal direction, perpendicular to the surface at which
!         boundary conditions are applied. this can be the IAXIS, JAXIS, or KAXIS
!         direction.
! i2 -    refers to the second direction, that is the first transverse direction.
!         Normally this can be the IAXIS or JAXIS direction.
! i3 -    refers to the third direction, that is the second transverse direction.
!         Normally this can be the JAXIS or KAXIS direction.

  ! Establish ia,ib,ja,jb,ka,kb to represent a range that covers only the guard
  ! cell part of regionData:
  
  ! Loop ranges in IAXIS, JAXIS, KAXIS direction
  ia = endPoints(LOW, IAXIS)
  ib = endPoints(HIGH,IAXIS)
  ja = endPoints(LOW, JAXIS)
  jb = endPoints(HIGH,JAXIS)
  ka = endPoints(LOW, KAXIS)
  kb = endPoints(HIGH,KAXIS)
  ! Limit loop range in the normal direction to only the guard cell part of the region
  if (face==LOW) then
     select case (axis)
     case(IAXIS)
        ib = ia - 1 + guard
     case(JAXIS)
        jb = ja - 1 + guard
     case(KAXIS)
        kb = ka - 1 + guard
     end select
  else
     select case (axis)
     case(IAXIS)
        ia = ib - guard + 1
     case(JAXIS)
        ja = jb - guard + 1
     case(KAXIS)
        ka = kb - guard + 1
     end select
  end if

  select case (gridDataStruct)

  case(CENTER)
     call gr_bcGetCoords_internal !implemented below

     do k = ka,kb
        zz = zCoord(k)
        do j = ja,jb
           yy  = yCoord (j)
           do i = ia,ib
              xx  = xCoord (i)
              select case (axis)
              case(IAXIS)
                 i1 = i - offsets(IAXIS)
              case(JAXIS)
                 i1 = j - offsets(JAXIS)
              case(KAXIS)
                 i1 = k - offsets(KAXIS)
              end select
              select case (secondDir)
              case(IAXIS)
                 i2 = i - offsets(IAXIS)
              case(JAXIS)
                 i2 = j - offsets(JAXIS)
              case(KAXIS)          !shouldn't...
                 i2 = k - offsets(KAXIS)
              end select
              select case (thirdDir)
              case(IAXIS)          !shouldn't...
                 i3 = i - offsets(IAXIS)
              case(JAXIS)
                 i3 = j - offsets(JAXIS)
              case(KAXIS)
                 i3 = k - offsets(KAXIS)
              end select

! If we get here, we are implementing a USER BC for cell-centered data.
! In this specific USER method, we want to basically fill gourd cells the same
! way in which cells are initialized in Simulation_initBlock.
! Copy code from Simulation_initBlock here, making appropriate changes changes.

              regionData(i1,i2,i3,1:NUNK_VARS)=0.0 ! By default - zero most variables

              ! Multiple species
              !
              do n=SPECIES_BEGIN,SPECIES_END
                 if (n.EQ.SPECIES_BEGIN) then
                    regionData(i1,i2,i3,n)=1.0e0-(NSPECIES-1)*sim_smallX
                 else
                    regionData(i1,i2,i3,n)=sim_smallX
                 end if
              enddo

              !Initialize cells.
              ! Compute the gas energy and set the gamma-values
              ! needed for the equation of  state.
              
              tempZone = sim_tempRef

              gh = sim_gravVector(1) * (xx-sim_xyzRef) + &
                   sim_gravVector(2) * (yy-sim_xyzRef) + &
                   sim_gravVector(3) * (zz-sim_xyzRef)

              !           print*,i,j,exp( sim_molarMass * gh / (sim_gasconstant * tempZone) )
              rhoZone = sim_densRef * exp( sim_molarMass * gh / (sim_gasconstant * tempZone) )
              presZone = sim_presRef * rhoZone / sim_densRef

              ekinZone = 0.0

              enerZone = presZone / (sim_gamma-1.)
              enerZone = enerZone / rhoZone
              enerZone = enerZone + ekinZone

              ! store the variables for the current zone the regionData array

              regionData(i1,i2,i3,DENS_VAR)=rhoZone
              regionData(i1,i2,i3,PRES_VAR)=presZone
              regionData(i1,i2,i3,ENER_VAR)=enerZone
              regionData(i1,i2,i3,EINT_VAR)=enerZone
              regionData(i1,i2,i3,TEMP_VAR)=tempZone
              regionData(i1,i2,i3,GAMC_VAR)=sim_gamma
              regionData(i1,i2,i3,GAME_VAR)=sim_gamma

           enddo
        enddo
     enddo
     call deallocateMem_internal !implemented below

! Now some special handling for velocities.

     do ivar = VELX_VAR,VELZ_VAR !Handle velocities like reflecting
        sign = 1
#ifdef VELX_VAR
        if ((axis==IAXIS).and.(ivar==VELX_VAR))sign=-1
#endif
#ifdef VELY_VAR
        if((axis==JAXIS).and.(ivar==VELY_VAR))sign=-1
#endif
#ifdef VELZ_VAR
        if((axis==KAXIS).and.(ivar==VELZ_VAR))sign=-1
#endif
        if(face==LOW) then
           k = 2*guard+1
           do i = 1,guard
              regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)*real(sign)
           end do
        else
           k = 2*guard+1
           do i = 1,guard
              regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)*real(sign)
           end do
        end if
     end do

! Face variables... do something here.

#if NFACE_VARS > 0
  case(FACEX,FACEY,FACEZ)
     if(face==LOW) then
        do i = 1,guard
           regionData(i,1:je,1:ke,:)= regionData(guard+1,1:je,1:ke,:)
        end do
     else
        k=guard
        if(isFaceVarNormalDir)k=k+1
        do i = 1,guard
           regionData(k+i,1:je,1:ke,:)= regionData(k,1:je,1:ke,:)
        end do
     end if
#endif
  end select


  return
  

! Some details about coordinated initialization - move itnot its own little subroutine.

contains
  subroutine gr_bcGetCoords_internal
    integer :: sizesGC(1:MDIM)


    sizesGC(1:MDIM) = blkLimitsGC(HIGH,1:MDIM)

    allocate(zCoord(sizesGC(KAXIS)))
    if (NDIM == 3) then
       call gr_extendedGetCellCoords(KAXIS, blockHandle,gr_myPE,CENTER,.true., zCoord, sizesGC(KAXIS))
    else
       zCoord = 0.0
    end if

    allocate(yCoord(sizesGC(JAXIS)))
    call gr_extendedGetCellCoords(JAXIS, blockHandle,gr_myPE,CENTER, .true., yCoord,  sizesGC(JAXIS))

    allocate(xCoord(sizesGC(IAXIS)))
    call gr_extendedGetCellCoords(IAXIS, blockHandle,gr_myPE,CENTER,.true., xCoord,  sizesGC(IAXIS))

  end subroutine gr_bcGetCoords_internal

  subroutine deallocateMem_internal
    if(allocated(xCoord)) deallocate(xCoord)
    if(allocated(yCoord)) deallocate(yCoord)
    if(allocated(zCoord)) deallocate(zCoord)
  end subroutine deallocateMem_internal

end subroutine Grid_bcApplyToRegionSpecialized
