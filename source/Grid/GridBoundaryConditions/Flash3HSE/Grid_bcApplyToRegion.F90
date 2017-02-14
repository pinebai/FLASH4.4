!!****if* source/Grid/GridBoundaryConditions/Flash3HSE/Grid_bcApplyToRegion
!!
!! NAME
!!  Grid_bcApplyToRegion
!!
!! SYNOPSIS
!!
!!  call Grid_bcApplyToRegion(integer(IN)  :: bcType,
!!                            integer(IN)  :: gridDataStruct,
!!                            integer(IN)  :: guard,
!!                            integer(IN)  :: axis,
!!                            integer(IN)  :: face,
!!                            real(INOUT)  :: regionData(:,:,:,:),
!!                            integer(IN)  :: regionSize(:),
!!                            logical(IN)  :: mask(:),
!!                            logical(OUT) :: applied,
!!                            integer(IN)  :: blockHandle,
!!                            integer(IN)  :: secondDir,
!!                            integer(IN)  :: thirdDir,
!!                            integer(IN)  :: endPoints(LOW:HIGH,MDIM),
!!                            integer(IN)  :: blkLimitsGC(LOW:HIGH,MDIM),
!!                   OPTIONAL,integer(IN)  :: idest)
!!
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
!!  One reason why information about direction and variable is
!!  included in this interface is because Velocities need to be
!!  treated specially for REFLECTING boundary conditions. if
!!  axis=IAXIS, then the variable VELX_VAR is treated differently,
!!  same with VELY_VAR if axis=JAXIS and VELZ_VAR if
!!  axis=KAXIS. All supported mesh packages extract the vector passed
!!  in through the argument "dataRow" from the appropriated blocks,
!!  and send it to this routine for boundary calculation. The
!!  PERIODIC boundary is calculated by default when the blocks are
!!  exchanging data with each other.
!!  One possible implementation of this interface passes handling of all other
!!  boundary condition types on to calls of the old style Grid_applyBCEdge,
!!  which is called for each of the variables in turn.  However, this implementation
!!  does not do that; it implements the handling of simple boundary condition types
!!  directly by modifying the regionData where it represents guard cells.
!!
!!  This routine supports simple boundary conditions that are applied strictly
!!  directionally, and have no need for other grid information such as the coordinates etc.
!!  Additional dummy arguments blockHandle, secondDir, thirdDir, endPoints, and blkLimitsGC
!!  are not needed for these simple kinds of BCs, but can be used for BC types that do
!!  need coordinate information etc.
!!  Currently supported boundary conditions include "OUTFLOW", "REFLECTING" and "DIODE".
!!  The "PERIODIC" boundary conditions are automatically applied in the process of filling
!!  the guard cells all over the domain, and therefore do not need to call this routine.
!!
!!  If the user wishes to apply different boundary conditions, they can either
!!  use the interface Grid_bcApplyToRegionSpecialized, or make a copy of this
!!  routine in their simulation directory and customize it.
!!
!!
!! ARGUMENTS 
!!
!!  bcType - the type of boundary condition being applied.
!!  gridDataStruct - the Grid dataStructure, should be given as
!!                   one of the constants CENTER, FACEX, FACEY, FACEZ
!!                   (or, with some Grid implementations, WORK).
!!  guard -    number of guard cells
!!  axis  - the dimension along which to apply boundary conditions,
!!          can take values of IAXIS, JAXIS and KAXIS
!!  face    -  can take values LOW and HIGH, defined in constants.h,
!!             to indicate whether to apply boundary on lowerface or 
!!             upperface
!!  regionData     : the extracted region from a block of permanent storage of the 
!!                   specified data structure. Its size is given by regionSize.
!!  regionSize     : regionSize(BC_DIR) contains the size of each row of data
!!                   in the regionData array.  For the common case of guard=4,
!!                   regionSize(BC_DIR) will be 8 for cell-centered data structures
!!                   (e.g., when gridDataStruct=CENTER) and either 8 or 9 for face-
!!                   centered data, depending on the direction given by axis.
!!                   regionSize(SECOND_DIR) contains the number of rows along the second
!!                   direction, and regionSize(THIRD_DIR) has the number of rows
!!                   along the third direction. regionSize(GRID_DATASTRUCT) contains the
!!                   number of variables in the data structure.
!!  mask - if present applies boundary conditions to only selected variables
!!  applied - is set true if this routine has handled the given bcType, otherwise it is 
!!            set to false.
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
!!          In other words, you nearly always want to ignore this optional
!!          argument.  As of FLASH 3.1, it is only used internally within the
!!          Grid unit by a Multigrid GridSolver implementation.
!!
!!  NOTES 
!!
!!   (1)      NOTE that the second index of the endPoints and
!!            blkLimitsGC arrays count the (IAXIS, JAXIS, KAXIS)
!!            directions in the usual order, not permuted as in
!!            regionSize.
!!
!!   (2)      The preprocessor symbols appearing in this description
!!            as well as in the dummy argument declarations (i.e.,
!!            all the all-caps token (other than IN and OUT)) are
!!            defined in constants.h.
!!
!!   (3)      This routine is common to all the mesh packages supported.
!!            The mesh packages extract the small arrays relevant to
!!            boundary condition calculations from their Grid data 
!!            structures. 
!!
!!   (4)      If users wish to apply a different boundary condition, 
!!            they should look at routine Grid_bcApplyToRegionSpecialized.
!!            Customization can occur by creating a nontrivial implementation
!!            of Grid_bcApplyToRegionSpecialized in the Simulation unit (preferred)
!!            or by replacing the default implementation of  Grid_bcApplyToRegion.
!!
!!  HISTORY
!!    (generic:)
!!    2007       Initial Grid_bcApplyToRegion logic    - Anshu Dubey
!!    2007       Cleanup and fixes                     - Klaus Weide, Dongwook Lee 
!!    2008       Tweaks to support Multigrid           - Klaus Weide
!!    2012       Implemented AXISYMMETRIC,EQTSYMMETRIC - Petros Tzeferacos
!!    2012       Tweaks for face variables at reflecting boundaries - Klaus Weide
!!    2013       Added support for GRIDBC_ZERO         - Klaus Weide
!!    (support for hydrostatic boundary conditions:)
!!    2007       HSE BC Implementation based on Zingale/Dursi/etc.
!!                                                     - Dean Townsley
!!    2008,2009  some fixes                            - Dean, Chris, Dongwook
!!    Dec-2009   Integrated into more general code     - Klaus Weide
!!***

subroutine Grid_bcApplyToRegion(bcType,gridDataStruct,&
          guard,axis,face,regionData,regionSize,mask,applied,&
     blockHandle,secondDir,thirdDir,endPoints,blkLimitsGC, idest)

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos
  use Grid_interface, ONLY : Grid_getDeltas
  use gr_bcInterface, ONLY : gr_bcMapBcType, gr_hseStep
  use Grid_data, ONLY : gr_meshMe, gr_geometry, gr_dirGeom
  use gr_bcHseData, ONLY : gr_bcHseDirection, gr_bcHseGravConst, HSE_FORWARD, HSE_BACKWARD, HSE_SETTEMP

  implicit none
  
  integer, intent(IN) :: bcType,axis,face,guard,gridDataStruct
  integer,dimension(REGION_DIM),intent(IN) :: regionSize
  real,dimension(regionSize(BC_DIR),&
       regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),&
       regionSize(STRUCTSIZE)),intent(INOUT)::regionData
  logical,intent(IN),dimension(regionSize(STRUCTSIZE)):: mask
  logical, intent(OUT) :: applied
  integer,intent(IN) :: blockHandle
  integer,intent(IN) :: secondDir,thirdDir
  integer,intent(IN),dimension(LOW:HIGH,MDIM) :: endPoints, blkLimitsGC
  integer,intent(IN),OPTIONAL:: idest

  integer :: i,j, k,ivar,je,ke,n,varCount,bcTypeActual
  logical :: isFace
  integer :: sign
  integer :: sizeGC, start, iend, step, direction
  integer :: velVarBcDir, velVarSecondDir, velVarThirdDir
  real, allocatable, dimension(:) :: cellCenterCoord
  real, allocatable, dimension(:) :: sumyRow, yeRow
  real    :: del(MDIM), deltaBcDir
  real, dimension(EOS_NUM) :: eosData
  real, dimension(NSPECIES) :: massFrac

  select case (bcType)
  case(REFLECTING, AXISYMMETRIC, EQTSYMMETRIC, OUTFLOW,DIODE,GRIDBC_MG_EXTRAPOLATE, &
       NEUMANN_INS)
     applied = .TRUE.           !will handle these types of BCs below
  case(GRIDBC_ZERO)
     applied = .TRUE.           !will handle this type of BCs below
  case(DIRICHLET)
     if ((gridDataStruct==CENTER) .or. (gridDataStruct==WORK)) then
        applied = .TRUE.           !will handle this types below
     else
        applied = .FALSE.       !This file does not implement Dirichlet for
        return                  !face variables etc.; RETURN immediately!
     end if
  case(HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL)
     if (gridDataStruct==CENTER) then
        applied = .FALSE.       ! should have been picked up by Flash2HSE
                                ! implementation if included.
        return                  ! RETURN immediately!
     else
        applied = .TRUE.           !will handle these types below (like OUTFLOW)
     end if
  case(HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL)
     if (gridDataStruct==CENTER) then
        if (axis == gr_bcHseDirection) then
           applied = .TRUE.     ! Handled below.
        else
           applied = .FALSE.       !Wrong direction, don't handle it!
           return                  !RETURN immediately.
        end if
     else
        applied = .TRUE.           !will handle these types below (like OUTFLOW)
     end if
  case default
     applied = .FALSE.
     return                     !RETURN immediately!
  end select


  je=regionSize(SECOND_DIR)
  ke=regionSize(THIRD_DIR)
  varCount=regionSize(STRUCTSIZE)
  isFace = (gridDataStruct==FACEX).and.(axis==IAXIS)
  isFace = isFace.or.((gridDataStruct==FACEY).and.(axis==JAXIS))
  isFace = isFace.or.((gridDataStruct==FACEZ).and.(axis==KAXIS))

  
!!  print*,'in applyBcRegion ',varCount,gridDataStruct,guard,axis,face

  do ivar = 1,varCount
     if(mask(ivar)) then
        call gr_bcMapBcType(bcTypeActual,bcType,ivar,gridDataStruct,axis,face,idest)
        sign = 1

        if (bcTypeActual == REFLECTING) then 
           if (gridDataStruct==CENTER) then
#ifdef VELX_VAR
              if ((axis==IAXIS).and.(ivar==VELX_VAR))sign=-1
#endif
#ifdef VELY_VAR
              if((axis==JAXIS).and.(ivar==VELY_VAR))sign=-1
#endif
#ifdef VELZ_VAR
              if((axis==KAXIS).and.(ivar==VELZ_VAR))sign=-1
#endif
#ifdef MAGX_VAR
              if((axis==IAXIS).and.(ivar==MAGX_VAR))sign=-1
#endif
#ifdef MAGY_VAR
              if((axis==JAXIS).and.(ivar==MAGY_VAR))sign=-1
#endif
#ifdef MAGZ_VAR
              if((axis==KAXIS).and.(ivar==MAGZ_VAR))sign=-1
#endif
           else if (gridDataStruct == FACEX .OR. &
                    gridDataStruct == FACEY .OR. &
                    gridDataStruct == FACEZ ) then
              if (isFace) then
#ifdef MAG_FACE_VAR
                 if(ivar==MAG_FACE_VAR ) sign = -1 
#endif
#ifdef MAGI_FACE_VAR
                 if(ivar==MAGI_FACE_VAR) sign = -1 
#endif
              end if
           end if
        
        else if (bcTypeActual == AXISYMMETRIC) then 
           if (gridDataStruct==CENTER) then
#ifdef VELX_VAR
              if ((axis==IAXIS).and.(ivar==VELX_VAR))sign=-1
#endif
#ifdef VELY_VAR
              if(ivar==VELY_VAR) then
                 if(axis==JAXIS) then
                    sign=-1
                 else ! axis==IAXIS or KAXIS
                    if (gr_dirGeom(JAXIS)==PHI_CYL)sign=-1
                 end if
              end if
#endif
#ifdef VELZ_VAR
              if(ivar==VELZ_VAR) then
                 if(axis==KAXIS) then
                    sign=-1
                 else ! axis==IAXIS or JAXIS 
                    if (gr_dirGeom(KAXIS)==PHI_CYL)sign=-1
                 end if
              end if
#endif
#ifdef MAGX_VAR
              if((axis==IAXIS).and.(ivar==MAGX_VAR))sign=-1
#endif
#ifdef MAGY_VAR
              if(ivar==MAGY_VAR) then
                 if(axis==JAXIS) then
                    sign=-1
                 else ! axis==IAXIS or KAXIS
                    if (gr_dirGeom(JAXIS)==PHI_CYL)sign=-1
                 end if
              end if
#endif
#ifdef MAGZ_VAR
              if(ivar==MAGZ_VAR) then
                 if(axis==KAXIS) then
                    sign=-1
                 else ! axis==IAXIS or JAXIS
                    if (gr_dirGeom(KAXIS)==PHI_CYL)sign=-1
                 end if
              end if
#endif
           else if (gridDataStruct == FACEX .OR. &
                    gridDataStruct == FACEY .OR. &
                    gridDataStruct == FACEZ ) then
              if (isFace) then
#ifdef MAG_FACE_VAR
                 if(ivar==MAG_FACE_VAR ) sign = -1 
#endif
#ifdef MAGI_FACE_VAR
                 if(ivar==MAGI_FACE_VAR) sign = -1 
#endif
              else if((gridDataStruct == FACEY .AND. &
                       gr_dirGeom(JAXIS)==PHI_CYL) .OR. &
                      (gridDataStruct == FACEZ .AND. &
                       gr_dirGeom(KAXIS)==PHI_CYL)) then
#ifdef MAG_FACE_VAR
                 if(ivar==MAG_FACE_VAR ) sign = -1 
#endif
#ifdef MAGI_FACE_VAR
                 if(ivar==MAGI_FACE_VAR) sign = -1 
#endif
              endif
           endif

        else if (bcTypeActual == EQTSYMMETRIC) then 
           if (gridDataStruct==CENTER) then
#ifdef VELX_VAR
              if ((axis==IAXIS).and.(ivar==VELX_VAR))sign=-1
#endif
#ifdef VELY_VAR
              if((axis==JAXIS).and.(ivar==VELY_VAR))sign=-1
#endif
#ifdef VELZ_VAR
              if((axis==KAXIS).and.(ivar==VELZ_VAR))sign=-1
#endif
#ifdef MAGX_VAR
              if((axis==JAXIS).and.(ivar==MAGX_VAR))sign=-1
              if((axis==KAXIS).and.(ivar==MAGX_VAR))sign=-1
#endif
#ifdef MAGY_VAR
              if((axis==IAXIS).and.(ivar==MAGY_VAR))sign=-1
              if((axis==KAXIS).and.(ivar==MAGY_VAR))sign=-1
#endif
#ifdef MAGZ_VAR
              if((axis==IAXIS).and.(ivar==MAGZ_VAR))sign=-1
              if((axis==JAXIS).and.(ivar==MAGZ_VAR))sign=-1
#endif
           else if (gridDataStruct == FACEX .OR. &
                    gridDataStruct == FACEY .OR. &
                    gridDataStruct == FACEZ ) then
              if (.NOT.isFace) then
#ifdef MAG_FACE_VAR
                 if(ivar==MAG_FACE_VAR ) sign = -1 
#endif
#ifdef MAGI_FACE_VAR
                 if(ivar==MAGI_FACE_VAR) sign = -1 
#endif
              end if
            end if

        else ! other bcTypeActual, including DIODE
           if (gridDataStruct==CENTER) then
#ifdef VELX_VAR
              if ((axis==IAXIS).and.(ivar==VELX_VAR))sign=-1
#endif
#ifdef VELY_VAR
              if((axis==JAXIS).and.(ivar==VELY_VAR))sign=-1
#endif
#ifdef VELZ_VAR
              if((axis==KAXIS).and.(ivar==VELZ_VAR))sign=-1
#endif
            end if
        end if

        
        if(face==LOW) then
           select case (bcTypeActual)
           case(REFLECTING)
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)*real(sign)
              end do

           case(AXISYMMETRIC)
              if (gr_geometry == CARTESIAN) call Driver_abortFlash("AXISYMMETRIC boundary only works with curvilinear coordinates")
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)*real(sign)
              end do

           case(EQTSYMMETRIC)
              if (gr_geometry == CARTESIAN) call Driver_abortFlash("EQTSYMMETRIC boundary only works with curvilinear coordinates")
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)*real(sign)
              end do
              
           case(DIRICHLET)
              do i = 1,guard
                 regionData(guard+1-i,1:je,1:ke,ivar)= (1-2*i)*regionData(guard+1,1:je,1:ke,ivar)
              end do
           case(GRIDBC_MG_EXTRAPOLATE)
              do i = 1,guard
                 regionData(guard+1-i,1:je,1:ke,ivar)= (1+i)*regionData(guard+1,1:je,1:ke,ivar) &
                      - i*regionData(guard+2,1:je,1:ke,ivar)
              end do

           case(GRIDBC_ZERO)
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= 0.0
              end do

           case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
                        HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL, &
                        NEUMANN_INS)
!!              print*,'since face was low',je,ke,ivar
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(guard+1,1:je,1:ke,ivar)
              end do
           case(DIODE)
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(guard+1,1:je,1:ke,ivar)
                 if (sign == -1) then
                    do n=1,ke
                       do j=1,je
                          regionData(i,j,n,ivar) = min(regionData(i,j,n,ivar),0.0)
                       end do
                    end do
                 end if
              end do
           case default
              print*,'boundary is',bcType
              call Driver_abortFlash("unsupported boundary condition on Lower Face")
           end select
           
        else  !(face==HIGH)
           
           select case (bcTypeActual)
           case(REFLECTING)
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)*real(sign)
              end do

           case(AXISYMMETRIC)
              if (gr_geometry == CARTESIAN) call Driver_abortFlash("AXISYMMETRIC boundary only works with curvilinear coordinates")
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)*real(sign)
              end do

           case(EQTSYMMETRIC)
              if (gr_geometry == CARTESIAN) call Driver_abortFlash("EQTSYMMETRIC boundary only works with curvilinear coordinates")
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)*real(sign)
              end do

           case(DIRICHLET)
              k=guard
              do i = 1,guard
                 regionData(k+i,1:je,1:ke,ivar)= (1-2*i)*regionData(k,1:je,1:ke,ivar)
              end do
           case(GRIDBC_MG_EXTRAPOLATE)
              k=guard
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(k+i,1:je,1:ke,ivar)= (1+i)*regionData(k,1:je,1:ke,ivar) &
                      - i*regionData(k-1,1:je,1:ke,ivar)
              end do

           case(GRIDBC_ZERO)
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= 0.0
              end do

           case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
                        HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL, &
                        NEUMANN_INS)
              k=guard
              if(isFace)k=k+1
!!              print*,'since face was high',k,je,ke,ivar
              do i = 1,guard
                 regionData(k+i,1:je,1:ke,ivar)= regionData(k,1:je,1:ke,ivar)
              end do
           case(DIODE)
              k=guard
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(k+i,1:je,1:ke,ivar)= regionData(k,1:je,1:ke,ivar)
                 if (sign == -1) then
                    do n = 1,ke
                       do j = 1,je
                          regionData(k+i,j,n,ivar) = max(regionData(k+i,j,n,ivar),0.0)
                       end do
                    end do
                 end if
              end do
           case default
              print*,'boundary is',bcType
              call Driver_abortFlash("unsupported boundary condition on Upper Face")
           end select
        end if
     end if
  end do

  ! Now we are done, unless we are handling Flash3-style HSE BCs.
  
  if (.NOT.((gridDataStruct==CENTER) .AND. &
       (bcType==HYDROSTATIC_NVOUT .OR. bcType==HYDROSTATIC_NVDIODE .OR. bcType==HYDROSTATIC_NVREFL))) then
     return                     ! RETURN in most cases - execution ends here!!
  end if

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!! FROM THIS POINT ON, CODE DEALS WITH FLASH3-STYLE HYDROSTATIC !!!!!!!!!!
  ! !!!!!!!!!!!! BOUNDARY CONDITIONS ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! So we are handling Flash3-style HSE BCs for UNK.
  ! We can assume at this point that (axis == gr_bcHseDirection) is .TRUE.
  select case (axis)
  case(IAXIS)
     velVarBcDir =     VELX_VAR
     velVarSecondDir = VELY_VAR
     velVarThirdDir =  VELZ_VAR
  case(JAXIS)
     velVarBcDir =     VELY_VAR
     velVarSecondDir = VELX_VAR
     velVarThirdDir =  VELZ_VAR
  case(KAXIS)
     velVarBcDir =     VELZ_VAR
     velVarSecondDir = VELX_VAR
     velVarThirdDir =  VELY_VAR
  end select

#if(1)
  ! get coordinates
  sizeGC = blkLimitsGC(HIGH,axis)
  allocate(cellCenterCoord(sizeGC))
  call gr_extendedGetCellCoords(axis, blockHandle, gr_meshMe, CENTER, .true., cellCenterCoord, sizeGC)
  ! assume this is uniform
  deltaBcDir = cellCenterCoord(2)-cellCenterCoord(1)
  deallocate(cellCenterCoord)
#else
  call Grid_getDeltas(blockHandle,del)
  deltaBcDir = del(axis)
#endif
  allocate(sumyRow(regionSize(BC_DIR)))
  allocate(yeRow(regionSize(BC_DIR)))

  do k = 1, regionSize(THIRD_DIR)
     do j = 1, regionSize(SECOND_DIR)

        !-------------------
        !  first take care of the velocities
        if (face==HIGH) then
           do i = 1,guard
              ! zero-gradient has been set above for everything
              ! (cannot know what the user has defined that needs to be propagated)
              if (bcType == HYDROSTATIC_NVREFL) then
                 regionData(regionSize(BC_DIR)-guard+i,j,k,velVarBcDir) = -regionData(regionSize(BC_DIR)-guard+1-i,j,k,velVarBcDir)
                 regionData(regionSize(BC_DIR)-guard+i,j,k,velVarSecondDir) =  &
                      regionData(regionSize(BC_DIR)-guard+1-i,j,k,velVarSecondDir)
                 regionData(regionSize(BC_DIR)-guard+i,j,k,velVarThirdDir) = &
                      regionData(regionSize(BC_DIR)-guard+1-i,j,k,velVarThirdDir)
              else if (bcType == HYDROSTATIC_NVDIODE) then
                 regionData(regionSize(BC_DIR)-guard+i,j,k,velVarBcDir) = &
                      max(0.0,regionData(regionSize(BC_DIR)-guard,j,k,velVarBcDir))
                 ! Velocity components velVarSecondDir and velVarThirdDir already set for OUTFLOW above
!!$              else if (bcType == HYDROSTATIC_NVOUT) then
                 ! All velocity components already set for OUTFLOW above
              endif
          enddo
        else if (face==LOW) then
           do i = 1,guard
              ! zero-gradient has been set above for everything
              ! (cannot know what the user has defined that needs to be propagated)
              if (bcType == HYDROSTATIC_NVREFL) then
                 regionData(i,j,k,velVarBcDir) = -regionData(2*guard+1-i,j,k,velVarBcDir)
                 regionData(i,j,k,velVarSecondDir) =  regionData(2*guard+1-i,j,k,velVarSecondDir)
                 regionData(i,j,k,velVarThirdDir) =  regionData(2*guard+1-i,j,k,velVarThirdDir)
              else if (bcType == HYDROSTATIC_NVDIODE) then
                 regionData(i,j,k,velVarBcDir) = min(0.0,regionData(guard+1,j,k,velVarBcDir))
                 ! Velocity components velVarSecondDir and velVarThirdDir already set for OUTFLOW above
!!$              else if (bcType == HYDROSTATIC_NVOUT) then
                 ! All velocity components already set for OUTFLOW above
              endif
          enddo
        endif
       
        if (mask(DENS_VAR) .or. mask(TEMP_VAR) .or. mask(PRES_VAR) &
             .or. mask(EINT_VAR) .or. mask(GAME_VAR) .or. mask(GAMC_VAR) &
             .or. mask(ENER_VAR)) then
           if (face==HIGH) then
              ! fill stuff out
              start = regionSize(BC_DIR)-guard+1
              iend  = regionSize(BC_DIR)
              step  = 1
              direction = HSE_FORWARD
           else
              start = guard
              iend  = 1
              step  = -1
              direction = HSE_BACKWARD
           endif
#ifdef SUMY_MSCALAR
           sumyRow(:) = regionData(:,j,k,SUMY_MSCALAR)
#else
           sumyRow(:) = 1.0
#endif
#ifdef YE_MSCALAR
           yeRow(:) = regionData(:,j,k,YE_MSCALAR)
#else
           yeRow(:) = 1.0
#endif
           do i = start, iend, step
              massFrac(1:NSPECIES) = regionData(i,j,k,SPECIES_BEGIN:SPECIES_END)
              ! density guess value
              regionData(i,j,k,DENS_VAR) = regionData(i-step,j,k,DENS_VAR)
              call gr_hseStep(regionData(:,j,k,DENS_VAR), &
                   regionData(:,j,k,TEMP_VAR), &
                   yeRow,   &
                   sumyRow, &
                   i,gr_bcHseGravConst, deltaBcDir, direction,2, HSE_SETTEMP, &
                   massFrac)

              ! now get all the eos stuff and fill it in
              eosData(EOS_DENS) = regionData(i,j,k,DENS_VAR)
              eosData(EOS_TEMP) = regionData(i,j,k,TEMP_VAR)
#ifdef SUMY_MSCALAR
              eosData(EOS_ABAR) = 1.0/regionData(i,j,k,SUMY_MSCALAR)
#ifdef YE_MSCALAR
              eosData(EOS_ZBAR) = regionData(i,j,k,YE_MSCALAR)*eosData(EOS_ABAR)
#endif
#endif
              call Eos(MODE_DENS_TEMP, 1, eosData, massFrac)
              regionData(i,j,k,PRES_VAR) = eosData(EOS_PRES)
              regionData(i,j,k,EINT_VAR) = eosData(EOS_EINT)
              regionData(i,j,k,GAME_VAR) = eosData(EOS_PRES)/(eosData(EOS_EINT)*eosData(EOS_DENS)) +1.0
              regionData(i,j,k,GAMC_VAR) = eosData(EOS_GAMC)

              ! Dongwook put in the previously missing factor 1/2 for kinetic energy in October 2009.
              regionData(i,j,k,ENER_VAR) = eosData(EOS_EINT) + 0.5*(regionData(i,j,k,VELX_VAR)**2 &
                   + regionData(i,j,k,VELY_VAR)**2 + regionData(i,j,k,VELZ_VAR)**2)

           enddo
        end if ! mask

     enddo ! second coord
  enddo ! third coord

  deallocate(sumyRow)
  deallocate(yeRow)

  return
end subroutine Grid_bcApplyToRegion
