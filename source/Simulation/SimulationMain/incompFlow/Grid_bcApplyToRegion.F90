!!****if* source/Simulation/SimulationMain/incompFlow/Grid_bcApplyToRegion
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
!!   Where masked(variables)
!!     If (face=LOW)  
!!       regionData(1:guard,:,:,variables) =  boundary values
!!     If (face=HIGH) 
!!       regionData(regionSize(BC_DIR)-guard+1:regionSize(BC_DIR),:,:,variables) =  boundary values
!!
!!  One reason why information about direction and variable is
!!  included in this interface is because velocities need to be
!!  treated specially for REFLECTING boundary conditions. if
!!  axis=IAXIS, then the variable VELX_VAR is treated differently,
!!  same with VELY_VAR if axis=JAXIS and VELZ_VAR if
!!  axis=KAXIS. All supported mesh packages extract the vector passed
!!  in through the argument "dataRow" from the appropriated blocks,
!!  and send it to this routine for boundary calculation. The
!!  PERIODIC boundary is calculated by default when the blocks are
!!  exchanging data with each other, and therefore will not be encountered
!!  by this routine.
!!  One possible implementation of this interface passes handling of all other
!!  boundary condition types on to calls of the old style Grid_applyBCEdge,
!!  which is called for each of the variables in turn.  However, this implementation
!!  does not do that; it implements the handling of simple boundary condition types
!!  directly by modifying the regionData where it represents guard cells.
!!
!!  This default version of this routine implements only simple boundary
!!  conditions that are applied strictly directionally and have no need for
!!  other grid information, such as cell coordinates, etc.  Currently
!!  supported simple boundary conditions include "OUTFLOW", "REFLECTING" and
!!  "DIODE".
!!  Additional dummy arguments blockHandle, secondDir, thirdDir, endPoints,
!!  and blkLimitsGC are not needed for these simple kinds of BCs, but can be
!!  used by alternative implementations for BC types that do need coordinate
!!  information, etc.
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
!!  axis  - the direction along which to apply boundary conditions,
!!          can take values of IAXIS, JAXIS and KAXIS
!!  face    -  can take values LOW and HIGH, defined in constants.h,
!!             to indicate whether to apply boundary on lowerface or 
!!             upperface
!!  regionData     : the extracted region from a block of permanent storage of the 
!!                   specified data structure. Its size is given by regionSize.
!!                     NOTE that the first three dimensions of this array do not necessarily
!!                     correspond to the (IAXIS, JAXIS, KAXIS) directions in this order;
!!                     rather, the axes are permuted such that the first index
!!                     of regionData always corresponds to the direction given by axis.
!!                     See regionSize for more information.
!!  regionSize     : regionSize(BC_DIR) contains the size of each row of data
!!                   in the regionData array.  With row we mean here an array slice
!!                   regionData(:,I2,I3,VAR), corresponding to cells that are situated
!!                   along a line in the 'axis' direction. For the common case of guard=4,
!!                   regionSize(BC_DIR) will be 8 for cell-centered data structures
!!                   (e.g., when gridDataStruct=CENTER) and either 8 or 9 for face-
!!                   centered data, depending on the direction given by axis.
!!                   regionSize(SECOND_DIR) contains the number of rows along the second
!!                   direction, and regionSize(THIRD_DIR) has the number of rows
!!                   along the third direction. regionSize(GRID_DATASTRUCT) contains the
!!                   number of variables in the data structure.
!!  mask - if present, boundary conditions are to be applied only to selected variables.
!!         However, an implementation of this interface may ignore the mask argument;
!!         a mask should be understood as a possible opportunity for optimization which
!!         an implementation may ignore.
!!         Specifying a mask does not mean that previous values of other variables in
!!         guard cells will be left undisturbed.
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
!!                         the sweep direction.  SecondDir and thirdDir give
!!                         the meaning of the second and third dimension,
!!                         respectively, of the regionData array.
!!                         This is not needed for simple boundary condition types
!!                         such as REFLECTIVE or OUTFLOW, It is provided for
!!                         convenience so that more complex boundary condition
!!                         can make use of it.
!!                         The values are currently fully determined by the sweep
!!                         direction axis as follows:
!!                          axis   |    secondDir       thirdDir
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
!!   (1)      NOTE that the second indices of the endPoints and
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
!!
!!    2007       Initial Grid_bcApplyToRegion logic    - Anshu Dubey
!!    2007       Cleanup and fixes                     - Klaus Weide, Dongwook Lee 
!!    2008       Tweaks to support Multigrid           - Klaus Weide
!!    2010-2012  Implement bc types for incomp. NS     - Marcos Vanella
!!    2012       Implemented AXISYMMETRIC,EQTSYMMETRIC - Petros Tzeferacos
!!    2012       Tweaks for face variables at reflecting boundaries - Klaus Weide
!!    2013       Added support for GRIDBC_ZERO         - Klaus Weide
!!    2015       Support for GRIDBC_EXTRAPOLATE_NSC    - Klaus Weide
!!
!!***

subroutine Grid_bcApplyToRegion(bcType,gridDataStruct,&
          guard,axis,face,regionData,regionSize,mask,applied,&
     blockHandle,secondDir,thirdDir,endPoints,blkLimitsGC, idest)

#include "constants.h"
#include "Flash.h"
#include "IncompNS.h"

  use Driver_interface, ONLY : Driver_abortFlash
  use gr_bcInterface, ONLY : gr_bcMapBcType
  use Grid_data, ONLY : gr_geometry, gr_dirGeom, &
       gr_smallrho, gr_smallE

!!$  use Grid_interface, ONLY : Grid_getDeltas,         &
!!$                             Grid_getBlkBoundBox,    &
!!$                             Grid_getBlkCenterCoords

  use Driver_data, ONLY : dr_simTime, dr_dt

#if NDIM == 2
  use IncompNS_data, ONLY : ins_invRe,ins_convvel,ins_predcorrflg,ins_alfa,uvel_x,vvel_x,wvel_x, &
                            uvel_y,vvel_y,wvel_y,ins_outflowgridChanged,ins_intschm, &
                            ins_tlevel
#elif NDIM == 3
  use IncompNS_data, ONLY : ins_invRe,ins_convvel,ins_predcorrflg,ins_alfa,uvel_x,vvel_x,wvel_x, &
                            uvel_y,vvel_y,wvel_y,uvel_z,vvel_z,wvel_z, ins_outflowgridChanged,   &
                            ins_rhoa, ins_intschm, ins_tlevel
#endif

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
  integer    :: sign
  real    :: smallP

  integer :: ia,ib,ja,jb,ka,kb

  real, dimension(MDIM)  :: del
!!$  real, dimension(MDIM)  :: coord,bsize
!!$  real ::  boundBox(2,MDIM)

  real :: alfadt
!!$  real :: xcell,xedge,ycell,yedge
!!$ 
!!$  real, save :: TLEVEL
!!$
!!$  integer :: countj


  select case (bcType)
  case(REFLECTING, AXISYMMETRIC, EQTSYMMETRIC, OUTFLOW,DIODE,GRIDBC_MG_EXTRAPOLATE, &
       GRIDBC_EXTRAPOLATE_NSC, &
       DIRICHLET)
     applied = .TRUE.           !will handle these types of BCs below
  case(NEUMANN_INS, NOSLIP_INS, SLIP_INS, INFLOW_INS, MOVLID_INS, OUTFLOW_INS) ! Incompressible solver BCs
     applied = .TRUE.           !will handle these types of BCs below
  case(GRIDBC_ZERO)
     applied = .TRUE.           !will handle this type of BCs below
  case(HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
       HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL)
     if (gridDataStruct==CENTER) then
        applied = .FALSE.       !should be picked up by Flash2HSE implementation
                                !or Flash3HSE implementation if included.
        return                  !RETURN immediately!
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


  !print*,'in applyBcRegion ',gridDataStruct,CENTER !,FACEX,FACEY,FACEZ 
  !print*,'in applyBcRegion ',varCount,guard,axis,face



  do ivar = 1,varCount
     if(mask(ivar)) then
        call gr_bcMapBcType(bcTypeActual,bcType,ivar,gridDataStruct,axis,face,idest)
        sign = 1


!!  Handle sign flip for reflective, axisymmetric and eqtsymmetric boundaries.
!!  Here n stands for normal, while p,n are tangent components. In curvilinear
!!  coords we denote p for poloidal and t for toiroidal
!!
!!   REFLECTIVE
!! Vn -> -Vn,  Bn -> -Bn
!! Vp ->  Vp,  Bp ->  Bp
!! Vt ->  Vt,  Bt ->  Bt
!!
!!   AXISYMMETRIC
!! Vn -> -Vn,  Bn -> -Bn
!! Vp ->  Vp,  Bp ->  Bp
!! Vt -> -Vt,  Bt -> -Bt
!!
!!   EQTSYMMETRIC
!! Vn -> -Vn,  Bn ->  Bn
!! Vp ->  Vp,  Bp -> -Bp
!! Vt ->  Vt,  Bt -> -Bt
 
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
           endif
        
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
           case(GRIDBC_EXTRAPOLATE_NSC)
              do i = 1,guard
                 regionData(guard+1-i,1:je,1:ke,ivar)= (1+i)*regionData(guard+1,1:je,1:ke,ivar) &
                      - i*regionData(guard+2,1:je,1:ke,ivar)
                 where (regionData(guard+1-i,1:je,1:ke,ivar)*regionData(guard+1,1:je,1:ke,ivar) < 0.0)
                    regionData(guard+1-i,1:je,1:ke,ivar)= 0.0
                 end where
#ifdef PRES_VAR
                 if (ivar==PRES_VAR) then
                    smallP = gr_smallRho * gr_smallE
                    regionData(guard+1-i,1:je,1:ke,ivar)= max(smallP,regionData(guard+1-i,1:je,1:ke,ivar))
                 end if
#endif
#ifdef ENER_VAR
                 if (ivar==ENER_VAR) then
                    regionData(guard+1-i,1:je,1:ke,ivar)= max(gr_smallE,regionData(guard+1-i,1:je,1:ke,ivar))
                 end if
#endif
#ifdef EINT_VAR
                 if (ivar==EINT_VAR) then
                    regionData(guard+1-i,1:je,1:ke,ivar)= max(gr_smallE,regionData(guard+1-i,1:je,1:ke,ivar))
                 end if
#endif
#ifdef DENS_VAR
                 if (ivar==DENS_VAR) then
                    regionData(guard+1-i,1:je,1:ke,ivar)= max(gr_smallRho,regionData(guard+1-i,1:je,1:ke,ivar))
                 end if
#endif
              end do

           case(GRIDBC_ZERO)
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= 0.0
              end do

           case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
                        HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL)
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


           ! The following cases correspond to Incompressible solver BC combinations:
           ! Staggered grid, velocities: face centered, pressure: cell centered. 
           case(NEUMANN_INS)  ! Neumann BC

              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
              end do


           case(NOSLIP_INS)

              if (gridDataStruct==WORK) then ! NEUMANN 

              k = 2*guard+1 
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
              end do  

              elseif (gridDataStruct==CENTER) then !NEUMANN BC FOR PRESSURE AND DELTAP
              select case(ivar)
              case (PRES_VAR,DELP_VAR,TVIS_VAR)
              k = 2*guard+1 
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
              end do                          
              end select


              else ! BOUNDARY CONDITIONS ON VELOCITIES               
              if (isFace) then ! Set to zero velocities normal to boundary, up to boundary
                 if(ivar == VELC_FACE_VAR) then
                 k = 2*guard+2
                 do i = 1,guard
                    regionData(i,1:je,1:ke,ivar)= -regionData(k-i,1:je,1:ke,ivar)
                 end do
                 regionData(guard+1,1:je,1:ke,ivar)= 0.
                 endif
              else             ! Use guardcells to set to zero velocities not normal to boundary, at boundary
                 k = 2*guard+1   
                 if(ivar == VELC_FACE_VAR) then                               
                 do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= -regionData(k-i,1:je,1:ke,ivar)
                 end do
                 endif
              endif

              endif


           case(SLIP_INS)


              if (gridDataStruct==WORK) then ! NEUMANN 

              k = 2*guard+1 
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
              end do  

              elseif (gridDataStruct==CENTER) then !NEUMANN BC FOR PRESSURE AND DELTAP
              select case(ivar)
              case (PRES_VAR,DELP_VAR,TVIS_VAR)
              k = 2*guard+1 
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
              end do                          
              end select


              else ! BOUNDARY CONDITIONS ON VELOCITIES               
              if (isFace) then ! Set to zero velocities normal to boundary, up to boundary
                 if((ivar == VELC_FACE_VAR)) then
                 k = 2*guard+2
                 do i = 1,guard
                    regionData(i,1:je,1:ke,ivar)= -regionData(k-i,1:je,1:ke,ivar)
                 end do
                 regionData(guard+1,1:je,1:ke,ivar)= 0.
                 endif
              else             ! Use guardcells to set to zero normal gradients of velocities not normal to boundary, at boundary
                 k = 2*guard+1   
                 if(ivar == VELC_FACE_VAR) then                               
                 do i = 1,guard !-1
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
                 end do
                 !regionData(guard,1:je,1:ke,ivar)= (1.+1./16.)*regionData(k-guard,1:je,1:ke,ivar)-3./16.*regionData(k-guard+1,1:je,1:ke,ivar)+&
                 !                                  3./16.*regionData(k-guard+2,1:je,1:ke,ivar)-1./16.*regionData(k-guard+3,1:je,1:ke,ivar)
                 endif
              endif

              endif


           case(MOVLID_INS)


              if (gridDataStruct==WORK) then ! NEUMANN 

              k = 2*guard+1 
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
              end do  

              elseif (gridDataStruct==CENTER) then !NEUMANN BC FOR PRESSURE AND DELTAP
              select case(ivar)
              case (PRES_VAR,DELP_VAR,TVIS_VAR)
              k = 2*guard+1 
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
              end do                          
              end select


              else ! BOUNDARY CONDITIONS ON VELOCITIES               
              if (isFace) then ! Set to zero velocities normal to boundary, up to boundary
                 if((ivar == VELC_FACE_VAR)) then 
                 k = 2*guard+2
                 do i = 1,guard
                    regionData(i,1:je,1:ke,ivar)= -regionData(k-i,1:je,1:ke,ivar)
                 end do
                 regionData(guard+1,1:je,1:ke,ivar)= 0.
                 endif
              else             ! Use guardcells to set to zero normal gradients of velocities not normal to boundary, at boundary
                 k = 2*guard+1   
                 if(ivar == VELC_FACE_VAR) then                               
                 do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= 2. - regionData(k-i,1:je,1:ke,ivar)
                 end do
                 endif
              endif

              endif


           case(INFLOW_INS)


              if (gridDataStruct==WORK) then ! NEUMANN 
              k = 2*guard+1 
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
              end do 

              elseif (gridDataStruct==CENTER) then !NEUMANN BC FOR PRESSURE AND DELTAP
              select case(ivar)
              case (PRES_VAR,DELP_VAR,TVIS_VAR)
              k = 2*guard+1 
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
              end do                          
              end select


              else ! BOUNDARY CONDITIONS ON VELOCITIES               

              if (isFace) then ! Set to zero velocities normal to boundary, up to boundary
                 if(ivar == VELC_FACE_VAR) then
                 do i = 1,guard+1
                    regionData(i,1:je,1:ke,ivar)= 1.
                 end do
                 endif
              else             ! Use guardcells to set to zero normal gradients of velocities not normal to boundary, at boundary
                 k = 2*guard+1   
                 if(ivar == VELC_FACE_VAR) then                               
                 do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)
                 end do
                 endif
              endif


!!$              ! Taylor Vortex velocity inflow + uconv 1. only on x Axis
!!$              if (axis .eq. IAXIS) then 
!!$              ! Get blocks dx, dy ,dz:
!!$              call Grid_getDeltas(blockHandle,del)
!!$              ! Get blocks coord and bsize
!!$              ! Bounding box:
!!$              call Grid_getBlkBoundBox(blockHandle,boundBox)
!!$              bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)
!!$              call Grid_getBlkCenterCoords(blockHandle,coord)
!!$                 
!!$              TLEVEL = ins_tlevel
!!$
!!$              write(*,*) 'ins_outflowgridChanged=',ins_outflowgridChanged
!!$
!!$              !if (.not. ins_outflowgridChanged) then
!!$              if (isFace) then ! Set to zero velocities normal to boundary, up to boundary
!!$                 if(ivar == VELC_FACE_VAR) then
!!$                 countj = 0
!!$                 do j = endpoints(LOW,JAXIS),endpoints(HIGH,JAXIS)
!!$                    countj = countj + 1
!!$                    ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
!!$                             real(j - NGUARD - 1)*del(JAXIS)  +  &
!!$                             0.5*del(JAXIS)
!!$                 do i = 1,guard+1
!!$                    xedge = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
!!$                            real(i - NGUARD - 1)*del(IAXIS)
!!$
!!$                    regionData(i,countj,1:ke,ivar)=-EXP(-2.0*ins_invRe*(TLEVEL)*1.)*        &
!!$                                                    COS(xedge-1.*(TLEVEL))*SIN(ycell) + 1.
!!$                 end do
!!$                 end do
!!$                 endif
!!$              else             ! Use guardcells to set to zero normal gradients of velocities not normal to boundary, at boundary
!!$                 k = 2*guard+1   
!!$                 if(ivar == VELC_FACE_VAR) then
!!$                 countj = 0                               
!!$                 do j = endpoints(LOW,JAXIS),endpoints(HIGH,JAXIS)
!!$                    countj = countj + 1
!!$                    yedge = coord(JAXIS) - bsize(JAXIS)/2.0 +   &
!!$                            real(j - NGUARD - 1)*del(JAXIS) 
!!$                 do i = 1,guard
!!$                    xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
!!$                            real(i - NGUARD - 1)*del(IAXIS) +   &
!!$                            0.5*del(IAXIS)
!!$
!!$                    regionData(i,countj,1:ke,ivar)= EXP(-2.0*ins_invRe*TLEVEL*1.)*        &
!!$                                                    SIN(xcell-1.*TLEVEL)*COS(yedge)
!!$                 end do
!!$                 end do
!!$                 endif
!!$              endif
!!$              !endif
!!$              endif
!!$
              endif


           case(OUTFLOW_INS) ! Convective outflow boundary condition. Only done for face == high  
           
              call Driver_abortFlash("OUTFLOW_INS boundary condition is not supported on Lower Face")

           case default
!!              print*,'boundary is',bcType
!!              call Driver_abortFlash("unsupported boundary condition on Lower Face")
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
           case(GRIDBC_EXTRAPOLATE_NSC)
              k=guard
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(k+i,1:je,1:ke,ivar)= (1+i)*regionData(k,1:je,1:ke,ivar) &
                      - i*regionData(k-1,1:je,1:ke,ivar)
                 where (regionData(k+i,1:je,1:ke,ivar)*regionData(k,1:je,1:ke,ivar) < 0.0)
                    regionData(k+i,1:je,1:ke,ivar)= 0.0
                 end where
#ifdef PRES_VAR
                 if (ivar==PRES_VAR) then
                    smallP = gr_smallRho * gr_smallE
                    regionData(k+i,1:je,1:ke,ivar)= max(smallP,regionData(k+i,1:je,1:ke,ivar))
                 end if
#endif
#ifdef ENER_VAR
                 if (ivar==ENER_VAR) then
                    regionData(k+i,1:je,1:ke,ivar)= max(gr_smallE,regionData(k+i,1:je,1:ke,ivar))
                 end if
#endif
#ifdef EINT_VAR
                 if (ivar==EINT_VAR) then
                    regionData(k+i,1:je,1:ke,ivar)= max(gr_smallE,regionData(k+i,1:je,1:ke,ivar))
                 end if
#endif
#ifdef DENS_VAR
                 if (ivar==DENS_VAR) then
                    regionData(k+i,1:je,1:ke,ivar)= max(gr_smallRho,regionData(k+i,1:je,1:ke,ivar))
                 end if
#endif
              end do

           case(GRIDBC_ZERO)
              k = 2*guard+1
              if(isFace)k=k+1
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= 0.0
              end do

           case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
                        HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL)
              k=guard
              if(isFace)k=k+1
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


           ! The following cases correspond to Incompressible solver BC combinations:
           ! Staggered grid, velocities: face centered, pressure: cell centered. 
           case(NEUMANN_INS)  ! Neumann BC

              if (gridDataStruct==WORK) then ! NEUMANN

              k = 2*guard+1
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do

              elseif (gridDataStruct==CENTER) then !NEUMANN BC FOR PRESSURE AND DELTAP
              select case(ivar)
              case (PRES_VAR,DELP_VAR,TVIS_VAR)
              k = 2*guard+1 
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do                          
              end select


              else ! BOUNDARY CONDITIONS ON VELOCITIES  

              k = 2*guard+1
              if(isFace) then
              ! First order down-wind for collocated var in the face:
              do i =1,guard
                 regionData(guard+1+i,1:je,1:ke,ivar)= regionData(guard+i,1:je,1:ke,ivar)
              enddo
              else
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do
              endif

              endif

           case(NOSLIP_INS)
 
              if (gridDataStruct==WORK) then ! NEUMANN 

              k = 2*guard+1 
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do                          

              elseif (gridDataStruct==CENTER) then !NEUMANN BC FOR PRESSURE AND DELTAP
              select case(ivar)
              case (PRES_VAR,DELP_VAR,TVIS_VAR)
              k = 2*guard+1 
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do                          
              end select


              else ! BOUNDARY CONDITIONS ON VELOCITIES               
              k = 2*guard+1
              if (isFace) then ! Set to zero velocities normal to boundary, up to boundary
                 if(ivar == VELC_FACE_VAR) then
                 k = k+1
                 do i = 1,guard
                    regionData(k-i,1:je,1:ke,ivar)= -regionData(i,1:je,1:ke,ivar)
                 end do
                 regionData(guard+1,1:je,1:ke,ivar)= 0.
                 endif
              else             ! Use guardcells to set to zero velocities not normal to boundary, at boundary
                 if(ivar == VELC_FACE_VAR) then                               
                 do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= -regionData(i,1:je,1:ke,ivar)
                 end do
                 endif
              endif

              endif           

           case(SLIP_INS)


              if (gridDataStruct==WORK) then ! NEUMANN 
              k = 2*guard+1 
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do                          

              elseif (gridDataStruct==CENTER) then !NEUMANN BC FOR PRESSURE AND DELTAP
              select case(ivar)
              case (PRES_VAR,DELP_VAR,TVIS_VAR)
           
              k = 2*guard+1 
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do                          
              end select


              else ! BOUNDARY CONDITIONS ON VELOCITIES               
              k = 2*guard+1
              if (isFace) then ! Set to zero velocities normal to boundary, up to boundary
                 if((ivar == VELC_FACE_VAR)) then 
                 k=k+1
                 do i = 1,guard
                    regionData(k-i,1:je,1:ke,ivar)= -regionData(i,1:je,1:ke,ivar)
                 end do
                 regiondata(guard+1,1:je,1:ke,ivar)= 0.
                 endif
              else             ! Use guardcells to set to zero velocities not normal to boundary, at boundary
                 if(ivar == VELC_FACE_VAR) then       
                 do i = 1,guard !-1
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
                 end do
                 !regionData(guard+1,1:je,1:ke,ivar)= (1.-1./16.)*regionData(guard,1:je,1:ke,ivar)+3./16.*regionData(guard-1,1:je,1:ke,ivar)-&
                 !                                  3./16.*regionData(guard-2,1:je,1:ke,ivar)+1./16.*regionData(guard-3,1:je,1:ke,ivar)
                 endif
              endif

              endif  


           case(MOVLID_INS)

              if (gridDataStruct==WORK) then ! NEUMANN 
              k = 2*guard+1 
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do                          

              elseif (gridDataStruct==CENTER) then !NEUMANN BC FOR PRESSURE AND DELTAP
              select case(ivar)
              case (PRES_VAR,DELP_VAR,TVIS_VAR)
           
              k = 2*guard+1 
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do                          
              end select


              else ! BOUNDARY CONDITIONS ON VELOCITIES               
              k = 2*guard+1
              if (isFace) then ! Set to zero velocities normal to boundary, up to boundary
                 if((ivar == VELC_FACE_VAR)) then 
                 k=k+1
                 do i = 1,guard+1
                    regionData(k-i,1:je,1:ke,ivar)= 0.
                 end do
                 endif
              else             ! Use guardcells to set to zero velocities not normal to boundary, at boundary
                 if(ivar == VELC_FACE_VAR) then       
                 do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= 2. - regionData(i,1:je,1:ke,ivar)
                 end do
                 endif
              endif

              endif  

           case(INFLOW_INS)

              if (gridDataStruct==WORK) then ! NEUMANN 
              k = 2*guard+1 
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do                          

              elseif (gridDataStruct==CENTER) then !NEUMANN BC FOR PRESSURE AND DELTAP
              select case(ivar)
              case (PRES_VAR,DELP_VAR,TVIS_VAR)
              k = 2*guard+1 
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do                          
              end select

              else ! BOUNDARY CONDITIONS ON VELOCITIES
              k = 2*guard+1               
              if (isFace) then ! Set to zero velocities normal to boundary, up to boundary
                 if(ivar == VELC_FACE_VAR) then
                 k=k+1
                 do i = 1,guard+1
                    regionData(k-i,1:je,1:ke,ivar)= -1.
                 end do
                 endif
              else             ! Use guardcells to set to zero normal gradients of velocities not normal to boundary, at boundary
                 if(ivar == VELC_FACE_VAR) then                               
                 do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
                 end do
                 endif
              endif

              endif




           case(OUTFLOW_INS) ! Convective outflow boundary condition.



              if (gridDataStruct==WORK) then ! NEUMANN 
              k = 2*guard+1 
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do                          

              elseif (gridDataStruct==CENTER) then !NEUMANN BC FOR PRESSURE AND DELTAP
              select case(ivar)
              case (PRES_VAR,DELP_VAR,TVIS_VAR)
              k = 2*guard+1 
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)
              end do                          
              end select

              else ! BOUNDARY CONDITIONS ON VELOCITIES - ONLY for 2nd ORDER STAGGERED GRIDS !

              alfadt = ins_alfa*dr_dt        
       
              ! Get blocks dx, dy ,dz:
              call Grid_getDeltas(blockHandle,del)              

              ! X direction:
              if (axis .eq. IAXIS) then               

                 !write(*,*) 'In Outflow',ins_predcorrflg

                 ja = endPoints(LOW,JAXIS)
                 jb = endPoints(HIGH,JAXIS)
                 ka = endPoints(LOW,KAXIS)
                 kb = endPoints(HIGH,KAXIS)

              if (ins_predcorrflg) then
              if (isFace) then 
                 if(ivar == VELC_FACE_VAR) then ! U velocities on X face grid
             
                   regionData(guard+1,1:je,1:ke,ivar) = uvel_x(guard+1,ja:jb,ka:kb,HIGH,blockHandle) -&
                    (ins_convvel(HIGH,axis)*alfadt/del(axis))* &
                    (uvel_x(guard+1,ja:jb,ka:kb,HIGH,blockHandle) - &
                     uvel_x(guard,ja:jb,ka:kb,HIGH,blockHandle))

                 endif

              elseif(gridDataStruct .eq. FACEY) then ! V velocities on X direction             
                 if(ivar == VELC_FACE_VAR) then                               

                   vvel_x(guard+1,ja:jb,ka:kb,HIGH,blockHandle) = vvel_x(guard,ja:jb,ka:kb,HIGH,blockHandle) - & 
                        (ins_convvel(HIGH,axis)*alfadt/del(axis))* &
                        (vvel_x(guard,ja:jb,ka:kb,HIGH,blockHandle) - &
                         vvel_x(guard-1,ja:jb,ka:kb,HIGH,blockHandle))

                   regionData(guard+1,1:je,1:ke,ivar) =  vvel_x(guard+1,ja:jb,ka:kb,HIGH,blockHandle)

                 endif
#if NDIM == 3 
              elseif(gridDataStruct .eq. FACEZ) then ! W velocities on X direction 
                 if(ivar == VELC_FACE_VAR) then

                   wvel_x(guard+1,ja:jb,ka:kb,HIGH,blockHandle) = wvel_x(guard,ja:jb,ka:kb,HIGH,blockHandle) - & 
                        (ins_convvel(HIGH,axis)*alfadt/del(axis))* &
                        (wvel_x(guard,ja:jb,ka:kb,HIGH,blockHandle) - &
                         wvel_x(guard-1,ja:jb,ka:kb,HIGH,blockHandle))

                   regionData(guard+1,1:je,1:ke,ivar) =  wvel_x(guard+1,ja:jb,ka:kb,HIGH,blockHandle)


                 endif
#endif
              endif

              else !ins_predcorrflg

              if (isFace) then 
!!$                 if(ivar == VELC_FACE_VAR) then ! U velocities on X face grid
!!$                 endif
              elseif(ins_outflowgridChanged) then ! Estrapolate linearly V or W from the interior when grid changes.
                 regionData(guard+1,1:je,1:ke,ivar) = 2.*regionData(guard,1:je,1:ke,ivar) - &
                                                         regionData(guard-1,1:je,1:ke,ivar)

              elseif(gridDataStruct .eq. FACEY) then ! V velocities on X direction             
                 if(ivar == VELC_FACE_VAR) then                               

                    regionData(guard+1,1:je,1:ke,ivar) =  vvel_x(guard+1,ja:jb,ka:kb,HIGH,blockHandle)

                 endif 

#if NDIM == 3
              elseif(gridDataStruct .eq. FACEZ) then ! W velocities on X direction 
                 if(ivar == VELC_FACE_VAR) then

                    regionData(guard+1,1:je,1:ke,ivar) =  wvel_x(guard+1,ja:jb,ka:kb,HIGH,blockHandle)

                 endif
#endif
              endif


              endif !ins_predcorrflg
              endif !IAXIS



              ! Y direction:
              if (axis .eq. JAXIS) then               

                 ia = endPoints(LOW,IAXIS)
                 ib = endPoints(HIGH,IAXIS)
                 ka = endPoints(LOW,KAXIS)
                 kb = endPoints(HIGH,KAXIS)

              if (ins_predcorrflg) then
              if (isFace) then 
                 if(ivar == VELC_FACE_VAR) then ! V velocities on Y face grid
             
                   regionData(guard+1,1:je,1:ke,ivar) = vvel_y(guard+1,ia:ib,ka:kb,HIGH,blockHandle) -&
                    (ins_convvel(HIGH,axis)*alfadt/del(axis))* &
                    (vvel_y(guard+1,ia:ib,ka:kb,HIGH,blockHandle) - &
                     vvel_y(guard,ia:ib,ka:kb,HIGH,blockHandle))

                 endif

              elseif(gridDataStruct .eq. FACEX) then ! U velocities on Y direction             
                 if(ivar == VELC_FACE_VAR) then                               

                   uvel_y(guard+1,ia:ib,ka:kb,HIGH,blockHandle) = uvel_y(guard,ia:ib,ka:kb,HIGH,blockHandle) - & 
                        (ins_convvel(HIGH,axis)*alfadt/del(axis))* &
                        (uvel_y(guard,ia:ib,ka:kb,HIGH,blockHandle) - &
                         uvel_y(guard-1,ia:ib,ka:kb,HIGH,blockHandle))

                   regionData(guard+1,1:je,1:ke,ivar) =  uvel_y(guard+1,ia:ib,ka:kb,HIGH,blockHandle)

                 endif

#if NDIM == 3 
              elseif(gridDataStruct .eq. FACEZ) then ! W velocities on Y direction 
                 if(ivar == VELC_FACE_VAR) then

                   wvel_y(guard+1,ia:ib,ka:kb,HIGH,blockHandle) = wvel_y(guard,ia:ib,ka:kb,HIGH,blockHandle) - & 
                        (ins_convvel(HIGH,axis)*alfadt/del(axis))* &
                        (wvel_y(guard,ia:ib,ka:kb,HIGH,blockHandle) - &
                         wvel_y(guard-1,ia:ib,ka:kb,HIGH,blockHandle))

                   regionData(guard+1,1:je,1:ke,ivar) =  wvel_y(guard+1,ia:ib,ka:kb,HIGH,blockHandle)


                 endif
#endif

              endif

              else !ins_predcorrflg

              if (isFace) then 
!!$                 if(ivar == VELC_FACE_VAR) then ! V velocities on Y face grid
!!$                 endif
              elseif(ins_outflowgridChanged) then ! Estrapolate linearly U or W from the interior when grid changes.
                 regionData(guard+1,1:je,1:ke,ivar) = 2.*regionData(guard,1:je,1:ke,ivar) - &
                                                         regionData(guard-1,1:je,1:ke,ivar)

              elseif(gridDataStruct .eq. FACEX) then ! U velocities on Y direction             
                 if(ivar == VELC_FACE_VAR) then                               

                    regionData(guard+1,1:je,1:ke,ivar) =  uvel_y(guard+1,ia:ib,ka:kb,HIGH,blockHandle)

                 endif 

#if NDIM == 3
              elseif(gridDataStruct .eq. FACEZ) then ! W velocities on Y direction 
                 if(ivar == VELC_FACE_VAR) then

                    regionData(guard+1,1:je,1:ke,ivar) =  wvel_y(guard+1,ia:ib,ka:kb,HIGH,blockHandle)

                 endif
#endif
              endif


              endif !ins_predcorrflg
              endif !JAXIS

#if NDIM == 3

              ! Z direction:
              if (axis .eq. KAXIS) then               

                 ia = endPoints(LOW,IAXIS)
                 ib = endPoints(HIGH,IAXIS)
                 ja = endPoints(LOW,JAXIS)
                 jb = endPoints(HIGH,JAXIS)

              if (ins_predcorrflg) then
              if (isFace) then 
                 if(ivar == VELC_FACE_VAR) then ! W velocities on Z face grid
             
                   regionData(guard+1,1:je,1:ke,ivar) = wvel_z(guard+1,ia:ib,ja:jb,HIGH,blockHandle) -&
                    (ins_convvel(HIGH,axis)*alfadt/del(axis))* &
                    (wvel_z(guard+1,ia:ib,ja:jb,HIGH,blockHandle) - &
                     wvel_z(guard,ia:ib,ja:jb,HIGH,blockHandle))

                 endif


              elseif(gridDataStruct .eq. FACEX) then ! U velocities on Z direction             
                 if(ivar == VELC_FACE_VAR) then                               

                   uvel_z(guard+1,ia:ib,ja:jb,HIGH,blockHandle) = uvel_z(guard,ia:ib,ja:jb,HIGH,blockHandle) - & 
                        (ins_convvel(HIGH,axis)*alfadt/del(axis))* &
                        (uvel_z(guard,ia:ib,ja:jb,HIGH,blockHandle) - &
                         uvel_z(guard-1,ia:ib,ja:jb,HIGH,blockHandle))

                   regionData(guard+1,1:je,1:ke,ivar) =  uvel_z(guard+1,ia:ib,ja:jb,HIGH,blockHandle)

                 endif


              elseif(gridDataStruct .eq. FACEY) then ! V velocities on Z direction 
                 if(ivar == VELC_FACE_VAR) then

                   vvel_z(guard+1,ia:ib,ja:jb,HIGH,blockHandle) = vvel_z(guard,ia:ib,ja:jb,HIGH,blockHandle) - & 
                        (ins_convvel(HIGH,axis)*alfadt/del(axis))* &
                        (vvel_z(guard,ia:ib,ja:jb,HIGH,blockHandle) - &
                         vvel_z(guard-1,ia:ib,ja:jb,HIGH,blockHandle))

                   regionData(guard+1,1:je,1:ke,ivar) =  vvel_z(guard+1,ia:ib,ja:jb,HIGH,blockHandle)


                 endif

              endif

              else !ins_predcorrflg

              if (isFace) then 
!!$                 if(ivar == VELC_FACE_VAR) then ! W velocities on Z face grid
!!$                 endif
              elseif(ins_outflowgridChanged) then ! Estrapolate linearly V or U from the interior when grid changes.
                 regionData(guard+1,1:je,1:ke,ivar) = 2.*regionData(guard,1:je,1:ke,ivar) - &
                                                         regionData(guard-1,1:je,1:ke,ivar)

              elseif(gridDataStruct .eq. FACEX) then ! U velocities on Z direction             
                 if(ivar == VELC_FACE_VAR) then                               

                    regionData(guard+1,1:je,1:ke,ivar) =  uvel_z(guard+1,ia:ib,ja:jb,HIGH,blockHandle)

                 endif 

              elseif(gridDataStruct .eq. FACEY) then ! V velocities on Z direction 
                 if(ivar == VELC_FACE_VAR) then

                    regionData(guard+1,1:je,1:ke,ivar) =  vvel_z(guard+1,ia:ib,ja:jb,HIGH,blockHandle)

                 endif

              endif


              endif !ins_predcorrflg
              endif !KAXIS



#endif



              endif ! Velocities

           case default
!!              print*,'boundary is',bcType
!!              call Driver_abortFlash("unsupported boundary condition on Upper Face")
           end select


        end if
     end if
  end do

  return
end subroutine Grid_bcApplyToRegion
