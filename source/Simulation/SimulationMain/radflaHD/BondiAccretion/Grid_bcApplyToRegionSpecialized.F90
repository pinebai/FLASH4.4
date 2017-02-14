!!****if* source/Simulation/SimulationMain/radflaHD/BondiAccretion/Grid_bcApplyToRegionSpecialized
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
!!   Where masked(variables)
!!     If (face=LOW)  
!!       regionData(1:guard,:,:,variables) =  boundary values
!!     If (face=HIGH) 
!!       regionData(regionSize(BC_DIR)-guard+1:regionSize(BC_DIR),:,:,variables) =  boundary values
!!
!!
!! ARGUMENTS 
!!
!! 1. BASIC ARGUMENTS
!!
!!    bcType - the type of boundary condition being applied.
!!    gridDataStruct - the Grid dataStructure, should be given as
!!                     one of the constants CENTER, FACEX, FACEY, FACEZ.
!!    guard -    number of guard cells
!!    axis  - the direction along which to apply boundary conditions,
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
!!    regionSize     : regionSize(BC_DIR) contains the size of each row of data
!!                     in the regionData array.  With row we mean here an array slice
!!                     regionData(:,I2,I3,VAR), corresponding to cells that are situated
!!                     along a line in the 'axis' direction. For the common case of guard=4,
!!                     (e.g., when gridDataStruct=CENTER) and either 8 or 9 for face-
!!                     centered data, depending on the direction given by axis.
!!                     regionSize(SECOND_DIR) contains the number of rows along the
!!                     second direction, and regionSize(THIRD_DIR) has the number of rows
!!                     along the third direction. (See also below under secondDir,thirdDir
!!                     for the meaning of second and third direction; and see also NOTE (1)
!!                     below.)
!!                     Finally, regionSize(GRID_DATASTRUCT) contains the
!!                     number of variables in the data structure.
!!  mask - if present, boundary conditions are to be applied only to selected variables.
!!         However, an implementation of this interface may ignore the mask argument;
!!         a mask should be understood as a possible opportunity for optimization which
!!         an implementation may ignore.
!!         Specifying a mask does not mean that previous values of other variables in
!!         guard cells will be left undisturbed.
!!    applied - is set true if this routine has handled the given bcType, otherwise it is 
!!              set to false.
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
!!                         the sweep direction.  SecondDir and thirdDir give
!!                         the meaning of the second and third dimension,
!!                         respectively, of the regionData array.
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
!!          In other words, an implementation can nearly always ignore this optional
!!          argument.  As of FLASH 3.0, it is only used internally within the
!!          Grid unit and is handled by the GridBoundaryConditions/Grid_bcApplyToRegion
!!          implementation. It is used within the Grid unit by a Multigrid GridSolver
!!          implementation which requires some special handling, but this is only
!!          applied to the WORK data structure.  The argument has been added to the
!!          Grid_bcApplyToRegionSpecialized interface for consistency with
!!          Grid_bcApplyToRegion.
!!
!! NOTES
!!
!! (1)        NOTE that the second indices of the endPoints and
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
!! SEE ALSO
!!
!!   Grid_bcApplyToRegion            
!!
!!***


subroutine Grid_bcApplyToRegionSpecialized(bcType,gridDataStruct,&
     guard,axis,face,regionData,regionSize,mask,applied,&
     blockHandle,secondDir,thirdDir,endPoints,blkLimitsGC, idest)

#include "constants.h"
#include "Flash.h"

  use Eos_interface, ONLY : Eos_arrayWrapped
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_bcInterface, ONLY : gr_bcMapBcType
  use Simulation_data, ONLY : sim_sb, sim_speedlt, sim_useMGD
  use Grid_data, ONLY : gr_meshMe, gr_geometry, gr_dirGeom, &
       gr_smallrho, gr_smallE
  use Simulation_data, ONLY : sim_smallP, sim_smallEele, &
                              sim_accretionRate, sim_rho_vac, sim_t_vac
  use sim_interface, ONLY : sim_computeAnaBondi

  implicit none

  integer, intent(IN) :: bcType,axis,face,guard,gridDataStruct
  integer,dimension(REGION_DIM),intent(IN) :: regionSize
  real,dimension(regionSize(BC_DIR),&
       regionSize(SECOND_DIR),&
       regionSize(THIRD_DIR),&
       regionSize(STRUCTSIZE)),intent(INOUT)::regionData
  target :: regionData
  logical,intent(IN),dimension(regionSize(STRUCTSIZE)):: mask
  logical, intent(OUT) :: applied
  integer,intent(IN) :: blockHandle
  integer,intent(IN) :: secondDir,thirdDir
  integer,intent(IN),dimension(LOW:HIGH,MDIM) :: endPoints, blkLimitsGC
  integer,intent(IN),OPTIONAL:: idest

  integer :: i,j, k,ivar,je,ke,n,varCount,bcTypeActual
  logical :: isFace
  integer :: sign
  real    :: smallP
  integer :: sizeGC
  real, allocatable, dimension(:) :: cellCenterCoord
  real    :: delR,drh, rim,rip,rgm,rgp 
  real    :: unused_vel
  real, parameter :: alphaDens = 1.0 ! 1.0 better than -1.5 based on resulting Mdot ?
  real,dimension(1:regionSize(SECOND_DIR),1:regionSize(THIRD_DIR)) :: rhoLeft
  real,dimension(1:regionSize(SECOND_DIR),1:regionSize(THIRD_DIR)) :: denom,a,rho0
  integer,dimension(LOW:HIGH,MDIM) :: eosRange
  real,pointer :: dataForEosPtr(:,:,:,:), data1DForEosPtr(:,:,:,:)
  real,target  :: dataForEos(regionSize(STRUCTSIZE),1,1,1)

  select case (bcType)
  case(USER_DEFINED)
     if (gridDataStruct .NE. CENTER) then
        applied = .FALSE.       !should be picked up by Grid_bcApplyToRegion implementation
        return                  !RETURN immediately!
     else
        applied = .TRUE.           !will handle this type below.
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

  
!!  print*,'in applyBcRegion ',varCount,gridDataStruct,WORK,guard,axis,face

  ! get coordinates
  sizeGC = blkLimitsGC(HIGH,axis)
  allocate(cellCenterCoord(sizeGC))
  call gr_extendedGetCellCoords(axis, blockHandle, gr_meshMe, CENTER, .true., cellCenterCoord, sizeGC)

  delR = cellCenterCoord(guard+1) - cellCenterCoord(guard)
  drh = 0.5*delR

  if(mask(DENS_VAR)) then
     if(face==LOW) then
        denom = gF(guard+2) - gF(guard+1)
        a = (regionData(guard+2,1:je,1:ke,DENS_VAR) - regionData(guard+1,1:je,1:ke,DENS_VAR)) / denom
        rho0 = (  gF(guard+2)*regionData(guard+1,1:je,1:ke,DENS_VAR)   &
                - gF(guard+1)*regionData(guard+2,1:je,1:ke,DENS_VAR) ) &
                / denom
        do i = 1,guard
           regionData(guard+1-i,1:je,1:ke,DENS_VAR)= rho0 + a*gF(guard+1-i)
           regionData(guard+1-i,1:je,1:ke,DENS_VAR)= max(gr_smallRho,regionData(guard+1-i,1:je,1:ke,DENS_VAR))
           forall (n=1:ke, j=1:je, (regionData(guard+1-i,j,n,DENS_VAR)*regionData(guard+1,j,n,DENS_VAR) < 0.0))
              regionData(guard+1-i,j,n,DENS_VAR) = 0.0
           end forall
        end do
     else  !(face==HIGH)
        k=guard
        if(isFace)k=k+1
#if(0)
        denom = gF(endPoints(LOW,axis)-1+k) - gF(endPoints(LOW,axis)-2+k)
        !!print*,'densDeno:', denom
        if (ANY(denom==0.0)) print*,'***WARNING*** denom should not be zero!'
        a = (regionData(k,1:je,1:ke,DENS_VAR) - regionData(k-1,1:je,1:ke,DENS_VAR)) / denom
        rho0 = (  gF(endPoints(LOW,axis)-1+k)*regionData(k-1,1:je,1:ke,DENS_VAR)   &
             - gF(endPoints(LOW,axis)-2+k)*regionData(k  ,1:je,1:ke,DENS_VAR) ) &
             / denom
        !!print*,'densrho0:', rho0
        !!print*,'dens a :', a
        do i = 1,guard
           regionData(k+i,1:je,1:ke,DENS_VAR)= rho0 + a*gF(endPoints(LOW,axis)-1+k+i)
           regionData(k+i,1:je,1:ke,DENS_VAR)= max(gr_smallRho,regionData(k+i,1:je,1:ke,DENS_VAR))
           !regionData(k+i,1:je,1:ke,DENS_VAR)= max(sim_rho_vac,regionData(k+i,1:je,1:ke,DENS_VAR))
        end do
#else
        do i = 1,guard
           do n=1,ke
              do j=1,je
                 call sim_computeAnaBondi(cellCenterCoord(endPoints(LOW,axis)-1+k+i),&
                      regionData(k+i,j,n,VELX_VAR), &
                      regionData(k+i,j,n,DENS_VAR) )
              end do
           end do
        end do
#endif
     end if
  end if

  do ivar = 1,varCount
     if (ivar==DENS_VAR) CYCLE  !handled above
     if(mask(ivar)) then
        call gr_bcMapBcType(bcTypeActual,bcType,ivar,gridDataStruct,axis,face,idest)

        if(face==LOW) then
           if (ivar==VELX_VAR) then
              bcTypeActual = 1001
           else if (ivar==DENS_VAR) then
              ! leave as USER_DEFINED
           else
#if defined(ERAD_VAR) && defined(TRAD_VAR)
              if (ivar==ERAD_VAR .OR. ivar==TRAD_VAR .OR. (NMASS_SCALARS==1 .AND. ivar==MASS_SCALARS_BEGIN)) then
              ! leave as USER_DEFINED
              elseif (ivar==FLLM_VAR .OR. ivar==FLXL_VAR) then
              ! leave as USER_DEFINED
              else
#endif
                 bcTypeActual = GRIDBC_EXTRAPOLATE_NSC
#if defined(ERAD_VAR) && defined(TRAD_VAR)
              end if
#endif
           end if

        else if(face==HIGH) then
           if (ivar==VELX_VAR) then
              bcTypeActual = 1000
           else if (ivar==DENS_VAR) then
              ! leave as USER_DEFINED
           else
#if defined(ERAD_VAR) && defined(TRAD_VAR)
              if (ivar==ERAD_VAR .OR. ivar==TRAD_VAR .OR. (NMASS_SCALARS==1 .AND. ivar==MASS_SCALARS_BEGIN)) then
                 if (.NOT. sim_useMGD) bcTypeActual = GRIDBC_ZERO
                 ! else leave as USER_DEFINED
              elseif (ivar==FLLM_VAR .OR. ivar==FLXL_VAR) then
                 ! leave as USER_DEFINED
              else
#endif
                 !bcTypeActual = GRIDBC_EXTRAPOLATE_NSC
                 ! leave as USER_DEFINED !!!
#if defined(ERAD_VAR) && defined(TRAD_VAR)
              end if
#endif
           endif
        end if

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
                 if (ivar==PRES_VAR) then
                    smallP = gr_smallRho * gr_smallE
                    regionData(guard+1-i,1:je,1:ke,ivar)= max(smallP,regionData(guard+1-i,1:je,1:ke,ivar))
                 else if (ivar==EINT_VAR .OR. ivar==ENER_VAR) then
                    regionData(guard+1-i,1:je,1:ke,ivar)= max(gr_smallE,regionData(guard+1-i,1:je,1:ke,ivar))
#ifdef EELE_VAR
                 else if (ivar==EELE_VAR) then
                    regionData(guard+1-i,1:je,1:ke,ivar)= max(sim_smallEele,regionData(guard+1-i,1:je,1:ke,ivar))
#endif
#ifdef PELE_VAR
                 else if (ivar==PELE_VAR) then
                    smallP = gr_smallRho * gr_smallE
                    regionData(guard+1-i,1:je,1:ke,ivar)= max(smallP,regionData(guard+1-i,1:je,1:ke,ivar))
#endif
                 else if (ivar==DENS_VAR) then
                    regionData(guard+1-i,1:je,1:ke,ivar)= max(gr_smallRho,regionData(guard+1-i,1:je,1:ke,ivar))
                 end if
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

           case(1001)
              rim = cellCenterCoord(guard+1) - drh 
              rip = cellCenterCoord(guard+1) + drh 
              do i = 1,guard
                 regionData(guard+1-i,1:je,1:ke,ivar)= &
                      regionData(guard+1,1:je,1:ke,ivar) * regionData(guard+1,1:je,1:ke,DENS_VAR)
                 rgm = cellCenterCoord(guard+1-i) - drh 
                 rgp = cellCenterCoord(guard+1-i) + drh 
                 regionData(guard+1-i,1:je,1:ke,ivar)= &
                      regionData(guard+1-i,1:je,1:ke,ivar) &
                      *(rim**2 + rim*rip + rip**2) / (rgm**2 + rgm*rgp + rgp**2)
                 regionData(guard+1-i,1:je,1:ke,ivar)= &
                      regionData(guard+1-i,1:je,1:ke,ivar) / regionData(guard+1-i,1:je,1:ke,DENS_VAR)
                 regionData(guard+1-i,1:je,1:ke,ivar)= &
                      min(0.0,regionData(guard+1-i,1:je,1:ke,ivar))
              end do
           case(1002)
              do i = 1,guard
                 rgp = cellCenterCoord(guard+1-i) + drh 
                 regionData(guard+1-i,1:je,1:ke,ivar)= sim_accretionRate &
                      / (2*PI * &
                         rgp**2 * (rho0 + a * rgp**alphaDens) )
                 regionData(guard+1-i,1:je,1:ke,ivar)= &
                      regionData(guard+1-i,1:je,1:ke,ivar) - regionData(guard+2-i,1:je,1:ke,ivar)
              end do
           case(USER_DEFINED)
              if (ivar == DENS_VAR) then
                 denom = gF(guard+2) - gF(guard+1)
                 a = (regionData(guard+2,1:je,1:ke,DENS_VAR) - regionData(guard+1,1:je,1:ke,DENS_VAR)) / denom
                 rho0 = (  gF(guard+2)*regionData(guard+1,1:je,1:ke,DENS_VAR)   &
                         - gF(guard+1)*regionData(guard+2,1:je,1:ke,DENS_VAR) ) &
                         / denom
              end if
              do i = 1,guard
                 if (ivar == DENS_VAR) then
                    regionData(guard+1-i,1:je,1:ke,DENS_VAR)= rho0 + a*gF(guard+1-i)
                    regionData(guard+1-i,1:je,1:ke,ivar)= max(gr_smallRho,regionData(guard+1-i,1:je,1:ke,ivar))
                 else
                    regionData(guard+1-i,1:je,1:ke,ivar)= (1+i)*regionData(guard+1,1:je,1:ke,ivar) &
                         - i*regionData(guard+2,1:je,1:ke,ivar) !like GRIDBC_MG_EXTRAPOLATE
                 end if
                 forall (n=1:ke, j=1:je, (regionData(guard+1-i,j,n,ivar)*regionData(guard+1,j,n,ivar) < 0.0))
                    regionData(guard+1-i,j,n,ivar) = 0.0
                 end forall
!!$                 where (regionData(guard+1-i,1:je,1:ke,ivar)*regionData(guard+1,1:je,1:ke,ivar) < 0.0)
!!$                    regionData(guard+1-i,1:je,1:ke,ivar) = 0.0
!!$                    if (ivar==PRES_VAR) regionData(guard+1-i,1:je,1:ke,ivar) = sim_smallP
!!$                 end where
#if(0)
                 regionData(i,1:je,1:ke,ivar)= regionData(guard+1,1:je,1:ke,ivar)
#endif
                 if (sign == -1) then
                    regionData(i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar) * &
                         (cellCenterCoord(guard+1) / cellCenterCoord(i))**2
                    do n=1,ke
                       do j=1,je
                          regionData(i,j,n,ivar) = max(-3.0e10,min(regionData(i,j,n,ivar),0.0))
                       end do
                    end do
                 else if (ivar==TRAD_VAR) then
                    regionData(i,1:je,1:ke,ivar)= &
                         (1.6e+05*3.99e+33/(16.*PI*sim_sb*(cellCenterCoord(i))**2))**0.25
                 else if (ivar==ERAD_VAR .OR. ivar==MASS_SCALARS_BEGIN) then
                    regionData(i,1:je,1:ke,ivar)=  &
                         1.6e+05*3.99e+33 &
                         /(4*PI*sim_speedlt*(cellCenterCoord(i))**2  &
                           * regionData(i,1:je,1:ke,DENS_VAR))
                 else if (ivar==FLLM_VAR .OR. ivar==FLXL_VAR) then

                    regionData(i,1:je,1:ke,ivar)=  &
                         0.5 * cellCenterCoord(i)  &
                           * 0.4 * regionData(i,1:je,1:ke,DENS_VAR)    * 3.0

                 end if
                 if (ivar==PRES_VAR) then
                    forall (n=1:ke, j=1:je, (regionData(guard+1-i,j,n,ivar) == 0.0))
                       regionData(guard+1-i,j,n,ivar) = sim_smallP
                    end forall
                 end if
              end do
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
                 if (ivar==PRES_VAR) then
                    smallP = gr_smallRho * gr_smallE
                    regionData(k+i,1:je,1:ke,ivar)= max(smallP,regionData(k+i,1:je,1:ke,ivar))
                 else if (ivar==EINT_VAR .OR. ivar==ENER_VAR) then
                    regionData(k+i,1:je,1:ke,ivar)= max(gr_smallE,regionData(k+i,1:je,1:ke,ivar))
#ifdef EELE_VAR
                 else if (ivar==EELE_VAR) then
                    regionData(k+i,1:je,1:ke,ivar)= max(sim_smallEele,regionData(k+i,1:je,1:ke,ivar))
#endif
#ifdef PELE_VAR
                 else if (ivar==PELE_VAR) then
                    smallP = gr_smallRho * gr_smallE
                    regionData(k+i,1:je,1:ke,ivar)= max(smallP,regionData(k+i,1:je,1:ke,ivar))
#endif
                 else if (ivar==DENS_VAR) then
                    regionData(k+i,1:je,1:ke,ivar)= max(gr_smallRho,regionData(k+i,1:je,1:ke,ivar))
                 end if
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

           case(1000)
              ! ignore coz already done!
           case(1001)
              k=guard
              if(isFace)k=k+1
              rim = cellCenterCoord(endPoints(LOW,axis)-1+k) - drh 
              rip = cellCenterCoord(endPoints(LOW,axis)-1+k) + drh 
              do i = 1,guard
                 regionData(k+i,1:je,1:ke,ivar)= &
                      regionData(k,1:je,1:ke,ivar) * regionData(k,1:je,1:ke,DENS_VAR)
                 rgm = cellCenterCoord(endPoints(LOW,axis)-1+k+i) - drh
                 rgp = cellCenterCoord(endPoints(LOW,axis)-1+k+i) + drh
                 regionData(k+i,1:je,1:ke,ivar)= &
                      regionData(k+i,1:je,1:ke,ivar) &
                      *(rim**2 + rim*rip + rip**2) / (rgm**2 + rgm*rgp + rgp**2)
                 regionData(k+i,1:je,1:ke,ivar)= &
                      regionData(k+i,1:je,1:ke,ivar) / regionData(k+i,1:je,1:ke,DENS_VAR)
                 regionData(k+i,1:je,1:ke,ivar)= &
                      min(0.0,regionData(k+i,1:je,1:ke,ivar))
              end do
!!$              print*,'DENS_VAR', regionData(1:k+guard,1:je,1:ke,DENS_VAR)
!!$              print*,'velx R :', regionData(1:k+guard,1:je,1:ke,ivar)
           case(1002)
              k=guard
              if(isFace)k=k+1
              do i = 1,1
                 rgm = cellCenterCoord(endPoints(LOW,axis)-1+k+i) - drh
                 rhoLeft = max(gr_smallRho,rho0 + a * rgm**alphaDens)
!                 rhoLeft = max(sim_rho_vac,rho0 + a * rgm**alphaDens)
                 regionData(k+i,1:je,1:ke,ivar)= sim_accretionRate &
                      / (2*PI * &
                         rgm**2 * rhoLeft )
                 regionData(k+i,1:je,1:ke,ivar)= &
                      regionData(k+i,1:je,1:ke,ivar) - regionData(k+i-1,1:je,1:ke,ivar)
              end do
              do i = 2,guard
                 regionData(k+i,1:je,1:ke,ivar)= regionData(k+i-1,1:je,1:ke,ivar)
              end do
             !! print*,'vx rho0:', rho0
             !! print*,'velx a :', a
             !! print*,'DENS_VAR', regionData(1:k+guard,1:je,1:ke,DENS_VAR)
             !! print*,'rgm,rho(rgm)=',rgm,(rho0 + a * rgm**alphaDens)
             !! print*,'velx R :', regionData(1:k+guard,1:je,1:ke,ivar)
           case(USER_DEFINED)
              k=guard
              if(isFace)k=k+1
              if (ivar == DENS_VAR) then
                 denom = gF(endPoints(LOW,axis)-1+k) - gF(endPoints(LOW,axis)-2+k)
                 !!print*,'densDeno:', denom
                 if (ANY(denom==0.0)) print*,'***WARNING*** denom should not be zero!'
                 a = (regionData(k,1:je,1:ke,DENS_VAR) - regionData(k-1,1:je,1:ke,DENS_VAR)) / denom
                 rho0 = (  gF(endPoints(LOW,axis)-1+k)*regionData(k-1,1:je,1:ke,DENS_VAR)   &
                         - gF(endPoints(LOW,axis)-2+k)*regionData(k  ,1:je,1:ke,DENS_VAR) ) &
                         / denom
                 !!print*,'densrho0:', rho0
                 !!print*,'dens a :', a
              end if
              do i = 1,guard
                 if (ivar == DENS_VAR) then
                    regionData(k+i,1:je,1:ke,DENS_VAR)= rho0 + a*gF(endPoints(LOW,axis)-1+k+i)
                    regionData(k+i,1:je,1:ke,ivar)= max(gr_smallRho,regionData(k+i,1:je,1:ke,ivar))
!                    regionData(k+i,1:je,1:ke,ivar)= max(sim_rho_vac,regionData(k+i,1:je,1:ke,ivar))
                 else
                    regionData(k+i,1:je,1:ke,ivar)= regionData(k,1:je,1:ke,ivar)
                 end if
                 if (sign == -1) then
                    do n = 1,ke
                       do j = 1,je
                          regionData(k+i,j,n,ivar) = min(regionData(k+i,j,n,ivar),0.0)
                       end do
                    end do
#if defined(TION_VAR)
                 else if (ivar==TION_VAR) then
                    regionData(k+i,1:je,1:ke,ivar)= 0.0
#endif
#if defined(TELE_VAR)
                 else if (ivar==TELE_VAR) then
                    regionData(k+i,1:je,1:ke,ivar)= sim_t_vac
#endif
#if defined(ERAD_VAR) && defined(TRAD_VAR)
                 else if (ivar==TRAD_VAR) then
                    regionData(k+i,1:je,1:ke,ivar)= (1.6e+05*3.99e+33 &
                         / (16.*PI*sim_sb*&
                         (cellCenterCoord(endPoints(LOW,axis)-1+k+i))**2))**0.25
                 else if (ivar==ERAD_VAR .OR. ivar==MASS_SCALARS_BEGIN) then
                    regionData(k+i,1:je,1:ke,ivar)=  &
                         1.6e+05*3.99e+33 &
                         /(4.*PI*sim_speedlt*(cellCenterCoord(endPoints(LOW,axis)-1+k+i))**2  &
                           * regionData(k+i,1:je,1:ke,DENS_VAR))
                 else if (ivar==FLLM_VAR .OR. ivar==FLXL_VAR) then

                    regionData(k+i,1:je,1:ke,ivar)=  &
                         0.5 * cellCenterCoord(endPoints(LOW,axis)-1+k+i)  &
                           * 0.4 * regionData(k+i,1:je,1:ke,DENS_VAR)    * 3.0

#endif
                 end if
              end do

           case default
!!              print*,'boundary is',bcType
!!              call Driver_abortFlash("unsupported boundary condition on Upper Face")
           end select
        end if
     end if
  end do

  if(face==LOW) then
     ! for now, not callin Eos here
  else
     eosRange(LOW,:) = 1
     eosRange(HIGH,IAXIS)  = 1
     eosRange(HIGH,2:MDIM) = 1 ! regionSize(2:MDIM)
     do n = 1,ke
        do j = 1,je
           k=guard
           if(isFace)k=k+1
!!$           print*,'bc guard PRES_VAR before:',
           do i = 1,guard
              dataForEosPtr => regionData(k+i:k+i,j:j,n:n,:)
              dataForEos(:,1,1,1) = dataForEosPtr(1,1,1,:)
              data1DForEosPtr => dataForEos(:,:,:,:)
              call Eos_arrayWrapped(MODE_DENS_TEMP_GATHER            ,eosRange,data1DForEosPtr,CENTER)
              if (sim_useMGD) then
                 call Eos_arrayWrapped(MODE_DENS_EI_MAT_GATHER_PRADSCALE,eosRange,data1DForEosPtr,CENTER)
              end if
               dataForEosPtr(1,1,1,:) = dataForEos(:,1,1,1)
           end do
        end do
     end do
  end if

  deallocate(cellCenterCoord)

  return
contains
  pure real function gFactor(ri,delRi,alpha)
    real,intent(IN) :: ri,delRi,alpha
    real :: rim,rip,delh
    delh = 0.5 * delRi
    rim  = ri - delh
    rip  = ri + delh
    gFactor = 3.0*(rip**(3.0+alpha)-rim**(3.0+alpha)) &
         / ((3.0+alpha)*(rip**3 - rim**3))
  end function gFactor
  real function gF(i)
    integer,intent(IN) :: i
    gF = gFactor(cellCenterCoord(i), delR, alphaDens)
    !!print*,'gF(i) at i=',i,', r(i)=',cellCenterCoord(i),' is', gF
  end function gF
end subroutine Grid_bcApplyToRegionSpecialized
