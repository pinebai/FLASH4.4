!!****if* source/Simulation/SimulationMain/FlameChannel/Grid_bcApplyToRegionSpecialized
!!
!! NAME
!!  Grid_bcApplyToRegionSpecialized
!!
!! SYNOPSIS
!!
!!  call Grid_bcApplyToRegionSpecialized(integer(IN)  :: bcType,
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
!!  exchanging data with each other.
!!  This routine currently passes handling of all other boundary condition
!!  types on to Grid_applyBCEdge, which is called for each of the variables
!!  in turn.
!!
!!  This routine supports only simple boundary conditions that are applied strictly
!!  directionally, and have no need for other grid information such as the coordinates etc.
!!  Currently supported boundary conditions are "OUTFLOW", "REFLECTING" and "DIODE".
!!  The "PERIODIC" boundary conditions are automatically applied in the process of filling
!!  the guard cells all over the domain, and therefore do not need to call this routine
!!
!!  If the user wishes to apply different boundary conditions, they can either
!!  use the interface Grid_bcApplyToRegionSpecializedSpecialized, or make a copy of this
!!  routine in their simulation directory and customize it.
!!
!!
!! ARGUMENTS 
!!
!!  bcType - the type of boundary conditions being applied
!!  gridDataStruct - the Grid dataStructure
!!  guard -    number of guard cells 
!!    axis  - the direction along which to apply boundary conditions,
!!          can take values of IAXIS, JAXIS and KAXIS
!!  face    -  can take values LOW and HIGH, defined in constants.h,
!!             to indicate whether to apply boundary on lowerface or 
!!             upperface
!!  regionData     : the extracted region from a block of permanent storage of the 
!!                   specified data structure. Its size is given by regionSize.
!!                   NOTE that the first three dimensions of this array do not necessarily
!!                   correspond to the (IAXIS, JAXIS, KAXIS) directions in this order;
!!                   rather, the axes are permuted such that the first index
!!                   of regionData always corresponds to the direction given by axis.
!!                   See regionSize for more information.
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
!!  applied - is set true if this routine has handled the given bcType, otherwise it is 
!!            set to false.
!!
!!  blockHandle - Handle for the block for which guard cells are to be filled.
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
!!                         such as REFLECTING or OUTFLOW, It is provided for
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
!!  HISTORY
!!
!!    2007       Initial Grid_bcApplyToRegionSpecialized logic    - Anshu Dubey
!!    2007       Cleanup and fixes                     - Klaus Weide, Dongwook Lee 
!!    2008       Tweaks to support Multigrid           - Klaus Weide
!!
!!***

subroutine Grid_bcApplyToRegionSpecialized(bcType,gridDataStruct,&
     guard,axis,face,regionData,regionSize,mask,applied,&
     blockHandle,secondDir,thirdDir,endPoints,blkLimitsGC, idest)

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  use Simulation_data, ONLY : sim_eosData_u, sim_eosData_b, sim_inflowVortex, &
                              sim_inflowVx, sim_cFrac, sim_neFrac, sim_sigT, &
                              sim_sigP, sim_sigVx, sim_sigVy, sim_sigVz, &
                              sim_ymin, sim_ymax, sim_zmin, sim_zmax, &
                              sim_vortexSize, sim_vortexStrength, &
                              sim_yctrVortex
  use Eos_interface, ONLY : Eos
  use Eos_data, ONLY : eos_maxNewton, eos_tol
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_bcInterface, ONLY : gr_bcMapBcType
  use Grid_data, ONLY : gr_meshMe

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
  integer :: startBC, endBC, startInt, endInt, interior
  logical :: isFace
  integer :: sign
  real :: kine, velx, vely, velz, vel2
  real :: vy_l, vy_c, vy_r, vz_l, vz_c, vz_r
  real :: target_vely, target_velz
  real :: y_vort1, y_vort2, z_vort, vortex_stream1, vortex_stream2
  real :: vortexSize, vortexStrength

  integer :: iter, iter2
  real, dimension(EOS_NUM) :: eosData, eosData_new
  real :: rho_guess, temp_guess, p_guess, eint_guess
  real :: rho_target, temp_target, p_target, eint_target
  real :: err, err2
  logical, dimension(EOS_VARS+1:EOS_NUM) :: eosMask

  real, allocatable, dimension(:) :: secondDirLeft, thirdDirLeft
  real, allocatable, dimension(:) :: secondDirRight, thirdDirRight
  real, allocatable, dimension(:) :: secondDirCenter, thirdDirCenter

  select case (bcType)
  case(USER_DEFINED)
     if ((gridDataStruct==CENTER).and.(axis==IAXIS)) then
        applied = .TRUE.           !will handle this types below
     else
        print*,"User BCs called for non-centered data or non x-axis"
        applied = .FALSE. 
     end if
!  case(OUTFLOW)
!     if ((gridDataStruct==CENTER).and.(axis==IAXIS).and.(face==HIGH)) then
!        applied = .TRUE.
!     else
!        applied = .FALSE.
!     endif
  case(REFLECTING)
     if ((gridDataStruct==CENTER).and.(axis==IAXIS).and.(face==LOW)) then
        applied = .TRUE.
     else
        applied = .FALSE.
     endif
  case default
     applied = .FALSE.
  end select
  if (.not. applied) return   !RETURN immediately!

  je=regionSize(SECOND_DIR)
  ke=regionSize(THIRD_DIR)
  varCount=regionSize(STRUCTSIZE)
  isFace = (gridDataStruct==FACEX).and.(axis==IAXIS)
  isFace = isFace.or.((gridDataStruct==FACEY).and.(axis==JAXIS))
  isFace = isFace.or.((gridDataStruct==FACEZ).and.(axis==KAXIS))

  ! get coordinates
  ! for axis == IAXIS, secondDir = JAXIS and thirdDir = KAXIS
  allocate(secondDirLeft( blkLimitsGC(HIGH,secondDir) ))
  allocate(thirdDirLeft( blkLimitsGC(HIGH,thirdDir) ))
  allocate(secondDirCenter( blkLimitsGC(HIGH,secondDir) ))
  allocate(thirdDirCenter( blkLimitsGC(HIGH,thirdDir) ))
  allocate(secondDirRight( blkLimitsGC(HIGH,secondDir) ))
  allocate(thirdDirRight( blkLimitsGC(HIGH,thirdDir) ))
  call gr_extendedGetCellCoords(secondDir, blockHandle, gr_meshMe, LEFT_EDGE, &
     .true., secondDirLeft, blkLimitsGC(HIGH,secondDir) )
  call gr_extendedGetCellCoords(thirdDir, blockHandle, gr_meshMe, LEFT_EDGE, &
     .true., thirdDirLeft, blkLimitsGC(HIGH,thirdDir) )
  call gr_extendedGetCellCoords(secondDir, blockHandle, gr_meshMe, CENTER, &
     .true., secondDirCenter, blkLimitsGC(HIGH,secondDir) )
  call gr_extendedGetCellCoords(thirdDir, blockHandle, gr_meshMe, CENTER, &
     .true., thirdDirCenter, blkLimitsGC(HIGH,thirdDir) )
  call gr_extendedGetCellCoords(secondDir, blockHandle, gr_meshMe, RIGHT_EDGE, &
     .true., secondDirRight, blkLimitsGC(HIGH,secondDir) )
  call gr_extendedGetCellCoords(thirdDir, blockHandle, gr_meshMe, RIGHT_EDGE, &
     .true., thirdDirRight, blkLimitsGC(HIGH,thirdDir) )

!!  print*,'in applyBcRegion ',varCount,gridDataStruct,guard,axis,face

   eosMask(:) = .false.
   eosMask(EOS_DPD) = .true.
   eosMask(EOS_DED) = .true.
   eosMask(EOS_DPT) = .true.

  ! set bounds on regionData
  if (face==LOW) then
     startBC = 1
     endBC = guard
     interior = endBC + 1
     startInt = interior
     endInt = regionSize(BC_DIR)
  else
     startBC = regionSize(BC_DIR) - guard + 1
     endBC = regionSize(BC_DIR)
     interior = startBC - 1
     startInt = 1
     endInt = interior
  endif

  ! initialize regionData
  do ivar = 1,varCount
     if (mask(ivar)) then
        call gr_bcMapBcType(bcTypeActual,bcType,ivar,gridDataStruct,axis,face,idest)
        if (bcTypeActual /= bcType) &
           call Driver_abortFlash("[Grid_bcApplyToRegionSpecialized] unexpected behavior")
        ! default is zero-gradient
        do i = startBC, endBC
           regionData(i,1:je,1:ke,ivar) = regionData(interior,1:je,1:ke,ivar)
        end do
     endif
  end do

  select case (bcType)
  case(USER_DEFINED)
  ! implement static pressure outflow
  if (face==LOW) then
     do k = 1, ke
        do j = 1, je

           eosData(:) = sim_eosData_b(:)
           ! guess initial density
           rho_guess = eosData(EOS_DENS)
           ! set p at infinity to initial ash state
           p_target = eosData(EOS_PRES)
           !p_target = regionData(interior,j,k,PRES_VAR)

           ! conserve energy (assuming zero-gradient y-z velocities)
           eint_target = regionData(interior,j,k,EINT_VAR) + &
              regionData(interior,j,k,PRES_VAR) / regionData(interior,j,k,DENS_VAR) - &
              p_target / rho_guess + 0.5 * regionData(interior,j,k,VELX_VAR)**2 * &
              ( 1.0 - regionData(interior,j,k,DENS_VAR) / rho_guess )
           err = abs(eint_target - eosData(EOS_EINT)) / eint_target
           !print*,"Looping over BC LOW (j, k) ",j,k
           !print*,"PRES,DENS,EINT",p_target,rho_guess,eint_target

           ! the double loop is necessary because using EOS mode MODE_DENS_PRES
           ! does not converge.
           iter = 0
           do while ( err > eos_tol )
              iter = iter + 1

              do iter2 = 1, eos_maxNewton

                 eosData(EOS_EINT) = eint_target
                 eosData(EOS_DENS) = rho_guess
                 ! assume inflow is fuel
                 eosData(EOS_ABAR) = sim_eosData_b(EOS_ABAR)
                 eosData(EOS_ZBAR) = sim_eosData_b(EOS_ZBAR)
                 call Eos(MODE_DENS_EI,1, eosData, mask=eosMask)

                 rho_guess = min(10.0*eosData(EOS_DENS), max(0.1*eosData(EOS_DENS), &
                    eosData(EOS_DENS) + (p_target - eosData(EOS_PRES)) / &
                    eosData(EOS_DPD)))

                 err2 = abs(rho_guess - eosData(EOS_DENS)) / rho_guess

                 if (err2 < eos_tol) exit
              end do

              if (err2 >= eos_tol) then
                 iter = eos_maxNewton + 1
                 exit
              endif

              ! conserve energy (assuming zero-gradient y-z velocities)
              eint_target = regionData(interior,j,k,EINT_VAR) + &
                 regionData(interior,j,k,PRES_VAR) / regionData(interior,j,k,DENS_VAR) - &
                 p_target / rho_guess + 0.5 * regionData(interior,j,k,VELX_VAR)**2 * &
                 ( 1.0 - regionData(interior,j,k,DENS_VAR) / rho_guess )

              err = abs(eint_target - eosData(EOS_EINT)) / eint_target

              if (iter > eos_maxNewton) exit
           end do

           if (iter > eos_maxNewton) &
              call Driver_abortFlash("[Grid_bcApplyToRegionSpecialized] Low face: Too many iterations")

           eosData(EOS_EINT) = eint_target
           eosData(EOS_DENS) = rho_guess
           eosData(EOS_ABAR) = sim_eosData_b(EOS_ABAR)
           eosData(EOS_ZBAR) = sim_eosData_b(EOS_ZBAR)
           call Eos(MODE_DENS_EI,1, eosData)

           ! conserve momentum
           velx = regionData(interior,j,k,DENS_VAR) * &
                  regionData(interior,j,k,VELX_VAR) / &
                  eosData(EOS_DENS)

           ! assume density doesn't change in y-z plane
           vely = regionData(interior,j,k,VELY_VAR)
           velz = regionData(interior,j,k,VELZ_VAR)
           kine = 0.5 * ( velx**2 + vely**2 + velz**2 )
           
           ! fill in solution
           do i = startBC, endBC
              regionData(i,j,k,DENS_VAR) = eosData(EOS_DENS)
              regionData(i,j,k,TEMP_VAR) = eosData(EOS_TEMP)
              regionData(i,j,k,PRES_VAR) = eosData(EOS_PRES)
              regionData(i,j,k,VELX_VAR) = velx
              regionData(i,j,k,VELY_VAR) = vely
              regionData(i,j,k,VELZ_VAR) = velz
              regionData(i,j,k,ENER_VAR) = eosData(EOS_EINT) + kine
              regionData(i,j,k,EINT_VAR) = eosData(EOS_EINT)
              regionData(i,j,k,GAMC_VAR) = eosData(EOS_GAMC)
              regionData(i,j,k,GAME_VAR) = eosData(EOS_PRES) / &
                 (eosData(EOS_EINT)*eosData(EOS_DENS)) + 1.0
           end do

        end do
     end do

  else ! face==HIGH ! SUBSONIC INFLOW

     do k = 1, ke
        do j = 1, je
           ! simplest subsonic inflow specifies velocities, mass fractions, and
           ! temperature while allows the density or pressure to be solved.
           ! non-reflecting subsonic inflow require variable velocities and/or
           ! temperature
           ! The current implementation assumes planar flame. Should revisit and
           ! consider non-zero derivatives in velocity and pressure in y-z
           ! plane.
           if (sim_inflowVortex) then

              ! allow temperature to vary, but damp it
              ! if sim_sigT = 0, then non-reflecting
              ! if sim_sigT = 1, then target maintained
              temp_target = (1.0-sim_sigT) * regionData(interior,j,k,TEMP_VAR) + &
                 sim_sigT * sim_eosData_u(EOS_TEMP)

              ! if this is steady-state 1D problem, then mass continuity
              ! and momentum conservation indicate that the pressure should
              ! be zero-gradient
              ! set p to match interior
              p_target = regionData(interior,j,k,PRES_VAR)
              ! initial guess at density
              rho_guess = regionData(interior,j,k,DENS_VAR)

              do iter = 1, eos_maxNewton

                 eosData(EOS_TEMP) = temp_target
                 eosData(EOS_DENS) = rho_guess
                 ! assume inflow is fuel
                 eosData(EOS_ABAR) = sim_eosData_u(EOS_ABAR)
                 eosData(EOS_ZBAR) = sim_eosData_u(EOS_ZBAR)
                 call Eos(MODE_DENS_TEMP,1, eosData, mask=eosMask)

                 rho_guess = min(10.0*eosData(EOS_DENS), max(0.1*eosData(EOS_DENS), &
                    eosData(EOS_DENS) + (p_target - eosData(EOS_PRES)) / &
                    eosData(EOS_DPD)))

                 err = abs(rho_guess - eosData(EOS_DENS)) / rho_guess

                 if (err < eos_tol) exit
              end do

              if (err >= eos_tol) &
                 call Driver_abortFlash("[Grid_bcApplyToRegionSpecialized] High face: Too many iterations")

              eosData(EOS_TEMP) = temp_target
              eosData(EOS_DENS) = rho_guess
              eosData(EOS_ABAR) = sim_eosData_u(EOS_ABAR)
              eosData(EOS_ZBAR) = sim_eosData_u(EOS_ZBAR)
              call Eos(MODE_DENS_TEMP,1, eosData)

!              ! allow velocity to vary, but damp it
!              ! if sim_sigV = 0, then non-reflecting
!              ! if sim_sigV = 1, then target maintained
!              ! get new velx from mass continuity
!              ! compute left edge
!              y_vort1 = secondDirLeft(j) - sim_yctrVortex
!              ! mirror about the periodic boundary
!              y_vort2 = secondDirLeft(j) - &
!                 ( sim_ymin - sim_yctrVortex + sim_ymax )
!              ! center of domain
!              z_vort = thirdDirLeft(k) - 0.5*(sim_zmax+sim_zmin)
!              vortex_stream1 = sim_vortexStrength * &
!                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!              ! counter-rotating, so make minus
!              vortex_stream2 = -sim_vortexStrength * &
!                 exp( -( y_vort2**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!
!              ! vy = d/dz(stream), vz = -d/dy(stream)
!              vy_l = -z_vort / sim_vortexSize**2 * &
!                 ( vortex_stream1 + vortex_stream2 )
!              vz_l = 1.0 / sim_vortexSize**2 * &
!                 ( y_vort1*vortex_stream1 + y_vort2*vortex_stream2 )
!
!              ! compute right edge
!              y_vort1 = secondDirRight(j) - sim_yctrVortex
!              ! mirror about the periodic boundary
!              y_vort2 = secondDirRight(j) - &
!                 ( sim_ymin - sim_yctrVortex + sim_ymax )
!              ! center of domain
!              z_vort = thirdDirRight(k) - 0.5*(sim_zmax+sim_zmin)
!              vortex_stream1 = sim_vortexStrength * &
!                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!              ! counter-rotating, so make minus
!              vortex_stream2 = -sim_vortexStrength * &
!                 exp( -( y_vort2**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!
!              vy_r = -z_vort / sim_vortexSize**2 * &
!                 ( vortex_stream1 + vortex_stream2 )
!              vz_r = 1.0 / sim_vortexSize**2 * &
!                 ( y_vort1*vortex_stream1 + y_vort2*vortex_stream2 )
!
!              ! compute cell-center
!              y_vort1 = secondDirCenter(j) - sim_yctrVortex
!              ! mirror about the periodic boundary
!              y_vort2 = secondDirCenter(j) - &
!                 ( sim_ymin - sim_yctrVortex + sim_ymax )
!              ! center of domain
!              z_vort = thirdDirCenter(k) - 0.5*(sim_zmax+sim_zmin)
!              vortex_stream1 = sim_vortexStrength * &
!                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!              ! counter-rotating, so make minus
!              vortex_stream2 = -sim_vortexStrength * &
!                 exp( -( y_vort2**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!
!              vy_c = -z_vort / sim_vortexSize**2 * &
!                 ( vortex_stream1 + vortex_stream2 )
!              vz_c = 1.0 / sim_vortexSize**2 * &
!                 ( y_vort1*vortex_stream1 + y_vort2*vortex_stream2 )

!** Only one vortex
              y_vort1 = secondDirLeft(j) - sim_yctrVortex
              z_vort = thirdDirLeft(k) - 0.5*(sim_zmax+sim_zmin)
              vortex_stream1 = sim_vortexStrength * &
                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )

              ! vy = d/dz (stream), vz = -d/dy (stream)
              vy_l = -z_vort / sim_vortexSize**2 * vortex_stream1
              vz_l = y_vort1 / sim_vortexSize**2 * vortex_stream1

              y_vort1 = secondDirRight(j) - sim_yctrVortex
              z_vort = thirdDirRight(k) - 0.5*(sim_zmax+sim_zmin)
              vortex_stream1 = sim_vortexStrength * &
                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )

              ! vy = d/dz (stream), vz = -d/dy (stream)
              vy_r = -z_vort / sim_vortexSize**2 * vortex_stream1
              vz_r = y_vort1 / sim_vortexSize**2 * vortex_stream1

              y_vort1 = secondDirCenter(j) - sim_yctrVortex
              z_vort = thirdDirCenter(k) - 0.5*(sim_zmax+sim_zmin)
              vortex_stream1 = sim_vortexStrength * &
                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )

              ! vy = d/dz (stream), vz = -d/dy (stream)
              vy_c = -z_vort / sim_vortexSize**2 * vortex_stream1
              vz_c = y_vort1 / sim_vortexSize**2 * vortex_stream1

              ! compute cell-average (simplified since density is constant)
              target_vely = (vy_l + 4.0*vy_c + vy_r) / 6.0
              target_velz = (vz_l + 4.0*vz_c + vz_r) / 6.0

              velx = regionData(interior,j,k,VELX_VAR) * &
                 regionData(interior,j,k,DENS_VAR) / eosData(EOS_DENS)
              velx = (1.0-sim_sigVx) * velx + sim_sigVx * sim_inflowVx
              vely = (1.0-sim_sigVy) * regionData(interior,j,k,VELY_VAR) + &
                 sim_sigVy * target_vely
              velz = (1.0-sim_sigVz) * regionData(interior,j,k,VELZ_VAR) + &
                 sim_sigVz * target_velz

              kine = 0.5 * ( velx**2 + vely**2 + velz**2 )
              ! fill in solution
              do i = startBC, endBC
                 regionData(i,j,k,DENS_VAR) = eosData(EOS_DENS)
                 regionData(i,j,k,TEMP_VAR) = eosData(EOS_TEMP)
                 regionData(i,j,k,PRES_VAR) = eosData(EOS_PRES)
                 regionData(i,j,k,VELX_VAR) = velx
                 regionData(i,j,k,VELY_VAR) = vely
                 regionData(i,j,k,VELZ_VAR) = velz
                 regionData(i,j,k,ENER_VAR) = eosData(EOS_EINT) + kine
                 regionData(i,j,k,EINT_VAR) = eosData(EOS_EINT)
                 regionData(i,j,k,GAMC_VAR) = eosData(EOS_GAMC)
                 regionData(i,j,k,GAME_VAR) = eosData(EOS_PRES) / &
                    (eosData(EOS_EINT)*eosData(EOS_DENS)) + 1.0
              end do

           else
              ! allow temperature to vary, but damp it
              ! if sim_sigT = 0, then non-reflecting
              ! if sim_sigT = 1, then target maintained
              temp_target = (1.0-sim_sigT) * regionData(interior,j,k,TEMP_VAR) + &
                 sim_sigT * sim_eosData_u(EOS_TEMP)

              ! if this is steady-state 1D problem, then mass continuity
              ! and momentum conservation indicate that the pressure should
              ! be zero-gradient
              ! set p to match interior
              p_target = regionData(interior,j,k,PRES_VAR)
              ! initial guess at density
              rho_guess = regionData(interior,j,k,DENS_VAR)

              do iter = 1, eos_maxNewton

                 eosData(EOS_TEMP) = temp_target
                 eosData(EOS_DENS) = rho_guess
                 ! assume inflow is fuel
                 eosData(EOS_ABAR) = sim_eosData_u(EOS_ABAR)
                 eosData(EOS_ZBAR) = sim_eosData_u(EOS_ZBAR)
                 call Eos(MODE_DENS_TEMP,1, eosData, mask=eosMask)

                 rho_guess = min(10.0*eosData(EOS_DENS), max(0.1*eosData(EOS_DENS), &
                    eosData(EOS_DENS) + (p_target - eosData(EOS_PRES)) / &
                    eosData(EOS_DPD)))

                 err = abs(rho_guess - eosData(EOS_DENS)) / rho_guess

                 if (err < eos_tol) exit
              end do

              if (err >= eos_tol) &
                 call Driver_abortFlash("[Grid_bcApplyToRegionSpecialized] High face: Too many iterations")

              eosData(EOS_TEMP) = temp_target
              eosData(EOS_DENS) = rho_guess
              eosData(EOS_ABAR) = sim_eosData_u(EOS_ABAR)
              eosData(EOS_ZBAR) = sim_eosData_u(EOS_ZBAR)
              call Eos(MODE_DENS_TEMP,1, eosData)

              ! allow velocity to vary, but damp it
              ! if sim_sigV = 0, then non-reflecting
              ! if sim_sigV = 1, then target maintained
              ! get new velx from mass continuity
              velx = regionData(interior,j,k,VELX_VAR) * &
                 regionData(interior,j,k,DENS_VAR) / eosData(EOS_DENS)
              velx = (1.0-sim_sigVx) * velx + sim_sigVx * sim_inflowVx
              vely = (1.0-sim_sigVy) * regionData(interior,j,k,VELY_VAR) ! +0
              velz = (1.0-sim_sigVz) * regionData(interior,j,k,VELZ_VAR) ! +0
              kine = 0.5 * ( velx**2 + vely**2 + velz**2 )

              ! fill in solution
              do i = startBC, endBC
                 regionData(i,j,k,DENS_VAR) = eosData(EOS_DENS)
                 regionData(i,j,k,TEMP_VAR) = eosData(EOS_TEMP)
                 regionData(i,j,k,PRES_VAR) = eosData(EOS_PRES)
                 regionData(i,j,k,VELX_VAR) = velx
                 regionData(i,j,k,VELY_VAR) = vely
                 regionData(i,j,k,VELZ_VAR) = velz
                 regionData(i,j,k,ENER_VAR) = eosData(EOS_EINT) + kine
                 regionData(i,j,k,EINT_VAR) = eosData(EOS_EINT)
                 regionData(i,j,k,GAMC_VAR) = eosData(EOS_GAMC)
                 regionData(i,j,k,GAME_VAR) = eosData(EOS_PRES) / &
                    (eosData(EOS_EINT)*eosData(EOS_DENS)) + 1.0
              end do

           end if

        end do
     end do

  endif

  ! modify the outflow to potentially specify a desired pressure
  case(OUTFLOW)

     ! zero-gradient by default
     if ( sim_sigP > 0.0 ) then

        ! implement non-reflecting outflow
        do k = 1, ke
           do j = 1, je
              ! guess initial density
              rho_guess = regionData(interior,j,k,DENS_VAR)
              ! set p at infinity to initial fuel state
              p_target = (1.0-sim_sigP) * regionData(interior,j,k,PRES_VAR) + &
                 sim_sigP * sim_eosData_u(EOS_PRES)

              ! conserve energy (assuming zero-gradient y-z velocities)
              eint_target = regionData(interior,j,k,EINT_VAR) + &
                 regionData(interior,j,k,PRES_VAR) / regionData(interior,j,k,DENS_VAR) - &
                 p_target / rho_guess + 0.5 * regionData(interior,j,k,VELX_VAR)**2 * &
                 ( 1.0 - regionData(interior,j,k,DENS_VAR) / rho_guess )

              ! the double loop is necessary because using EOS mode MODE_DENS_PRES
              ! does not converge.
              do iter = 1, eos_maxNewton

                 do iter2 = 1, eos_maxNewton

                    eosData(EOS_EINT) = eint_target
                    eosData(EOS_DENS) = rho_guess
                    ! assume inflow is fuel
                    eosData(EOS_ABAR) = sim_eosData_u(EOS_ABAR)
                    eosData(EOS_ZBAR) = sim_eosData_u(EOS_ZBAR)
                    call Eos(MODE_DENS_EI,1, eosData, mask=eosMask)

                    rho_guess = min(10.0*eosData(EOS_DENS), max(0.1*eosData(EOS_DENS), &
                       eosData(EOS_DENS) + (p_target - eosData(EOS_PRES)) / &
                       eosData(EOS_DPD)))

                    err2 = abs(rho_guess - eosData(EOS_DENS)) / rho_guess

                    if (err2 < eos_tol) exit
                 end do

                 if (err2 >= eos_tol) then
                    err = err2
                    exit
                 endif

                 ! conserve energy (assuming zero-gradient y-z velocities)
                 eint_target = regionData(interior,j,k,EINT_VAR) + &
                    regionData(interior,j,k,PRES_VAR) / regionData(interior,j,k,DENS_VAR) - &
                    p_target / rho_guess + 0.5 * regionData(interior,j,k,VELX_VAR)**2 * &
                    ( 1.0 - regionData(interior,j,k,DENS_VAR) / rho_guess )

                 err = abs(eint_target - eosData(EOS_EINT)) / eint_target

                 if (err < eos_tol) exit
              end do

              if (err >= eos_tol) &
                 call Driver_abortFlash("[Grid_bcApplyToRegionSpecialized] High face: Too many iterations")

              eosData(EOS_EINT) = eint_target
              eosData(EOS_DENS) = rho_guess
              eosData(EOS_ABAR) = sim_eosData_u(EOS_ABAR)
              eosData(EOS_ZBAR) = sim_eosData_u(EOS_ZBAR)
              call Eos(MODE_DENS_EI,1, eosData)

              ! conserve momentum
              velx = regionData(interior,j,k,DENS_VAR) * &
                     regionData(interior,j,k,VELX_VAR) / &
                     eosData(EOS_DENS)

              ! assume density doesn't change in y-z plane
              vely = regionData(interior,j,k,VELY_VAR)
              velz = regionData(interior,j,k,VELZ_VAR)
              kine = 0.5 * ( velx**2 + vely**2 + velz**2 )

              ! fill in solution
              do i = startBC, endBC
                 regionData(i,j,k,DENS_VAR) = eosData(EOS_DENS)
                 regionData(i,j,k,TEMP_VAR) = eosData(EOS_TEMP)
                 regionData(i,j,k,PRES_VAR) = eosData(EOS_PRES)
                 regionData(i,j,k,VELX_VAR) = velx
                 regionData(i,j,k,VELY_VAR) = vely
                 regionData(i,j,k,VELZ_VAR) = velz
                 regionData(i,j,k,ENER_VAR) = eosData(EOS_EINT) + kine
                 regionData(i,j,k,EINT_VAR) = eosData(EOS_EINT)
                 regionData(i,j,k,GAMC_VAR) = eosData(EOS_GAMC)
                 regionData(i,j,k,GAME_VAR) = eosData(EOS_PRES) / &
                    (eosData(EOS_EINT)*eosData(EOS_DENS)) + 1.0
              end do

           end do
        end do

     endif

  case (REFLECTING)
     ! implement reflecting with zero velocities (velocity nodes)
     ! and specified temperature 
     do k = 1, ke
        do j = 1, je

           eosData(:) = sim_eosData_b(:)
           ! guess initial density
           rho_guess = eosData(EOS_DENS)
           ! set p at infinity to initial ash state
           p_target = eosData(EOS_PRES)
           !p_target = regionData(interior,j,k,PRES_VAR)

           ! conserve energy (assuming zero-gradient y-z velocities)
           eint_target = regionData(interior,j,k,EINT_VAR) + &
              regionData(interior,j,k,PRES_VAR) / regionData(interior,j,k,DENS_VAR) - &
              p_target / rho_guess + 0.5 * regionData(interior,j,k,VELX_VAR)**2
           err = abs(eint_target - eosData(EOS_EINT)) / eint_target
           !print*,"Looping over BC LOW (j, k) ",j,k
           !print*,"PRES,DENS,EINT",p_target,rho_guess,eint_target

           ! the double loop is necessary because using EOS mode MODE_DENS_PRES
           ! does not converge.
           iter = 0
           do while ( err > eos_tol )
              iter = iter + 1

              do iter2 = 1, eos_maxNewton

                 eosData(EOS_EINT) = eint_target
                 eosData(EOS_DENS) = rho_guess
                 ! assume inflow is fuel
                 eosData(EOS_ABAR) = sim_eosData_b(EOS_ABAR)
                 eosData(EOS_ZBAR) = sim_eosData_b(EOS_ZBAR)
                 call Eos(MODE_DENS_EI,1, eosData, mask=eosMask)

                 rho_guess = min(10.0*eosData(EOS_DENS), max(0.1*eosData(EOS_DENS), &
                    eosData(EOS_DENS) + (p_target - eosData(EOS_PRES)) / &
                    eosData(EOS_DPD)))

                 err2 = abs(rho_guess - eosData(EOS_DENS)) / rho_guess

                 if (err2 < eos_tol) exit
              end do

              if (err2 >= eos_tol) then
                 iter = eos_maxNewton + 1
                 exit
              endif

              ! conserve energy (assuming zero-gradient y-z velocities)
              eint_target = regionData(interior,j,k,EINT_VAR) + &
                 regionData(interior,j,k,PRES_VAR) / regionData(interior,j,k,DENS_VAR) - &
                 p_target / rho_guess + 0.5 * regionData(interior,j,k,VELX_VAR)**2 

              err = abs(eint_target - eosData(EOS_EINT)) / eint_target

              if (iter > eos_maxNewton) exit
           end do

           if (iter > eos_maxNewton) &
              call Driver_abortFlash("[Grid_bcApplyToRegionSpecialized] Low face: Too many iterations")

           eosData(EOS_EINT) = eint_target
           eosData(EOS_DENS) = rho_guess
           eosData(EOS_ABAR) = sim_eosData_b(EOS_ABAR)
           eosData(EOS_ZBAR) = sim_eosData_b(EOS_ZBAR)
           call Eos(MODE_DENS_EI,1, eosData)

           ! don't allow mass to flow off the grid
           velx = 0.0 
           vely = regionData(interior,j,k,VELY_VAR)
           velz = regionData(interior,j,k,VELZ_VAR)
           kine = 0.5 * (velx**2 + vely**2 + velz**2)
        
           ! fill in solution
           do i = startBC, endBC
              regionData(i,j,k,DENS_VAR) = eosData(EOS_DENS)
              regionData(i,j,k,TEMP_VAR) = eosData(EOS_TEMP)
              regionData(i,j,k,PRES_VAR) = eosData(EOS_PRES)
              regionData(i,j,k,VELX_VAR) = velx
              regionData(i,j,k,VELY_VAR) = vely
              regionData(i,j,k,VELZ_VAR) = velz
              regionData(i,j,k,ENER_VAR) = eosData(EOS_EINT) + kine
              regionData(i,j,k,EINT_VAR) = eosData(EOS_EINT)
              regionData(i,j,k,GAMC_VAR) = eosData(EOS_GAMC)
              regionData(i,j,k,GAME_VAR) = eosData(EOS_PRES) / &
                 (eosData(EOS_EINT)*eosData(EOS_DENS)) + 1.0
           end do

        end do
     end do
 
  end select ! bcType 

  deallocate(secondDirLeft)
  deallocate(thirdDirLeft)
  deallocate(secondDirCenter)
  deallocate(thirdDirCenter)
  deallocate(secondDirRight)
  deallocate(thirdDirRight)

  return
end subroutine Grid_bcApplyToRegionSpecialized
