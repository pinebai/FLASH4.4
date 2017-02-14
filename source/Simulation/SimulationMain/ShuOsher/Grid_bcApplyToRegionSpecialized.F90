!!****if* source/Simulation/SimulationMain/ShuOsher/Grid_bcApplyToRegionSpecialized
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
!!    regionSize     : regionSize(BC_DIR) contains the size of the each row and
!!                     regionSize(SECOND_DIR) contains the number of along the second
!!                     direction, and regionSize(THIRD_DIR) has the number of rows
!!                     along the third direction. regionSize(GRID_DATASTRUCT) contains the
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
!!          In other words, you can nearly always ignore this optional
!!          argument.  As of FLASH 3.0, it is only used internally within the
!!          Grid unit and is handled by the GridBoundaryConditions/Grid_bcApplyToRegion
!!          implementation. The argument has been added to the
!!          Grid_bcApplyToRegionSpecialized interface for consistency with
!!          Grid_bcApplyToRegion.
!!
!!
!! 2. ADDITIONAL ARGUMENTS
!!
!!  blockHandle - the local block number
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
!!
!!  blkLimitsGC - the starting and endpoint of the whole block including
!!                the guard cells, as returned by Grid_getBlkIndexLimits.
!!
!!  NOTES 
!!
!!            This routine is common to all the mesh packages supported.
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
     guard,axis,face,regionData,regionSize,mask,&
     applied,blockHandle,secondDir,thirdDir,endPoints,blkLimitsGC, idest)

#include "constants.h"
#include "Flash.h"

  use Simulation_data, ONLY : sim_rhoLeft, sim_rhoRight, &
                              sim_pLeft, sim_pRight, &
                              sim_uLeft, sim_uRight, &
                              sim_gamma, sim_aRho,  sim_fRho
  use Grid_data, ONLY :gr_meshMe,gr_domainBC
  use Driver_interface, ONLY : Driver_abortFlash
  use Driver_data, ONLY      : dr_simTime

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

  logical,dimension(regionSize(STRUCTSIZE)) :: localMask
  integer :: i,j, k,ivar,ind,je,ke,n,varCount
  real    :: sign
  real, allocatable, dimension(:) :: xCoord
  integer :: sizeX,istat,xOffset
  real    :: tempRho

  !! Set this true if user boundary condition is needed
  applied=.true.

  je=regionSize(SECOND_DIR)
  ke=regionSize(THIRD_DIR)

  varCount=regionSize(STRUCTSIZE)

  xOffset = endPoints(LOW,axis)

  !! Get x-Coordinate
  sizeX = blkLimitsGC(HIGH,axis)
  allocate(xCoord(sizeX))
  call gr_extendedGetCellCoords(axis, blockHandle, gr_meshMe, CENTER, .true., xCoord, sizeX)


   do ivar = 1,varCount
      if(mask(ivar)) then
         sign = 1

           if (gridDataStruct==CENTER) then
#ifdef VELX_VAR
              if (ivar==VELX_VAR)sign=-1.
#endif
#ifdef VELY_VAR
              if (ivar==VELY_VAR)sign=-1.
#endif
#ifdef VELZ_VAR
              if (ivar==VELZ_VAR)sign=-1.
#endif
           endif

        if(face==LOW) then ! apply BC on the bottom boundary at x=xmin

           select case (bcType)
           case(REFLECTING)
              k = 2*guard+1
              do i = 1,guard
                 regionData(i,1:je,1:ke,ivar)= regionData(k-i,1:je,1:ke,ivar)*sign
              end do

           case(OUTFLOW)
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

           case(USER_DEFINED)
              ! User BC at x=xmin
              k = 2*guard+1
              do i = 1,guard
                 do j=1,je
                    if (gridDataStruct==CENTER) then
                       if (ivar == DENS_VAR) regionData(i,j,1:ke,ivar) = sim_rhoLeft
                       if (ivar == VELX_VAR) regionData(i,j,1:ke,ivar) = sim_uLeft
                       if (ivar == VELY_VAR) regionData(i,j,1:ke,ivar) = 0.
                       if (ivar == VELZ_VAR) regionData(i,j,1:ke,ivar) = 0.
                       if (ivar == PRES_VAR) regionData(i,j,1:ke,ivar) = sim_pLeft
                       if (ivar == GAMC_VAR) regionData(i,j,1:ke,ivar) = sim_gamma
                       if (ivar == GAME_VAR) regionData(i,j,1:ke,ivar) = sim_gamma
                       if (ivar == EINT_VAR) regionData(i,j,1:ke,ivar) = &
                            sim_pLeft/(sim_rhoLeft*(sim_gamma-1.0))
                       if (ivar == ENER_VAR) regionData(i,j,1:ke,ivar) = &
                            sim_pLeft/(sim_rhoLeft*(sim_gamma-1.0)) + 0.5*sim_uLeft**2
                    endif
                 end do
              end do
           case default
              print*,'boundary is',bcType
              call Driver_abortFlash("unsupported boundary condition on Lower Face")
           end select

        else !apply BC on the top boundary at x=xmax

           select case (bcType)
           case(REFLECTING)
              k = 2*guard+1
              do i = 1,guard
                 regionData(k-i,1:je,1:ke,ivar)= regionData(i,1:je,1:ke,ivar)*sign
              end do

           case(OUTFLOW)
              k=guard
              do i = 1,guard
                 regionData(k+i,1:je,1:ke,ivar)= regionData(k,1:je,1:ke,ivar)
              end do

           case(DIODE)
              k=guard
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

           case(USER_DEFINED)
              ! User BC at x=xmax
              k = guard
              do i = 1,guard
                 do j=1,je
                    if (gridDataStruct==CENTER) then
                       tempRho = sim_rhoRight*(1.0+sim_aRho*sin(sim_fRho*xCoord(i+k+xOffset-1)))
                       if (ivar == DENS_VAR) regionData(k+i,j,1:ke,ivar) = tempRho
                       if (ivar == VELX_VAR) regionData(k+i,j,1:ke,ivar) = sim_uRight
                       if (ivar == VELY_VAR) regionData(k+i,j,1:ke,ivar) = 0.
                       if (ivar == VELZ_VAR) regionData(k+i,j,1:ke,ivar) = 0.
                       if (ivar == PRES_VAR) regionData(k+i,j,1:ke,ivar) = sim_pRight
                       if (ivar == GAMC_VAR) regionData(k+i,j,1:ke,ivar) = sim_gamma
                       if (ivar == GAME_VAR) regionData(k+i,j,1:ke,ivar) = sim_gamma
                       if (ivar == EINT_VAR) regionData(k+i,j,1:ke,ivar) = &
                            sim_pRight/(tempRho*(sim_gamma-1.0))
                       if (ivar == ENER_VAR) regionData(k+i,j,1:ke,ivar) = &
                            sim_pRight/(tempRho*(sim_gamma-1.0)) + 0.5*sim_uRight**2
                    endif
                 end do
              end do
           case default
              print*,'boundary is',bcType
              call Driver_abortFlash("unsupported boundary condition on Upper Face")
           end select

        end if
     end if
  end do

  !! deallocate
  deallocate(xCoord)


  return

end subroutine Grid_bcApplyToRegionSpecialized
