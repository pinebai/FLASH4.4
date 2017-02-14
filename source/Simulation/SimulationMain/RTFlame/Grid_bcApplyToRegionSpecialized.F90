! Dean Townsley 2009
!
! See source/Grid/GridBoundaryConditions/OneRow/Grid_bcApplyToRegionSpecialized.F90
!     for API documentation
!
! This version is specialized for hydrostatic equilibrium boundary conditions
! similar to those discussed in Zingale et al (2002ApJS..143..539Z)

subroutine Grid_bcApplyToRegionSpecialized(bcType,gridDataStruct,&
     guard,axis,face,regionData,regionSize,mask,applied,&
     blockHandle,secondDir,thirdDir,endPoints,blkLimitsGC, idest)

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  use Grid_data, ONLY :gr_meshMe,gr_domainBC
  use Driver_interface, ONLY : Driver_abortFlash
  use hse_interface
  use Simulation_data, ONLY : sim_grav, HSE_FORWARD, HSE_BACKWARD, HSE_SETTEMP
  use Eos_interface, ONLY : Eos

  implicit none

  ! this subroutine is define OUTSIDE of the Eos Unit
  ! this is the callback for whoever requested the externalAbarZbar implementation
  interface eos_externalComputeAbarZbar
     subroutine eos_externalComputeAbarZbar( solnScalars, abarData, zbarData)
        implicit none
        real, dimension(:,:), intent(in) :: solnScalars
        real, dimension(:), intent(out)  :: abarData, zbarData
     end subroutine
  end interface


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

  integer :: i,j,k, sizeGC, start, end, step, direction, lowi, hii
  real, allocatable, dimension(:) :: cellCenterCoord
  real    :: deltax
  real, dimension(EOS_NUM) :: eosData

  real, allocatable, dimension(:,:) :: scalars
  real, allocatable, dimension(:) :: region_abar, region_zbar, region_sumy, region_ye


!=====================================================================

  ! we are going to assume that gravity is along the IAXIS direction
  ! apply our hydrostatic in both directions
  if (axis/=IAXIS) then
     applied = .false.
     return
  endif
  applied = .true.

  ! get coordinates
  sizeGC = blkLimitsGC(HIGH,axis)
  allocate(cellCenterCoord(sizeGC))
  call gr_extendedGetCellCoords(axis, blockHandle, gr_meshMe, CENTER, .true., cellCenterCoord, sizeGC)
  ! assume this is uniform
  deltax = cellCenterCoord(2)-cellCenterCoord(1)
  ! done with this
  deallocate(cellCenterCoord)

  ! allocate temparary arrays once per call
  ! allocate arrays to hold extra eos data need by HSE stepper
  allocate(region_abar(sizeGC))
  allocate(region_zbar(sizeGC))
  allocate(region_ye(sizeGC))
  allocate(region_sumy(sizeGC))
  allocate(scalars(UNK_VARS_END-SPECIES_BEGIN+1,sizeGC))

  do k = 1, regionSize(THIRD_DIR)
     do j = 1, regionSize(SECOND_DIR)

        !-------------------
        !  first take care of the velocities
        if (face==HIGH) then
           do i = 1,guard
              ! zero-gradient everything (cannot know what the user has defined that needs to be propagated)
              regionData(regionSize(BC_DIR)-guard+i,j,k,:) = regionData(regionSize(BC_DIR)-guard,j,k,:)
              if (bcType == REFLECTING) then
                 regionData(regionSize(BC_DIR)-guard+i,j,k,VELX_VAR) = -regionData(regionSize(BC_DIR)-guard+1-i,j,k,VELX_VAR)
                 regionData(regionSize(BC_DIR)-guard+i,j,k,VELY_VAR) = regionData(regionSize(BC_DIR)-guard+1-i,j,k,VELY_VAR)
                 regionData(regionSize(BC_DIR)-guard+i,j,k,VELZ_VAR) = regionData(regionSize(BC_DIR)-guard+1-i,j,k,VELZ_VAR)
              else if (bcType == OUTFLOW) then
                 regionData(regionSize(BC_DIR)-guard+i,j,k,VELX_VAR) = regionData(regionSize(BC_DIR)-guard,j,k,VELX_VAR)
                 regionData(regionSize(BC_DIR)-guard+i,j,k,VELY_VAR) = regionData(regionSize(BC_DIR)-guard,j,k,VELY_VAR)
                 regionData(regionSize(BC_DIR)-guard+i,j,k,VELZ_VAR) = regionData(regionSize(BC_DIR)-guard,j,k,VELZ_VAR)
              else if (bcType == DIODE) then
                 regionData(regionSize(BC_DIR)-guard+i,j,k,VELX_VAR) = max(0.0,regionData(regionSize(BC_DIR)-guard,j,k,VELX_VAR))
                 regionData(regionSize(BC_DIR)-guard+i,j,k,VELY_VAR) = regionData(regionSize(BC_DIR)-guard,j,k,VELY_VAR)
                 regionData(regionSize(BC_DIR)-guard+i,j,k,VELZ_VAR) = regionData(regionSize(BC_DIR)-guard,j,k,VELZ_VAR)
              endif
          enddo
        else if (face==LOW) then
           do i = 1,guard
              ! zero-gradient everything (cannot know what the user has defined that needs to be propagated)
              regionData(i,j,k,:)     = regionData(guard+1,j,k,:)
              if (bcType == REFLECTING) then
                 regionData(i,j,k,VELX_VAR) = -regionData(2*guard+1-i,j,k,VELX_VAR)
                 regionData(i,j,k,VELY_VAR) =  regionData(2*guard+1-i,j,k,VELY_VAR)
                 regionData(i,j,k,VELZ_VAR) =  regionData(2*guard+1-i,j,k,VELZ_VAR)
              else if (bcType == OUTFLOW) then
                 regionData(i,j,k,VELX_VAR) = regionData(guard+1,j,k,VELX_VAR)
                 regionData(i,j,k,VELY_VAR) = regionData(guard+1,j,k,VELY_VAR)
                 regionData(i,j,k,VELZ_VAR) = regionData(guard+1,j,k,VELZ_VAR)
              else if (bcType == DIODE) then
                 regionData(i,j,k,VELX_VAR) = min(0.0,regionData(guard+1,j,k,VELX_VAR))
                 regionData(i,j,k,VELY_VAR) = regionData(guard+1,j,k,VELY_VAR)
                 regionData(i,j,k,VELZ_VAR) = regionData(guard+1,j,k,VELZ_VAR)
              endif
          enddo
        endif
       
        if (face==HIGH) then
           ! fill stuff out
           start = regionSize(BC_DIR)-guard+1
           end   = regionSize(BC_DIR)
           step  = 1
           direction = HSE_FORWARD
           lowi = start-2
           hii  = end
        else
           start = guard
           end   = 1
           step  = -1
           direction = HSE_BACKWARD
           lowi = end
           hii  = start+2
        endif

        ! do HSE calculation if we need a variable that it determines
        if (mask(DENS_VAR) .or. mask(TEMP_VAR) .or. mask(PRES_VAR) &
             .or. mask(EINT_VAR) .or. mask(GAME_VAR) .or. mask(GAMC_VAR) &
             .or. mask(ENER_VAR)) then
           ! we have to reformat array of scalars because the main data arrays
           ! which are normally fed to eos_externalAbarZbar() are in a different
           ! format than the regiodData passed into this subroutine
           ! ( the spatial indices come first in the main data arrays)
           do i = lowi, hii
              scalars(:,i) = regionData(i,j,k,SPECIES_BEGIN:UNK_VARS_END)
           enddo
           call eos_externalComputeAbarZbar(scalars(:,lowi:hii),region_abar(lowi:hii),region_zbar(lowi:hii))
           do i = lowi, hii
              region_sumy(i) = 1.0e0/region_abar(i)
              region_ye(i) = region_sumy(i)*region_zbar(i)
           enddo
           do i = start, end, step
              ! density guess value
              regionData(i,j,k,DENS_VAR) = regionData(i-step,j,k,DENS_VAR)
              call sim_hse_step(regionData(:,j,k,DENS_VAR), &
                                regionData(:,j,k,TEMP_VAR), &
                                region_ye,   &
                                region_sumy, &
                                i,sim_grav, deltax, direction,2, HSE_SETTEMP)
      
              ! now get all the eos stuff and fill it in
              eosData(EOS_DENS) = regionData(i,j,k,DENS_VAR)
              eosData(EOS_TEMP) = regionData(i,j,k,TEMP_VAR)
              eosData(EOS_ABAR) = region_abar(i)
              eosData(EOS_ZBAR) = region_zbar(i)
              call Eos(MODE_DENS_TEMP, 1, eosData)
              regionData(i,j,k,PRES_VAR) = eosData(EOS_PRES)
              regionData(i,j,k,EINT_VAR) = eosData(EOS_EINT)
              regionData(i,j,k,GAME_VAR) = eosData(EOS_PRES)/(eosData(EOS_EINT)*eosData(EOS_DENS)) +1.0
              regionData(i,j,k,GAMC_VAR) = eosData(EOS_GAMC)
              regionData(i,j,k,ENER_VAR) = eosData(EOS_EINT) + regionData(i,j,k,VELX_VAR)**2 &
                                         + regionData(i,j,k,VELY_VAR)**2 + regionData(i,j,k,VELZ_VAR)**2
           enddo
        endif

     enddo ! second coord
  enddo ! third coord

  deallocate(region_abar)
  deallocate(region_zbar)
  deallocate(region_ye)
  deallocate(region_sumy)
  deallocate(scalars)

  return
end subroutine Grid_bcApplyToRegionSpecialized
