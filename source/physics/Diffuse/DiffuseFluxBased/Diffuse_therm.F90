!!****if* source/physics/Diffuse/DiffuseFluxBased/Diffuse_therm
!!
!! NAME
!!
!!  Diffuse_therm
!!
!!
!! SYNOPSIS
!!
!! call Diffuse_therm(integer(IN) :: sweepDir, 
!!                    integer(IN) :: igeom, 
!!                    integer(IN) :: blockID,
!!                    integer(IN) :: numCells,
!!                    integer(IN) :: blkLimits(2,MDIM),
!!                    integer(IN) :: blkLimitsGC(2,MDIM),
!!                    real(IN)    :: leftCoords(MAXCELLS or numCells),
!!                    real(IN)    :: rightCoords(MAXCELLS or numCells),
!!                    real(INOUT) :: temp_flx(NFLUXES,:,:,:),
!!                    real(IN)    :: areaLeft(:,:,:))
!!
!!
!! DESCRIPTION
!!
!!  Diffuse_therm alters the energy flux to account for heat losses through
!!  thermal diffusion.  This is an explicit method, so a timestep limiter
!!  will be required.  Stability is guaranteed for 
!! 
!!                           dx**2 
!!                 dt < .5* -------
!!                             D
!!
!!  where D is the diffusion coefficient, related to the isochoric conductivity (sigma)
!!  and specific heat at constant volume (c_v) by
!!
!!               D = sigma/(rho*c_v)
!!
!!  Fluxes in PPM are stored at the zone boundaries, the temperature, etc. are
!!  at the zone centers.  So the flux {F = -sigma * grad(T)} is
!! 
!!                      sigma_i + sigma_{i-1}     T_i - T_{i-1}
!!             F_i = - ----------------------- * --------------- 
!!                               2                      dx
!!
!!  and then the contribution to the energy equation (evolution of rho*E),
!!  where E is the energy/gram is - div(F) or
!! 
!!                        F_{i+1} - F_i
!!             dE_i  = - ---------------      
!!                             dx    
!!
!! 
!!  together, these give
!! 
!!                     sigma_{i+1} + sigma_i     T_{i+1} - T_i
!!            dE_i =  ----------------------- * ---------------   -
!!                               2                   dx**2
!!
!!                     sigma_i + sigma_{i-1}     T_i - T_{i-1}
!!                    ----------------------- * ---------------   -
!!                               2                   dx**2
!!
!!  which reduces to the standard expression for diffusion if sigma is constant.
!!                     
!!
!! NOTES 
!!
!!    This routine is used by the PPM implementation of the Hydro unit.
!!    Other implementations of Hydro, in particular MHD and RHD, may have their
!!    own mechanisms for handling diffusive effects that bypass Diffuse_therm
!!    (and the Diffuse code unit in general), or they may lack support for
!!    diffusive effects.
!!
!!    This routine computes the heat fluxes and adds them to the energy flux
!!    returned from hydro_1d.  It is to be called after calling hydro_1d on a 
!!    block.  The updated energy fluxes are then used in hy_updateSoln to 
!!    produce the updated energy.
!!
!!    Boundary conditions on the temperature are handled through the GridBoundaryConditions
!!    subunit, not in this routine.  Ex: setting a reflecting boundary will reflect the
!!    temperature into the guard cells, which are used in computing the flux
!!    at the block boundaries.
!!
!!
!!
!! ARGUMENTS
!!
!!  sweepDir    -   the current sweep direction
!!
!!  igeom       -   the geometry flag for the current sweep direction
!!
!!  blockID     -   the block number to operate on
!!  numCells    -   the number of Cells along the sweep direction
!!  leftCoords  -   Coordinates of the left edge of the zones
!!  rightCoords -   Coordinates of the right edge of the zones
!!
!!  temp_flx    -   Temporary storage for fluxes in sweep direction
!!  areaLeft    -   Cell face areas at the left (smaller) cell side
!!  blkLimits -   endpoints of block indices without including gcells
!!  blkLimitsGC -   endpoints of block indices including gcells
!!
!! HISTORY
!!
!!  Apparently this subroutine used to be called therm_explicit.
!!***

#ifdef DEBUG_ALL
#define DEBUG_DIFFUSE
#endif
subroutine Diffuse_therm(sweepDir, igeom, blockID,numCells,blkLimits,blkLimitsGC,&
                           leftCoords,rightCoords,&
                           temp_flx, areaLeft)


  use Diffuse_data, ONLY : diffusion_cutoff_density, useDiffuse, useDiffuseTherm,&
       thermal_diff_method, diff_geometricMeanDiff, diff_geometry, &
       diff_scaleFactThermFlux
  use Conductivity_interface, ONLY: Conductivity
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  
  implicit none
#include "constants.h"  
#include "Flash.h"
  
  integer, intent(IN) :: sweepDir, igeom, blockID, numCells
  integer, dimension(2,MDIM), intent (IN) :: blkLimitsGC, blkLimits
  
  
  
#ifdef FIXEDBLOCKSIZE
  real, intent(INOUT), DIMENSION(NFLUXES,                   &
                               GRID_ILO_GC:GRID_IHI_GC,     &
                               GRID_JLO_GC:GRID_JHI_GC,     &
                               GRID_KLO_GC:GRID_KHI_GC) ::  &
                               temp_flx
  real, intent(IN), DIMENSION(                   &
                               GRID_ILO_GC:GRID_IHI_GC,     &
                               GRID_JLO_GC:GRID_JHI_GC,     &
                               GRID_KLO_GC:GRID_KHI_GC) ::  &
                               areaLeft
  real,intent(IN), DIMENSION(MAXCELLS) :: leftCoords ,rightCoords
  real,dimension(MAXCELLS):: cond, condi, thermflux,dx, dxinv
#else
  real, intent(INOUT), DIMENSION(NFLUXES,                   &
                               blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                               temp_flx
  real, intent(IN), DIMENSION(                   &
                               blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                               areaLeft
  
  real, intent(IN), DIMENSION(numCells) :: leftCoords ,rightCoords
  real,dimension(numCells) :: cond, condi, thermflux,dx, dxinv
#endif

!! Do not use this implementation for the MHD and Hydro unsplit solvers.
#ifndef FLASH_HYDRO_UNSPLIT
  real, pointer, DIMENSION(:,:,:,:) :: solnData
  

! storage for the opacities and fluxes

  real :: diff_coeff


! storage for the 1d mass fractions in a zone
  real :: massfrac(NSPECIES), cond_zone, xtemp, xdens

! loop indices
  integer :: i, j, k, n

  integer :: tempVar
  logical :: supported, temp

  real :: dx41, dx40, dx31, dx30, dx21, dx20


  if(.not.useDiffuse) return
  if(.not.useDiffuseTherm) return !! simple return if turned off

#ifdef DEBUG_DIFFUSE

  if((sweepDir==SWEEP_Y).and.(NDIM<2))&
       call Driver_abortFlash("Diffuse_therm : wrong dimensionality for SWEEP_Y")
  if((sweepDir==SWEEP_Z).and.(NDIM<3))&
       call Driver_abortFlash("Diffuse_therm : wrong dimensionality for SWEEP_Z")

#endif

!  useFaceAreas = diff_useCellAreasForFluxes

!-----------------------------------------------------------------------------
! x-sweep
!-----------------------------------------------------------------------------


  call Grid_getBlkPtr(blockID,solnData)
  if (sweepDir .EQ. SWEEP_X) then

! compute dx for this block

     do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1
        dx(i) = rightCoords(i) - leftCoords(i)

        if (igeom == XYZ) then  
           dxinv(i) = 1.e0/dx(i)
        else if (igeom ==RAD_CYL) then 
           dx31 = (rightCoords(i)**3 - leftCoords(i)**3)/3.
           dx30 = (rightCoords(i-1)**3 - leftCoords(i-1)**3)/3.
           dx21 = (rightCoords(i)**2 - leftCoords(i)**2)/2.
           dx20 = (rightCoords(i-1)**2 - leftCoords(i-1)**2)/2.

           dxinv(i) = 1.e0/(dx31/dx21 - dx30/dx20)
        else if (igeom == RAD_SPH) then 
           dx41 = (rightCoords(i)**4 - leftCoords(i)**4)/4.
           dx40 = (rightCoords(i-1)**4 - leftCoords(i-1)**4)/4.
           dx31 = (rightCoords(i)**3 - leftCoords(i)**3)/3.
           dx30 = (rightCoords(i-1)**3 - leftCoords(i-1)**3)/3.

           dxinv(i) = 1.e0/(dx41/dx31 - dx40/dx30)
        endif

     enddo

     tempVar = TEMP_VAR
#ifdef E3_VAR
#ifdef ERAD_VAR
#ifdef TRAD_VAR
     tempVar = TRAD_VAR
#endif
#endif
#endif

     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           
! compute the isochoric conductivity -- we need the conductivities in the guard cells 
! immediately adjacent to the block we are operating on
           do i = blkLimits(LOW,IAXIS)-1, blkLimits(HIGH,IAXIS)+1

              xtemp = solnData(tempVar,i,j,k)
              xdens = solnData(DENS_VAR,i,j,k)

! load the mass fractions
              do n = 1, NSPECIES
                 massfrac(n) = solnData(SPECIES_BEGIN-1+n,i,j,k)
              enddo
              call Conductivity(xtemp, xdens, massfrac, cond_zone, diff_coeff, 2)
              
              cond(i) = cond_zone
           enddo

           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1
              if (diff_geometricMeanDiff) then
                 condi(i) = 2.*(cond(i)*cond(i-1))/(cond(i)+cond(i-1))
              else
                 condi(i) = 0.5*(cond(i)+cond(i-1))
              endif
           enddo
           
! compute the fluxes
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1
              
              if (thermal_diff_method == 1) then
                 thermflux(i) = - condi(i) * & 
                      dxinv(i)*(solnData(tempVar,i  ,j,k) -  & 
                      solnData(tempVar,i-1,j,k))

              else
                 
                 thermflux(i) = -0.5*(0.25*cond(i)/solnData(tempVar,i,j,k)**3 + &
                                      0.25*cond(i-1)/solnData(tempVar,i-1,j,k)**3) * &
                                      dxinv(i)*(solnData(tempVar,i,j,k)**4 - &
                                                solnData(tempVar,i-1,j,k)**4)
              endif

! if the density is too low, turn off thermal diffusion
              if (solnData(DENS_VAR,i,j,k) < diffusion_cutoff_density) then
                 thermflux(i) = 0.0
              else
                 thermflux(i) = thermflux(i) * diff_scaleFactThermFlux
              end if

! update the temporary fluxes from hydro to include the thermal flux,
! these fluxes are then passed on to update_soln
!
! FLASH 1.6 change -- the internal energy flux also gets the heat flux
!

#ifdef ERAD_FLUX
              temp_flx(ERAD_FLUX,i,j,k) = temp_flx(ERAD_FLUX,i,j,k) + areaLeft(i,j,k)*thermflux(i)
#endif

#ifdef E3_FLUX
              temp_flx(E3_FLUX,i,j,k) = &
                   temp_flx(E3_FLUX,i,j,k) + areaLeft(i,j,k)*thermflux(i)
#endif

              temp_flx(E_FLUX,i,j,k) = temp_flx(E_FLUX,i,j,k) + areaLeft(i,j,k)*thermflux(i)

              temp_flx(EINT_FLUX,i,j,k) = &
                   temp_flx(EINT_FLUX,i,j,k) + areaLeft(i,j,k)*thermflux(i)

           enddo

        enddo
     enddo

!-----------------------------------------------------------------------------
! y-sweep
!-----------------------------------------------------------------------------
  elseif (sweepDir .EQ. SWEEP_Y )then

! compute dx for this block -- in the y direction
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
        dx(j) = rightCoords(j) - leftCoords(j)
        dxinv(j) = 1.e0/dx(j)
     enddo

     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

! compute the isochoric conductivity -- we need the conductivities in the guard cells 
! immediately adjacent to the block we are operating on
              do j = blkLimits(LOW,JAXIS)-1,blkLimits(HIGH,JAXIS)+1

              xtemp = solnData(TEMP_VAR,i,j,k)
              xdens = solnData(DENS_VAR,i,j,k)

! load the mass fractions
              do n = 1, NSPECIES
                 massfrac(n) = solnData(SPECIES_BEGIN-1+n,i,j,k)
              enddo

              call Conductivity(xtemp, xdens, massfrac, cond_zone, diff_coeff, 2)

              cond(j) = cond_zone

           enddo


           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
              if (diff_geometricMeanDiff) then
                condi(j) = 2.*(cond(j)*cond(j-1))/(cond(j)+cond(j-1))
              else
                condi(j) = 0.5*(cond(j)+cond(j-1))
              endif
           enddo
               
! compute the fluxes
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1

              if (thermal_diff_method == 1) then
                 thermflux(j) = - condi(j) * & 
                      dxinv(j)*(solnData(TEMP_VAR,i,j  ,k) -  & 
                                solnData(TEMP_VAR,i,j-1,k))
              else
                 thermflux(j) = -0.5*(0.25*cond(j)/solnData(TEMP_VAR,i,j,  k)**3 + &
                                      0.25*cond(j-1)/solnData(TEMP_VAR,i,j-1,k)**3) * &
                                      dxinv(j)*(solnData(TEMP_VAR,i,j  ,k)**4 - &
                                                solnData(TEMP_VAR,i,j-1,k)**4)
              endif


! if the density is too low, turn off thermal diffusion

              if (solnData(DENS_VAR,i,j,k) < diffusion_cutoff_density) then
                 thermflux(j) = 0.0
              else
                 thermflux(j) = thermflux(j) * diff_scaleFactThermFlux
              end if


! update the temporary fluxes from hydro to include the thermal flux,
! these fluxes are then passed on to update_soln
              temp_flx(E_FLUX,i,j,k) = temp_flx(E_FLUX,i,j,k) + areaLeft(i,j,k)*thermflux(j)

              temp_flx(EINT_FLUX,i,j,k) = &
                   temp_flx(EINT_FLUX,i,j,k) + areaLeft(i,j,k)*thermflux(j)
           enddo
                  
        enddo
     enddo

!-----------------------------------------------------------------------------
! z-sweep
!-----------------------------------------------------------------------------
  elseif (sweepDir .EQ. SWEEP_Z ) then


! compute dx for this block -- in the z direction
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
        dx(k) = rightCoords(k) - leftCoords(k)
        dxinv(k) = 1.e0/dx(k)
     enddo

     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

! compute the isochoric conductivity -- we need the conductivities in the guard cells 
! immediately adjacent to the block we are operating on
           do k = blkLimits(LOW,KAXIS)-1,blkLimits(HIGH,KAXIS)+1

              xtemp = solnData(TEMP_VAR,i,j,k)
              xdens = solnData(DENS_VAR,i,j,k)

! load the mass fractions
              do n = 1, NSPECIES
                 massfrac(n) = solnData(SPECIES_BEGIN-1+n,i,j,k)
              enddo

              call Conductivity(xtemp, xdens, massfrac, cond_zone, diff_coeff, 2)

              cond(k) = cond_zone
              
           enddo

           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
              if (diff_geometricMeanDiff) then
                condi(k) = 2.*(cond(k)*cond(k-1))/(cond(k)+cond(k-1))
              else
                condi(k) = 0.5*(cond(k)+cond(k-1))
              endif
           enddo
               
! compute the fluxes
           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1

              if (thermal_diff_method == 1) then
                 thermflux(k) = - condi(k) * & 
                      dxinv(k)*(solnData(TEMP_VAR,i,j,k  ) -  & 
                                solnData(TEMP_VAR,i,j,k-1))
              else
                 thermflux(k) = -0.5*(0.25*cond(k)/solnData(TEMP_VAR,i,j,k  )**3 + &
                                      0.25*cond(k-1)/solnData(TEMP_VAR,i,j,k-1)**3) * &
                                      dxinv(k)*(solnData(TEMP_VAR,i,j,k  )**4 - &
                                                solnData(TEMP_VAR,i,j,k-1)**4)
              endif

! if the density is too low, turn off thermal diffusion

              if (solnData(DENS_VAR,i,j,k) < diffusion_cutoff_density) then
                 thermflux(k) = 0.0
              else
                 thermflux(k) = thermflux(k) * diff_scaleFactThermFlux
              end if

! update the temporary fluxes from hydro to include the thermal flux,
! these fluxes are then passed on to update_soln
              temp_flx(E_FLUX,i,j,k) = temp_flx(E_FLUX,i,j,k) + areaLeft(i,j,k)*thermflux(k)

              temp_flx(EINT_FLUX,i,j,k) = &
                   temp_flx(EINT_FLUX,i,j,k) + areaLeft(i,j,k)*thermflux(k)
           enddo
                  
        enddo
     enddo
         
  endif
  call Grid_releaseBlkPtr(blockID,solnData)
#endif
!! End of #ifndef FLASH_HYDRO_UNSPLIT
  return
end subroutine Diffuse_therm

      

