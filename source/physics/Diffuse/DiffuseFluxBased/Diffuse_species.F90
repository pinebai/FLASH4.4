!!****if* source/physics/Diffuse/DiffuseFluxBased/Diffuse_species
!!
!! NAME
!!
!!  Diffuse_species
!!
!!
!! SYNOPSIS
!!
!! Diffuse_species(integer(IN) :sweepDir, 
!!                 integer(IN) : igeom, 
!!                 integer(IN) : blockID,
!!                 integer(IN) : numCells,
!!                 integer(:,:)(IN):blkLimits, 
!!                 integer(:,:)(IN):blkLimitsGC
!!                 real(:)(IN):leftCoords,
!!                 real(:)(IN):rightCoords,
!!                 real(:,:)(INOUT) :temp_flx, 
!!                 real(:,:)(INOUT) :temp_fly, 
!!                 real(:,:)(INOUT) :temp_flz)
!!
!!
!! DESCRIPTION
!!
!!  Diffuse_species alters the mass-fraction flux to account for 
!!  material diffusion.  This is an explicit method, so a timestep limiter
!!  will be required.  Stability is guaranteed for
!!
!!                           dx**2
!!                 dt < .5* -------
!!                             D
!!
!!  where D is the diffusion coefficient.
!!
!!  Fluxes in PPM are stored at the zone boundaries; mass fractions are 
!!  at the zone centers.  So the flux {F = -sigma * grad(T)} is
!!
!!                      D_i + D_{i-1}     T_i - T_{i-1}
!!             F_i = - -------------- * ---------------
!!                           2                  dx
!!
!! NOTES
!!
!!    This routine computes the species fluxes and adds them to the flux
!!    returned from hydro.  It is to be called after calling hydro_1d on a
!!    block.  The updated energy fluxes are then used in update_soln to
!!    produce the updated energy.
!!
!!    Boundary conditions are handled through tot_bnd,
!!    not in this routine.  Ex: setting a reflecting boundary will reflect the
!!    species into the guard cells, which are used in computing the flux
!!    at the block boundaries.
!!
!!
!! ARGUMENTS
!!
!!  sweepDir        the current sweep direction
!!
!!  igeom           the geometry flag for the current sweep direction
!!
!!  blockID        the block number to operate on
!!  numCells        the number of Cells along the sweep direction
!!  leftCoords      Coordinates of the left edge of the zones
!!  rightCoords     Coordinates of the right edge of the zones
!!
!!  tempflx         Temporary storage for flux along first direction
!!  tempfly         Temporary storage for flux along second direction
!!  tempflz         Temporary storage for flux along third direction
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_DIFFUSE
#endif
subroutine Diffuse_species(sweepDir, igeom, blockID,numCells,blkLimits,blkLimitsGC,&
                           leftCoords,rightCoords,&
                           temp_flx, temp_fly, temp_flz)


  use Diffuse_data, ONLY : useDiffuse, useDiffuseSpecies, diff_geometricMeanDiff,&
                           diffusion_cutoff_density
  use MassDiffusivity_interface, ONLY : MassDiffusivity
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
    
  
  implicit none
#include "constants.h"  
#include "Flash.h"

  
  
  integer, intent(IN) :: sweepDir, igeom, blockID, numCells
  integer, dimension(2,MDIM), intent (IN) :: blkLimitsGC, blkLimits
  

#ifdef FIXEDBLOCKSIZE
  real, intent(OUT), DIMENSION(NFLUXES,                   &
                               GRID_ILO_GC:GRID_IHI_GC,     &
                               GRID_JLO_GC:GRID_JHI_GC,     &
                               GRID_KLO_GC:GRID_KHI_GC) ::  &
                               temp_flx, temp_fly, temp_flz
  real,intent(IN), DIMENSION(MAXCELLS) :: leftCoords ,rightCoords
  real,dimension(MAXCELLS) :: cond, condi, thermflux,dx, dxinv
  real, dimension(MAXCELLS) :: temporary, &
                               rho, mass_diff, mass_diffi
  real, dimension(MAXCELLS,NSPECIES) :: xn
#else
  real, intent(OUT), DIMENSION(NFLUXES,                   &
                               blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                               temp_flx, temp_fly, temp_flz
  
  real, intent(IN), DIMENSION(numCells) :: leftCoords ,rightCoords
  real,dimension(numCells):: cond, condi, thermflux,dx, dxinv
  real, dimension(numCells) :: temporary, &
                               rho, mass_diff, mass_diffi
  real, dimension(numCells,NSPECIES) :: xn
#endif

  real, pointer, DIMENSION(:,:,:,:) :: solnData

  real                  :: dx41, dx40, dx31, dx30, dx21, dx20
    
  integer               :: i, j, k, ii, n
  real                  :: xtemp, xdens, mass_diffusivity_zone
  real, dimension(NSPECIES)   :: abar, massfrac, mflux
  real                  :: D
  logical :: supported,temp




! storage for the opacities and fluxes

  real :: diff_coeff

! storage for the 1d mass fractions in a zone
  real :: cond_zone



     
  if(.not.useDiffuse) return
  if(.not.useDiffuseSpecies) return  !! If diffuse species is not wanted just return

#ifdef DEBUG_DIFFUSE
  
  
  if(.NOT.useDiffuseSpecies) call Driver_abortFlash("Diffuse_species, in wrong place")
  if((sweepDir==SWEEP_Y).and.(NDIM<2))&
       call Driver_abortFlash("species : wrong dimensionality")
  if((sweepDir==SWEEP_Z).and.(NDIM<3))&
       call Driver_abortFlash("species : wrong dimensionality")
  
#endif
  


  call Grid_getBlkPtr(blockID,solnData)

  !-----------------------------------------------------------------------------
  ! x-sweep
  !-----------------------------------------------------------------------------
  if (sweepDir == SWEEP_X)then
     ! compute dx for this block
          
     do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1
        dx(i) = rightCoords(i) - leftCoords(i)
        
        if (igeom == XYZ) then 
           dxinv(i) = 1.e0/dx(i)
        else if (igeom == RAD_CYL) then 
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
     
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           
           rho(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)) = &
                solnData(DENS_VAR, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), j, k)
           
           do ii= 1, NSPECIES
              xn(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),ii) = &
                   solnData(SPECIES_BEGIN + ii -1, &
                   blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), j, k)
              
           enddo
           
           ! compute the mass diffusivity -- we need the mass diffusivity
           ! in the guard cells immediately adjacent to the block 
           ! we are operating on
           do i = blkLimits(LOW,IAXIS)-1, blkLimits(HIGH,IAXIS)+1
              ! data doens't point to anthing: KMO
              !              xtemp = data(TEMP_VAR, i, j, k)
              !              xdens = data(DENS_VAR, i, j, k)
              xtemp = solnData(TEMP_VAR, i, j, k)
              xdens = solnData(DENS_VAR, i, j, k)
              
              ! load the mass fractions
              do n = 1, NSPECIES
                 massfrac(n) = solnData(SPECIES_BEGIN-1+n, i, j, k)
              enddo
              
              call MassDiffusivity(xtemp,xdens,massfrac,mass_diffusivity_zone)
              
              xn(i,:) = xn(i,:)
              
              mass_diff(i) = mass_diffusivity_zone*xdens
           enddo
           
           
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1
              if (diff_geometricMeanDiff) then
                 mass_diffi(i) = 2.*(mass_diff(i)*mass_diff(i-1))/&
                      (mass_diff(i)+mass_diff(i-1))
              else
                 mass_diffi(i) = 0.5*(mass_diff(i)+mass_diff(i-1))
              endif
           enddo
           
           ! compute the fluxes
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1
              
              do ii = 1, NSPECIES
                 mflux(ii) = - mass_diffi(i) * dxinv(i)*( xn(i,ii)  - xn(i-1,ii) )
                 
                 if (igeom == XYZ) then

                    temp_flx(SPECIES_FLUX_BEGIN+ii-1,i,j,k)   =  &
                         temp_flx(SPECIES_FLUX_BEGIN+ii-1,i,j,k) &
                         + mflux(ii)
                 elseif (igeom == RAD_CYL) then
                    temp_flx(SPECIES_FLUX_BEGIN+ii-1,i,j,k)   =  &
                         temp_flx(SPECIES_FLUX_BEGIN+ii-1,i,j,k) &
                         + mflux(ii)*abs(leftCoords(i))
                 elseif (igeom == RAD_SPH) then
                    temp_flx(SPECIES_FLUX_BEGIN+ii-1,i,j,k)   =  & 
                         temp_flx(SPECIES_FLUX_BEGIN+ii-1,i,j,k) &
                         + mflux(ii)*leftCoords(i)*leftCoords(i)
                 endif
              enddo
           enddo
        enddo
     enddo
     
     !-----------------------------------------------------------------------------
     ! y-sweep
     !-----------------------------------------------------------------------------
  elseif (sweepDir .EQ. SWEEP_Y) then
     
     dx = rightCoords - leftCoords
     dxinv = 1.e0/dx
     
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           rho(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)) = &
                solnData(DENS_VAR, i, blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), k)

           do ii= 1, NSPECIES
              xn(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),ii) = &
                   solnData(SPECIES_BEGIN + ii -1, i, &
                   blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), k)
              
           enddo
           
           ! compute the mass diffusivity -- we need the mass diffusivity
           ! in the guard cells immediately adjacent to the block 
           ! we are operating on
           do j = blkLimits(LOW,JAXIS)-1,blkLimits(HIGH,JAXIS)+1
              xtemp = solnData(TEMP_VAR, i, j, k)
              xdens = solnData(DENS_VAR, i, j, k)
              
              ! load the mass fractions
              do n = 1, NSPECIES
                 massfrac(n) = solnData(SPECIES_BEGIN-1+n, i, j, k)
              enddo
              
              call MassDiffusivity(xtemp,xdens,massfrac,mass_diffusivity_zone)
              
              xn(j,:) = xn(j,:)
              
              mass_diff(j) = mass_diffusivity_zone*xdens
           enddo
           
           
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
              if (diff_geometricMeanDiff) then
                 mass_diffi(j) = 2.*(mass_diff(j)*mass_diff(j-1))/&
                      (mass_diff(j)+mass_diff(j-1))
              else
                 mass_diffi(j) = 0.5*(mass_diff(j)+mass_diff(j-1))
              endif
           enddo
           
           ! compute the fluxes
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
              
              do ii = 1, NSPECIES
                 mflux(ii) = - mass_diffi(j) * dxinv(j)*( xn(j,ii)  - xn(j-1,ii) )
                 if (igeom == 0) then
                    temp_fly(SPECIES_FLUX_BEGIN+ii-1,i,j,k)   =  &
                         temp_fly(SPECIES_FLUX_BEGIN+ii-1,i,j,k)  &
                         + mflux(ii) 
                 elseif (igeom == 1) then
                    temp_fly(SPECIES_FLUX_BEGIN+ii-1,i,j,k)   =  &
                         temp_fly(SPECIES_FLUX_BEGIN+ii-1,i,j,k)  &
                         + abs(leftCoords(j))*mflux(ii) 
                 elseif (igeom == 2) then
                    temp_fly(SPECIES_FLUX_BEGIN+ii-1,i,j,k)   =  &
                         temp_fly(SPECIES_FLUX_BEGIN+ii-1,i,j,k)  &
                         + leftCoords(j)*leftCoords(j)*mflux(ii) 
                 endif
              enddo
              
           enddo
        enddo
     enddo
     
     !-----------------------------------------------------------------------------
     ! z-sweep
     !-----------------------------------------------------------------------------
  elseif (sweepDir .EQ. SWEEP_Z) then
     
     
     dx = rightCoords - leftCoords
     dxinv = 1.e0/dx
     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           rho(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) = &
                solnData(DENS_VAR, i, j, blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
           
           do ii= 1, NSPECIES
              xn(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS),ii) = &
                   solnData(SPECIES_BEGIN + ii -1, i, j, &
                   blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
              
           enddo
           ! compute the mass diffusivity -- we need the mass diffusivity
           ! in the guard cells immediately adjacent to the block 
           ! we are operating on
           do k = blkLimits(LOW,KAXIS)-1,blkLimits(HIGH,KAXIS)+1
              xtemp = solnData(TEMP_VAR, i, j, k)
              xdens = solnData(DENS_VAR, i, j, k)
              
              ! load the mass fractions
              do n = 1, NSPECIES
                 massfrac(n) = solnData(SPECIES_BEGIN-1+n, i, j, k)
              enddo
              
              call MassDiffusivity(xtemp,xdens,massfrac,mass_diffusivity_zone)
              
              xn(k,:) = xn(k,:)
              
              mass_diff(k) = mass_diffusivity_zone*xdens
           enddo
           
           
           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
              if (diff_geometricMeanDiff) then
                 mass_diffi(k) = 2.*(mass_diff(k)*mass_diff(k-1))/&
                      (mass_diff(k)+mass_diff(k-1))
              else
                 mass_diffi(k) = 0.5*(mass_diff(k)+mass_diff(k-1))
              endif
           enddo
           
           ! compute the fluxes
           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
              
              do ii = 1, NSPECIES
                 mflux(ii) = -mass_diffi(k)* dxinv(k)*&
                      ( xn(k,ii)  - xn(k-1,ii) )
                 if (igeom == 0) then
                    temp_flz(SPECIES_FLUX_BEGIN+ii-1,i,j,k)   =  &
                         temp_fly(SPECIES_FLUX_BEGIN+ii-1,i,j,k)  &
                         + mflux(ii)
                 elseif (igeom == 1) then
                    temp_flz(SPECIES_FLUX_BEGIN+ii-1,i,j,k)   =  &
                         temp_fly(SPECIES_FLUX_BEGIN+ii-1,i,j,k)  &
                         + abs(leftCoords(j))*mflux(ii)
                 elseif (igeom == 2) then
                    temp_flz(SPECIES_FLUX_BEGIN+ii-1,i,j,k)   =  &
                         temp_fly(SPECIES_FLUX_BEGIN+ii-1,i,j,k)  &
                         + leftCoords(j)*leftCoords(j)*mflux(ii)
                 endif
              enddo
              
           enddo
        enddo
     enddo
     
  endif
  call Grid_releaseBlkPtr(blockID,solnData)
  
  return
end subroutine Diffuse_species
