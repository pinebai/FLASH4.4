!!****if* source/physics/Diffuse/DiffuseFluxBased/Diffuse_visc
!!
!! NAME
!!
!!  Diffuse_visc
!!
!!
!! SYNOPSIS
!!
!! call Diffuse_visc (integer(IN) :: sweepDir, 
!!                    integer(IN) :: igeom, 
!!                    integer(IN) :: blockID,
!!                    integer(IN) :: numCells,
!!                    integer(IN) :: blkLimits(2,MDIM),
!!                    integer(IN) :: blkLimitsGC(2,MDIM),
!!                    real(IN)    :: leftCoords(MAXCELLS or numCells),
!!                    real(IN)    :: rightCoords(MAXCELLS or numCells),
!!                    real(INOUT) :: temp_flx(NFLUXES,:,:,:),
!!                    real(IN)    :: areaLeft(:,:,:),
!!                    real(IN)    :: secondCoord(MAXCELLS or numCells),
!!                    real(IN)    :: thirdCoord(MAXCELLS or numCells))
!!
!!
!! DESCRIPTION
!!
!!
!! Diffuse_visc alters the velocity fluxes for a block to account for viscosity.
!!
!! Energy fluxes (for E_FLUX and, if defined, EINT_FLUX) are also adjusted
!! accordingly.
!!
!! This is an explicit method, so a timestep limiter will be required.  
!! Stability is guaranteed for 
!! 
!!                           dx**2 
!!                 dt < .5* -------
!!                            nu
!!
!! Fluxes in PPM are stored at the zone boundaries, the temperature, etc. are
!! at the zone centers.  So the velocity fluxes (more exactly, momentum fluxes)
!! acquire terms like  F = -nu * grad(v), which is approximated as
!! 
!!                       nu_i + nu_{i-1}      v(xyz)_i - v(xyz)_{i-1}
!!         F(xyz)_i = - -----------------  * ------------------------   .
!!                              2                      dx
!!
!!
!! This routine computes the velocity and energy fluxes from
!! viscosity and adds them to the corresponding fluxes returned from
!! hydro_1d.  It is to be called after calling hydro_1d on a block.
!! The updated energy fluxes are then used in hy_ppm_updateSoln to
!! produce the updated velocities and energy.
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
!!  temp_flx    -   Temporary storage for fluxes along sweep direction
!!  areaLeft    -   Cell face areas at the left (smaller) cell side
!!
!!  secondCoord -  for an x sweep: y coordinate; for a y sweep: x coord; for a z sweep: x coord
!!  thirdCoord  -  for an x sweep: z coordinate; for a y sweep: z coord; for a z sweep: y coord
!!  blkLimits -   endpoints of block indices without including gcells
!!  blkLimitsGC -   endpoints of block indices including gcells
!!
!!
!! NOTES
!!
!!    This routine is used by the PPM implementation of the Hydro unit.
!!    Other implementations of Hydro, in particular MHD and RHD, may have their
!!    own mechanisms for handling diffusive effects that bypass Diffuse_visc
!!    (and the Diffuse code unit in general), or they may lack support for
!!    diffusive effects.
!!
!! SEE ALSO
!!
!!  Diffuse_therm
!!
!! HISTORY
!!
!!  Apparently this subroutine used to be called visc_explicit.
!!  Visc_explicit started out as a very lightly modified version of therm_explicit.
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_DIFFUSE
#endif
subroutine Diffuse_visc(sweepDir, igeom, blockID,numCells,blkLimits,blkLimitsGC,&
                           leftCoords,rightCoords,&
                           temp_flx, areaLeft, secondCoord, thirdCoord)


  use Diffuse_data, ONLY : useDiffuse, useDiffuseVisc, diff_geometricMeanDiff
  use Viscosity_interface, ONLY : Viscosity
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
  real,intent(IN), DIMENSION(MAXCELLS) :: leftCoords ,rightCoords, &
       secondCoord, thirdCoord
  real, dimension(MAXCELLS) :: dx, dxinv, u, ut, utt, rho, visc, visci,viscflux
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
  real, intent(IN), DIMENSION(numCells) :: leftCoords ,rightCoords, &
       secondCoord, thirdCoord
  real, dimension(numCells) :: dx, dxinv, u, ut, utt, rho, visc, visci,viscflux
#endif

!! Do not use this implementation for the MHD and Hydro unsplit solvers.
#ifndef FLASH_HYDRO_UNSPLIT

  real, pointer, DIMENSION(:,:,:,:) :: solnData
  
  real :: diff_coeff
  integer :: i,j,k,n
  ! storage for the 1d mass fractions in a zone
  real :: massfrac(NSPECIES), xtemp, xdens, viscDynamic, viscUnusedHere
  real :: twodxInv, twodyInv, twodzInv, enerFlux
  
  real :: dx41, dx40, dx31, dx30, dx21, dx20




  ! check the geometry -- currently, only cartesian geometry is supported

  if(.not.useDiffuse) return
  if(.not.useDiffuseVisc) return !! simply return if turned off
  
#ifdef DEBUG_DIFFUSE
  
  ! check the geometry -- currently, only cartesian, 2-d cylindrical, and
  ! 1-d spherical geometry is supported
  

  if(.NOT.useDiffuseVisc) call Driver_abortFlash("explicit Diffuse_visc, in wrong place")
  if((sweepDir==SWEEP_Y).and.(NDIM<2))call Driver_abortFlash("Diffuse_visc: invalid sweep direction for dimensionality")
  if((sweepDir==SWEEP_Z).and.(NDIM<3))call Driver_abortFlash("Diffuse_visc: invalid sweep direction for dimensionality")
#endif

!  useFaceAreas = diff_useCellAreasForFluxes

!----------------------------------------------------------------------------
! x-sweep
!-----------------------------------------------------------------------------
  call Grid_getBlkPtr(blockID, solnData)
  if (sweepDir == SWEEP_X) then
     
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
        if (NDIM > 2) twodzInv = 1.0/(thirdCoord(k+1)-thirdCoord(k-1))
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           if (NDIM > 1) twodyInv = 1.0/(secondCoord(j+1)-secondCoord(j-1))
           
           u(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)) = &
                solnData(VELX_VAR, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), j, k)
           ut(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)) = &
                solnData(VELY_VAR, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), j, k)
           utt(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)) = &
                solnData(VELZ_VAR, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), j, k)
           rho(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)) = &
                solnData(DENS_VAR, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), j, k)


           do i = blkLimits(LOW,IAXIS)-1, blkLimits(HIGH,IAXIS)+1

              xtemp = solnData(TEMP_VAR, i, j, k)
              xdens = solnData(DENS_VAR, i, j, k)
              
              ! load the mass fractions
              do n = 1, NSPECIES
                 massfrac(n) = solnData(SPECIES_BEGIN-1+n, i, j, k)
              enddo
              
              call Viscosity(xtemp,xdens,massfrac,viscDynamic, viscUnusedHere)
              
              
              visc(i) = viscDynamic
           enddo
           
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1
              if (diff_geometricMeanDiff) then
                 visci(i) = 2.*(visc(i)*visc(i-1))/(visc(i)+visc(i-1))
              else
                 visci(i) = 0.5*(visc(i)+visc(i-1))
              endif
           enddo
           
           ! compute the fluxes
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1

              viscflux(i) =  dxinv(i)*( u(i)  - u(i-1)  ) * 4
              if (NDIM > 1) then
                 viscflux(i) = viscflux(i) - twodyInv * &
                      ( solnData(VELY_VAR, i-1, j+1, k) + solnData(VELY_VAR, i, j+1, k) &
                      - solnData(VELY_VAR, i-1, j-1, k) - solnData(VELY_VAR, i, j-1, k) )
              end if
              if (NDIM > 2) then
                 viscflux(i) = viscflux(i) - twodzInv * &
                      ( solnData(VELZ_VAR, i-1, j, k+1) + solnData(VELZ_VAR, i, j, k+1) &
                      - solnData(VELZ_VAR, i-1, j, k-1) - solnData(VELZ_VAR, i, j, k-1) )
              end if
              viscflux(i) =  -visci(i) * viscflux(i) / 3.0
              enerFlux = viscflux(i) * 0.5 * (u(i-1)+u(i))
              temp_flx(U_FLUX,i,j,k)   = temp_flx(U_FLUX,i,j,k) + areaLeft(i,j,k)*viscflux(i)


              if (NDIM > 1) then
                 viscflux(i) =  dxinv(i)*( ut(i)  - ut(i-1)  )
                 viscflux(i) = viscflux(i) + twodyInv * &
                      ( solnData(VELX_VAR, i-1, j+1, k) + solnData(VELX_VAR, i, j+1, k) &
                      - solnData(VELX_VAR, i-1, j-1, k) - solnData(VELX_VAR, i, j-1, k) ) * 0.5
                 viscflux(i) =  -visci(i) * viscflux(i)
                 enerFlux = enerFlux + viscflux(i) * 0.5 * (ut(i-1)+ut(i))
                 temp_flx(UT_FLUX,i,j,k)   = temp_flx(UT_FLUX,i,j,k) + areaLeft(i,j,k)*viscflux(i)
              end if
              
              
              viscflux(i) =  -visci(i) * dxinv(i)*( utt(i)  - utt(i-1)*rho(i-1)/rho(i))
              
              if (NDIM > 2) then
                 viscflux(i) =  dxinv(i)*( utt(i)  - utt(i-1)  )
                 viscflux(i) = viscflux(i) + twodzInv * &
                      ( solnData(VELX_VAR, i-1, j, k+1) + solnData(VELX_VAR, i, j, k+1) &
                      - solnData(VELX_VAR, i-1, j, k-1) - solnData(VELX_VAR, i, j, k-1) ) * 0.5
                 viscflux(i) =  -visci(i) * viscflux(i)
                 enerFlux = enerFlux + viscflux(i) * 0.5 * (utt(i-1)+utt(i))
                 temp_flx(UTT_FLUX,i,j,k)   = temp_flx(UTT_FLUX,i,j,k) + areaLeft(i,j,k)*viscflux(i)
              end if
              temp_flx(E_FLUX,i,j,k)   = temp_flx(E_FLUX,i,j,k) + areaLeft(i,j,k)*enerFlux
#ifdef EINT_FLUX
              temp_flx(EINT_FLUX,i,j,k)   = temp_flx(EINT_FLUX,i,j,k) + areaLeft(i,j,k)*enerFlux
!!$              if (enerFlux > 0) then
!!$                 print *,'enerFlux OK!!',enerFlux
!!$              elseif (enerFlux < 0) then
!!$                 print *,'enerFlux ????',enerFlux
!!$              else
!!$                 print *,'enerFlux NULL',enerFlux
!!$              end if
#endif
              
              
           enddo
        enddo
     enddo
     
     !-----------------------------------------------------------------------------
     ! y-sweep
     !-----------------------------------------------------------------------------
  elseif (sweepDir .EQ. SWEEP_Y) then
     
     ! compute dx for this block -- in the y direction

     dx = rightCoords - leftCoords
     dxinv(:blkLimits(HIGH,JAXIS)+1) = 1.e0/dx(:blkLimits(HIGH,JAXIS)+1)
     
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        if (NDIM > 2) twodzInv = 1.0/(thirdCoord(k+1)-thirdCoord(k-1))
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           twodxInv = 1.0/(secondCoord(i+1)-secondCoord(i-1))
           
           u(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)) = &
                solnData(VELY_VAR, i, blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), k)
           ut(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)) = &
                solnData(VELX_VAR, i, blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), k)
           utt(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)) = &
                solnData(VELZ_VAR, i, blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), k)
           rho(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)) = &
                solnData(DENS_VAR, i, blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), k)
           


           
           ! compute the viscosity -- we need the viscosity in the guard cells
           ! immediately adjacent to the block we are operating on
           do j = blkLimits(LOW,JAXIS)-1,blkLimits(HIGH,JAXIS)+1
 
              xtemp = solnData(TEMP_VAR, i, j, k)
              xdens = solnData(DENS_VAR, i, j, k)
              
              ! load the mass fractions
              do n = 1, NSPECIES
                 massfrac(n) = solnData(SPECIES_BEGIN-1+n, i, j, k)
              enddo
              
              call Viscosity(xtemp,xdens,massfrac,viscDynamic,viscUnusedHere)
              
              visc(j) = viscDynamic
           enddo
           
           
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
              if (diff_geometricMeanDiff) then
                 visci(j) = 2.*(visc(j)*visc(j-1))/(visc(j)+visc(j-1))
              else
                 visci(j) = 0.5*(visc(j)+visc(j-1))
              endif
           enddo
           
           ! compute the fluxes
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
              
              viscflux(j) =  dxinv(j)*( u(j)  - u(j-1)  ) * 4
              viscflux(j) = viscflux(j) - twodxInv * &
                   ( solnData(VELX_VAR, i+1, j-1, k) + solnData(VELX_VAR, i+1, j, k) &
                   - solnData(VELX_VAR, i-1, j-1, k) - solnData(VELX_VAR, i-1, j, k) )
              if (NDIM > 2) then
                 viscflux(j) = viscflux(j) - twodzInv * &
                      ( solnData(VELZ_VAR, i, j-1, k+1) + solnData(VELZ_VAR, i, j, k+1) &
                      - solnData(VELZ_VAR, i, j-1, k-1) - solnData(VELZ_VAR, i, j, k-1) )
              end if
              viscflux(j) =  -visci(j) * viscflux(j) / 3.0
              enerFlux = viscflux(j) * 0.5 * (u(j-1)+u(j))
              temp_flx(U_FLUX,i,j,k)   = temp_flx(U_FLUX,i,j,k) + areaLeft(i,j,k)*viscflux(j)

              viscflux(j) =  dxinv(j) * ( ut(j)  - ut(j-1)  )
              viscflux(j) = viscflux(j) + twodxInv * &
                   ( solnData(VELY_VAR, i+1, j-1, k) + solnData(VELY_VAR, i+1, j, k) &
                   - solnData(VELY_VAR, i-1, j-1, k) - solnData(VELY_VAR, i-1, j, k) ) * 0.5
              viscflux(j) =  -visci(j) * viscflux(j)
              enerFlux = enerFlux + viscflux(j) * 0.5 * (ut(j-1)+ut(j))
              temp_flx(UT_FLUX,i,j,k)   = temp_flx(UT_FLUX,i,j,k) + areaLeft(i,j,k)*viscflux(j)

              if (NDIM > 2) then
                 viscflux(j) = dxinv(j) * ( utt(j) - utt(j-1) )
                 viscflux(j) = viscflux(j) + twodzInv * &
                      ( solnData(VELY_VAR, i, j-1, k+1) + solnData(VELY_VAR, i, j, k+1) &
                      - solnData(VELY_VAR, i, j-1, k-1) - solnData(VELY_VAR, i, j, k-1) ) * 0.5
                 viscflux(j) =  -visci(j) * viscflux(j)
                 enerFlux = enerFlux + viscflux(j) * 0.5 * (utt(j-1)+utt(j))
                 temp_flx(UTT_FLUX,i,j,k)   = temp_flx(UTT_FLUX,i,j,k) + areaLeft(i,j,k)*viscflux(j)
              end if
              temp_flx(E_FLUX,i,j,k)   = temp_flx(E_FLUX,i,j,k) + areaLeft(i,j,k)*enerFlux
#ifdef EINT_FLUX
              temp_flx(EINT_FLUX,i,j,k)   = temp_flx(EINT_FLUX,i,j,k) + areaLeft(i,j,k)*enerFlux
#endif
           enddo
        enddo
     enddo
     
     !-----------------------------------------------------------------------------
     ! z-sweep
     !-----------------------------------------------------------------------------
  else
     
     ! compute dx for this block -- in the z direction
     
     dx = rightCoords - leftCoords
     dxinv(:blkLimits(HIGH,KAXIS)+1) = 1.e0/dx(:blkLimits(HIGH,KAXIS)+1)
     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        twodyInv = 1.0/(thirdCoord(j+1)-thirdCoord(j-1))
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           twodxInv = 1.0/(secondCoord(i+1)-secondCoord(i-1))
           u(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) = &
                solnData(VELZ_VAR, i, j, blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
           ut(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) = &
                solnData(VELX_VAR, i, j,  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
           utt(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) = &
                solnData(VELY_VAR, i, j, blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
           rho(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) = &
                solnData(DENS_VAR, i, j, blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
           

           
           ! compute the viscosity -- we need the viscosity in the guard cells
           ! immediately adjacent to the block we are operating on
           do k = blkLimits(LOW,KAXIS)-1,blkLimits(HIGH,KAXIS)+1
              
              xtemp = solnData(TEMP_VAR, i, j, k)
              xdens = solnData(DENS_VAR, i, j, k)
              
              ! load the mass fractions
              do n = 1, NSPECIES
                 massfrac(n) = solnData(SPECIES_BEGIN-1+n, i, j, k)
              enddo
              
              call Viscosity(xtemp,xdens,massfrac,viscDynamic, viscUnusedHere)
              
              visc(k) = viscDynamic
              
           enddo
           
           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
              if (diff_geometricMeanDiff) then
                 visci(k) = 2.*(visc(k)*visc(k-1))/(visc(k)+visc(k-1))
              else
                 visci(k) = 0.5*(visc(k)+visc(k-1))
              endif
           enddo
           
           ! compute the fluxes
           do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
              
              viscflux(k) =  dxinv(k)*( u(k)  - u(k-1)  ) * 4
              viscflux(k) = viscflux(k) - twodxInv * &
                   ( solnData(VELX_VAR, i+1, j, k-1) + solnData(VELX_VAR, i+1, j, k) &
                   - solnData(VELX_VAR, i-1, j, k-1) - solnData(VELX_VAR, i-1, j, k) )
              viscflux(k) = viscflux(k) - twodyInv * &
                   ( solnData(VELY_VAR, i, j+1, k-1) + solnData(VELY_VAR, i, j+1, k) &
                   - solnData(VELY_VAR, i, j-1, k-1) - solnData(VELY_VAR, i, j-1, k) )
              viscflux(k) =  -visci(k) * viscflux(k) / 3.0
              enerFlux = viscflux(k) * 0.5 * (u(k-1)+u(k))
              temp_flx(U_FLUX,i,j,k)   = temp_flx(U_FLUX,i,j,k) + areaLeft(i,j,k)*viscflux(k)

              viscflux(k) =  dxinv(k) * ( ut(k)  - ut(k-1)  )
              viscflux(k) = viscflux(k) + twodxInv * &
                   ( solnData(VELZ_VAR, i+1, j, k-1) + solnData(VELZ_VAR, i+1, j, k) &
                   - solnData(VELZ_VAR, i-1, j, k-1) - solnData(VELZ_VAR, i-1, j, k) ) * 0.5
              viscflux(k) =  -visci(k) * viscflux(k)
              enerFlux = enerFlux + viscflux(k) * 0.5 * (ut(k-1)+ut(k))
              temp_flx(UT_FLUX,i,j,k)   = temp_flx(UT_FLUX,i,j,k) + areaLeft(i,j,k)*viscflux(k)

              viscflux(k) =  dxinv(k) * ( utt(k) - utt(k-1) )
              viscflux(k) = viscflux(k) + twodyInv * &
                   ( solnData(VELZ_VAR, i, j+1, k-1) + solnData(VELZ_VAR, i, j+1, k) &
                   - solnData(VELZ_VAR, i, j-1, k-1) - solnData(VELZ_VAR, i, j-1, k) ) * 0.5
              viscflux(k) =  -visci(k) * viscflux(k)
              enerFlux = enerFlux + viscflux(k) * 0.5 * (utt(k-1)+utt(k))
              temp_flx(UTT_FLUX,i,j,k)   = temp_flx(UTT_FLUX,i,j,k) + areaLeft(i,j,k)*viscflux(k)
              temp_flx(E_FLUX,i,j,k)   = temp_flx(E_FLUX,i,j,k) + areaLeft(i,j,k)*enerFlux
#ifdef EINT_FLUX
              temp_flx(EINT_FLUX,i,j,k)   = temp_flx(EINT_FLUX,i,j,k) + areaLeft(i,j,k)*enerFlux
#endif
           enddo
        enddo
     enddo
     
  endif
  
  call Grid_releaseBlkPtr(blockID,solnData)
#endif
!! End of #ifndef FLASH_HYDRO_UNSPLIT
  return
end subroutine Diffuse_visc
