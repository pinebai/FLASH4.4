!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_getElectricFields
!!
!! NAME
!!
!!  hy_uhd_getElecticFields
!!
!!
!! SYNOPSIS
!!
!!  hy_uhd_getElectricFields( integer(IN) :: blockID,
!!                            integer(IN) :: blkLimits(2,MDIM),
!!                            integer(IN) :: blkLimitsGC(2,MDIM),
!!                            integer(IN) :: del(MDIM),
!!                            real(IN)    :: flx(:,:,:,:),
!!                            real(IN)    :: fly(:,:,:,:),
!!                            real(IN)    :: flz(:,:,:,:))
!!
!! ARGUMENTS
!!
!!   blockID     - a local block ID
!!   blkLimits   - an array that holds the lower and upper indices of the section
!!                 of block without the guard cells
!!   blkLimitsGC - an array that holds the lower and upper indices of the section
!!                 of block with the guard cells
!!   del         - deltas in {x,y,z} directions
!!   flx         - face flux for x direction
!!   fly         - face flux for y direction
!!   flz         - face flux for z direction
!!
!! DESCRIPTION
!!
!!   This routine calculates electric fields based on the high-order Godunov fluxes.
!!   A choice of higher-order (formally third order) interpolation scheme to construct
!!   electric fields is also available by using Taylor expansion. 
!!   In flash.par file, the runtime parameter "E_modification" should be declared .TRUE.
!!   in order to use this modified electric field construction scheme.
!!
!!   Note that this routine is only required for multidimensional problems, and not
!!   needed for 1D simulations.
!!
!! REFERENCES
!!
!!  * Balsara & Spicer, JCP, 149:270-292, 1999
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!
!!***

!!REORDER(4):E
!!REORDER(4):U
!!REORDER(4): fl[xyz], B[xyz]

Subroutine hy_uhd_getElectricFields( blockID,blkLimits,blkLimitsGC,del,flx,fly,flz)

  use Hydro_data,       ONLY : hy_E_modification, hy_E_upwind
  use hy_uhd_interface, ONLY : hy_uhd_getFluxDeriv
  use Grid_interface,   ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr
  use hy_uhd_slopeLimiters, ONLY : signum

  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !!$ ----- Arguments type declaration ----------------------------
  integer, intent(IN)  :: blockID
  integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
  real,    intent(IN), dimension(MDIM)  :: del
  
#ifdef FIXEDBLOCKSIZE
       real, DIMENSION(NFLUXES,                  &
                       GRID_ILO_GC:GRID_IHI_GC,  &
                       GRID_JLO_GC:GRID_JHI_GC,  &
                       GRID_KLO_GC:GRID_KHI_GC), &
                       intent(IN) :: flx,fly,flz
#else
       real, DIMENSION(NFLUXES,                                         &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), &
                       intent(IN) :: flx,fly,flz
#endif
  !!$----------------------------------------------------------------


#if FLASH_NEDGE_VAR > 0
#if NDIM > 1

  integer :: i0,imax,j0,jmax,k0,kmax,kbeg,kend,i,j,k
  real    :: dx, dy, dz
  real    :: velx, vely, velz, magx, magy, magz
  real, pointer, dimension(:,:,:,:) :: Bx,By,Bz,E,U

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) &
       :: dxG,dyF,dx2G,dy2F
#if NDIM == 3
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) &
       :: dzF,dxH,dz2F,dx2H,dyH,dzG,dy2H,dz2G
#endif
#else
  real, dimension(blkLimitsGC(HIGH,IAXIS),   &
                  blkLimitsGC(HIGH,JAXIS),   &
                  blkLimitsGC(HIGH,KAXIS)) ::&
                  dxG,dyF,dx2G,dy2F
#if NDIM == 3
  real, dimension(blkLimitsGC(HIGH,IAXIS),   &
                  blkLimitsGC(HIGH,JAXIS),   &
                  blkLimitsGC(HIGH,KAXIS)) ::&
                  dzF,dxH,dz2F,dx2H,dyH,dzG,dy2H,dz2G
#endif
#endif

  real :: velxP,velxN,velyP,velyN,velzP,velzN
  real :: totVel,factor,eps=1.e-1,eps0=1.e-4,eps1=tiny(1.)
  real :: soundSpeed,pres,dens,gamc

  dx = del(DIR_X)
  dy = del(DIR_Y)

  i0   = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)
  j0   = blkLimits(LOW, JAXIS)
  jmax = blkLimits(HIGH,JAXIS)
  if (NDIM == 2) then
     k0   = 1
     kmax = 0
     kbeg = 1
     kend = 1
  elseif (NDIM == 3) then
     dz = del(DIR_Z)
     k0   = blkLimits(LOW, KAXIS)
     kmax = blkLimits(HIGH,KAXIS)
     kbeg = k0
     kend = kmax+1
  endif


  call Grid_getBlkPtr(blockID,E,SCRATCH)
  call Grid_getBlkPtr(blockID,U,CENTER)


  E(EZ_SCRATCH_GRID_VAR,:,:,:) = 0.
#if NDIM == 3
  E(EY_SCRATCH_GRID_VAR,:,:,:) = 0.
  E(EX_SCRATCH_GRID_VAR,:,:,:) = 0.
#endif


  if (hy_E_modification) then
     !!=====================================================================!!
     !! We compute cell-cornered electric fields using cell-centered fluxes !!
     !! by considering high order Taylor expansions to include derivatives  !!
     !!=====================================================================!!

#if NDIM == 3
     do k=k0-1,kmax+1
#elif NDIM == 2
     do k=1,1
#endif
        do j=j0-1,jmax+1
           do i=i0-1,imax+1
              call hy_uhd_getFluxDeriv(i,j,k,blkLimitsGC,DIR_X,DIR_Y,flx,dyF(i,j,k),dy2F(i,j,k))
              call hy_uhd_getFluxDeriv(i,j,k,blkLimitsGC,DIR_Y,DIR_X,fly,dxG(i,j,k),dx2G(i,j,k))
#if NDIM == 3
              call hy_uhd_getFluxDeriv(i,j,k,blkLimitsGC,DIR_X,DIR_Z,flx,dzF(i,j,k),dz2F(i,j,k))
              call hy_uhd_getFluxDeriv(i,j,k,blkLimitsGC,DIR_Z,DIR_X,flz,dxH(i,j,k),dx2H(i,j,k))

              call hy_uhd_getFluxDeriv(i,j,k,blkLimitsGC,DIR_Z,DIR_Y,flz,dyH(i,j,k),dy2H(i,j,k))
              call hy_uhd_getFluxDeriv(i,j,k,blkLimitsGC,DIR_Y,DIR_Z,fly,dzG(i,j,k),dz2G(i,j,k))
#endif
           enddo
        enddo
     enddo

     !!=========================================!!
     !! Now we are ready to get Electric fields !!
     !!=========================================!!

     do k=k0,kmax+1
        do j=j0,jmax+1
           do i=i0,imax+1

              if (.not. hy_E_upwind) then

                 E(EZ_SCRATCH_GRID_VAR, i,j,k) = &
                      .25*(- flx(F07MAGY_FLUX,i,  j-1,k  )+.5*dyF(i,  j-1,k)+.125*dy2F(i,  j-1,k) &
                           - flx(F07MAGY_FLUX,i,  j,  k  )-.5*dyF(i,  j,  k)+.125*dy2F(i,  j,  k) &
                           + fly(F06MAGX_FLUX,i-1,j,  k  )+.5*dxG(i-1,j,  k)+.125*dx2G(i-1,j,  k) &
                           + fly(F06MAGX_FLUX,i,  j,  k  )-.5*dxG(i,  j,  k)+.125*dx2G(i,  j,  k))
#if NDIM == 3
                 E(EY_SCRATCH_GRID_VAR,i,j,k) = &
                      .25*( flx(F08MAGZ_FLUX,i,  j,  k-1)+.5*dzF(i,  j,k-1)+.125*dz2F(i,  j,k-1) &
                           +flx(F08MAGZ_FLUX,i,  j,  k  )-.5*dzF(i,  j,k  )+.125*dz2F(i,  j,k  ) &
                           -flz(F06MAGX_FLUX,i-1,j,  k  )+.5*dxH(i-1,j,k  )+.125*dx2H(i-1,j,k  ) &
                           -flz(F06MAGX_FLUX,i,  j,  k  )-.5*dxH(i,  j,k  )+.125*dx2H(i,  j,k  ))

                 E(EX_SCRATCH_GRID_VAR,i,j,k) = &
                      .25*( flz(F07MAGY_FLUX,i,  j-1,k  )+.5*dyH(i,j-1,k  )+.125*dy2H(i,j-1,k  ) &
                           +flz(F07MAGY_FLUX,i,  j,  k  )-.5*dyH(i,j,  k  )+.125*dy2H(i,j,  k  ) &
                           -fly(F08MAGZ_FLUX,i,  j,  k-1)+.5*dzG(i,j,  k-1)+.125*dz2G(i,j,  k-1) &
                           -fly(F08MAGZ_FLUX,i,  j,  k  )-.5*dzG(i,j,  k  )+.125*dz2G(i,j,  k  ))
#endif
              else !if (hy_E_upwind)

                 ! check local velocities for proper upwindings: 
                 velx = (U(VELX_VAR,i,j,k)+U(VELX_VAR,i-1,j,k)+U(VELX_VAR,i,j-1,k)+U(VELX_VAR,i-1,j-1,k))*0.25
                 vely = (U(VELY_VAR,i,j,k)+U(VELY_VAR,i-1,j,k)+U(VELY_VAR,i,j-1,k)+U(VELY_VAR,i-1,j-1,k))*0.25

                 ! check the local sound speed
                 pres = (U(PRES_VAR,i,j,k)+U(PRES_VAR,i-1,j,k)+U(PRES_VAR,i,j-1,k)+U(PRES_VAR,i-1,j-1,k))*0.25
                 dens = (U(DENS_VAR,i,j,k)+U(DENS_VAR,i-1,j,k)+U(DENS_VAR,i,j-1,k)+U(DENS_VAR,i-1,j-1,k))*0.25
                 gamc = (U(GAMC_VAR,i,j,k)+U(GAMC_VAR,i-1,j,k)+U(GAMC_VAR,i,j-1,k)+U(GAMC_VAR,i-1,j-1,k))*0.25
                 soundspeed = sqrt(gamc*pres/dens)

                 totVel = sqrt(velx**2 + vely**2)

                 if (totVel/soundSpeed .le. eps0) then
                    factor = 0.25
                    velxP = 1.
                    velxN = 1.
                    velyP = 1.
                    velyN = 1.
                 else
                    if (abs(velx)/max(totVel,eps1) .le. eps) velx=0.
                    if (abs(vely)/max(totVel,eps1) .le. eps) vely=0.

                    velxP=0.5*(1.+signum(velx))*abs(signum(velx))
                    velxN=0.5*(1.-signum(velx))*abs(signum(velx))
                    velyP=0.5*(1.+signum(vely))*abs(signum(vely))
                    velyN=0.5*(1.-signum(vely))*abs(signum(vely))

                    if (velx*vely == 0.) then
                       factor = 1.
                    else
                       factor = 0.5
                    endif
                 endif


                 E(EZ_SCRATCH_GRID_VAR, i,j,k) = &
                      factor* (  velyP*(-flx(F07MAGY_FLUX,i,  j-1,k  )+.5*dyF(i,  j-1,k)+.125*dy2F(i,  j-1,k)) &
                                +velyN*(-flx(F07MAGY_FLUX,i,  j,  k  )-.5*dyF(i,  j,  k)+.125*dy2F(i,  j,  k)) &
                                +velxP*( fly(F06MAGX_FLUX,i-1,j,  k  )+.5*dxG(i-1,j,  k)+.125*dx2G(i-1,j,  k)) &
                                +velxN*( fly(F06MAGX_FLUX,i,  j,  k  )-.5*dxG(i,  j,  k)+.125*dx2G(i,  j,  k)))
#if NDIM == 3
                 ! check local velocities for proper upwindings
                 velx = (U(VELX_VAR,i,j,k)+U(VELX_VAR,i-1,j,k)+U(VELX_VAR,i,j,k-1)+U(VELX_VAR,i-1,j,k-1))*0.25
                 velz = (U(VELZ_VAR,i,j,k)+U(VELZ_VAR,i-1,j,k)+U(VELZ_VAR,i,j,k-1)+U(VELZ_VAR,i-1,j,k-1))*0.25

                 ! check the local sound speed
                 pres = (U(PRES_VAR,i,j,k)+U(PRES_VAR,i-1,j,k)+U(PRES_VAR,i,j,k-1)+U(PRES_VAR,i-1,j,k-1))*0.25
                 dens = (U(DENS_VAR,i,j,k)+U(DENS_VAR,i-1,j,k)+U(DENS_VAR,i,j,k-1)+U(DENS_VAR,i-1,j,k-1))*0.25
                 gamc = (U(GAMC_VAR,i,j,k)+U(GAMC_VAR,i-1,j,k)+U(GAMC_VAR,i,j,k-1)+U(GAMC_VAR,i-1,j,k-1))*0.25
                 soundspeed = sqrt(gamc*pres/dens)

                 totVel = sqrt(velx**2 + velz**2)

                 if (totVel/soundSpeed .le. eps0) then
                    factor = 0.25
                    velxP = 1.
                    velxN = 1.
                    velzP = 1.
                    velzN = 1.
                 else
                    if (abs(velx)/max(totVel,eps1) .le. eps) velx=0.
                    if (abs(velz)/max(totVel,eps1) .le. eps) velz=0.

                    velxP=0.5*(1.+signum(velx))*abs(signum(velx))
                    velxN=0.5*(1.-signum(velx))*abs(signum(velx))
                    velzP=0.5*(1.+signum(velz))*abs(signum(velz))
                    velzN=0.5*(1.-signum(velz))*abs(signum(velz))

                    if (velx*velz == 0.) then
                       factor = 1.
                    else
                       factor = 0.5
                    endif
                 endif

                 E(EY_SCRATCH_GRID_VAR,i,j,k) = &
                      factor* (  velzP*( flx(F08MAGZ_FLUX,i,  j,  k-1)+.5*dzF(i,  j,k-1)+.125*dz2F(i,  j,k-1)) &
                                +velzN*( flx(F08MAGZ_FLUX,i,  j,  k  )-.5*dzF(i,  j,k  )+.125*dz2F(i,  j,k  )) &
                                +velxP*(-flz(F06MAGX_FLUX,i-1,j,  k  )+.5*dxH(i-1,j,k  )+.125*dx2H(i-1,j,k  )) &
                                +velxN*(-flz(F06MAGX_FLUX,i,  j,  k  )-.5*dxH(i,  j,k  )+.125*dx2H(i,  j,k  )))

                 ! check local velocities for proper upwindings
                 vely = (U(VELY_VAR,i,j,k)+U(VELY_VAR,i,j-1,k)+U(VELY_VAR,i,j,k-1)+U(VELY_VAR,i,j-1,k-1))*0.25
                 velz = (U(VELZ_VAR,i,j,k)+U(VELZ_VAR,i,j-1,k)+U(VELZ_VAR,i,j,k-1)+U(VELZ_VAR,i,j-1,k-1))*0.25

                 ! check the local sound speed
                 pres = (U(PRES_VAR,i,j,k)+U(PRES_VAR,i,j-1,k)+U(PRES_VAR,i,j,k-1)+U(PRES_VAR,i,j-1,k-1))*0.25
                 dens = (U(DENS_VAR,i,j,k)+U(DENS_VAR,i,j-1,k)+U(DENS_VAR,i,j,k-1)+U(DENS_VAR,i,j-1,k-1))*0.25
                 gamc = (U(GAMC_VAR,i,j,k)+U(GAMC_VAR,i,j-1,k)+U(GAMC_VAR,i,j,k-1)+U(GAMC_VAR,i,j-1,k-1))*0.25
                 soundspeed = sqrt(gamc*pres/dens)

                 totVel = sqrt(vely**2 + velz**2)

                 if (totVel/soundSpeed .le. eps0) then
                    factor = 0.25
                    velyP = 1.
                    velyN = 1.
                    velzP = 1.
                    velzN = 1.
                 else
                    if (abs(vely)/max(totVel,eps1) .le. eps) vely=0.
                    if (abs(velz)/max(totVel,eps1) .le. eps) velz=0.

                    velyP=0.5*(1.+signum(vely))*abs(signum(vely))
                    velyN=0.5*(1.-signum(vely))*abs(signum(vely))
                    velzP=0.5*(1.+signum(velz))*abs(signum(velz))
                    velzN=0.5*(1.-signum(velz))*abs(signum(velz))

                    if (vely*velz == 0.) then
                       factor = 1.
                    else
                       factor = 0.5
                    endif
                 endif

                 E(EX_SCRATCH_GRID_VAR,i,j,k) = &
                      factor* (  velyP*( flz(F07MAGY_FLUX,i,  j-1,k  )+.5*dyH(i,j-1,k  )+.125*dy2H(i,j-1,k  )) &
                                +velyN*( flz(F07MAGY_FLUX,i,  j,  k  )-.5*dyH(i,j,  k  )+.125*dy2H(i,j,  k  )) &
                                +velzP*(-fly(F08MAGZ_FLUX,i,  j,  k-1)+.5*dzG(i,j,  k-1)+.125*dz2G(i,j,  k-1)) &
                                +velzN*(-fly(F08MAGZ_FLUX,i,  j,  k  )-.5*dzG(i,j,  k  )+.125*dz2G(i,j,  k  )))
#endif
              endif !end of if (.not. hy_E_upwind) then

           enddo
        enddo
     enddo

  else !if hy_E_modification is turned off (Balsara and Spicer's simple arithmetic averaging)
     
     do k=k0,kmax+1
        do j=j0,jmax+1
           do i=i0,imax+1
              if (.not. hy_E_upwind) then

                 E(EZ_SCRATCH_GRID_VAR,i,j,k) = &
                      .25*(-flx(F07MAGY_FLUX,i,  j-1,k  )  &
                           -flx(F07MAGY_FLUX,i,  j,  k  )  &
                           +fly(F06MAGX_FLUX,i-1,j,  k  )  &
                           +fly(F06MAGX_FLUX,i,  j,  k  ) )
#if NDIM == 3
                 E(EY_SCRATCH_GRID_VAR,i,j,k) = &
                      .25*( flx(F08MAGZ_FLUX,i,  j,  k-1)  &
                           +flx(F08MAGZ_FLUX,i,  j,  k  )  &
                           -flz(F06MAGX_FLUX,i-1,j,  k  )  &
                           -flz(F06MAGX_FLUX,i,  j,  k  ) )

                 E(EX_SCRATCH_GRID_VAR, i,j,k) = &
                      .25*( flz(F07MAGY_FLUX,i,  j-1,k  )  &
                           +flz(F07MAGY_FLUX,i,  j,  k  )  &
                           -fly(F08MAGZ_FLUX,i,  j,  k-1)  &
                           -fly(F08MAGZ_FLUX,i,  j,  k  ) )
#endif
              else !! hy_E_upwind true

                 ! check local velocities for proper upwindings: 
                 velx = (U(VELX_VAR,i,j,k)+U(VELX_VAR,i-1,j,k)+U(VELX_VAR,i,j-1,k)+U(VELX_VAR,i-1,j-1,k))*0.25
                 vely = (U(VELY_VAR,i,j,k)+U(VELY_VAR,i-1,j,k)+U(VELY_VAR,i,j-1,k)+U(VELY_VAR,i-1,j-1,k))*0.25

                 ! check the local sound speed
                 pres = (U(PRES_VAR,i,j,k)+U(PRES_VAR,i-1,j,k)+U(PRES_VAR,i,j-1,k)+U(PRES_VAR,i-1,j-1,k))*0.25
                 dens = (U(DENS_VAR,i,j,k)+U(DENS_VAR,i-1,j,k)+U(DENS_VAR,i,j-1,k)+U(DENS_VAR,i-1,j-1,k))*0.25
                 gamc = (U(GAMC_VAR,i,j,k)+U(GAMC_VAR,i-1,j,k)+U(GAMC_VAR,i,j-1,k)+U(GAMC_VAR,i-1,j-1,k))*0.25
                 soundspeed = sqrt(gamc*pres/dens)

                 totVel = sqrt(velx**2 + vely**2)

                 if (totVel/soundSpeed .le. eps0) then
                    factor = 0.25
                    velxP = 1.
                    velxN = 1.
                    velyP = 1.
                    velyN = 1.
                 else
                    if (abs(velx)/max(totVel,eps1) .le. eps) velx=0.
                    if (abs(vely)/max(totVel,eps1) .le. eps) vely=0.

                    velxP=0.5*(1.+signum(velx))*abs(signum(velx))
                    velxN=0.5*(1.-signum(velx))*abs(signum(velx))
                    velyP=0.5*(1.+signum(vely))*abs(signum(vely))
                    velyN=0.5*(1.-signum(vely))*abs(signum(vely))

                    if (velx*vely == 0.) then
                       factor = 1.
                    else
                       factor = 0.5
                    endif
                 endif


                 E(EZ_SCRATCH_GRID_VAR, i,j,k) = &
                      factor* (  velyP*(-flx(F07MAGY_FLUX,i,  j-1,k  )) &
                                +velyN*(-flx(F07MAGY_FLUX,i,  j,  k  )) &
                                +velxP*( fly(F06MAGX_FLUX,i-1,j,  k  )) &
                                +velxN*( fly(F06MAGX_FLUX,i,  j,  k  )))
#if NDIM == 3
                 ! check local velocities for proper upwindings
                 velx = (U(VELX_VAR,i,j,k)+U(VELX_VAR,i-1,j,k)+U(VELX_VAR,i,j,k-1)+U(VELX_VAR,i-1,j,k-1))*0.25
                 velz = (U(VELZ_VAR,i,j,k)+U(VELZ_VAR,i-1,j,k)+U(VELZ_VAR,i,j,k-1)+U(VELZ_VAR,i-1,j,k-1))*0.25

                 ! check the local sound speed
                 pres = (U(PRES_VAR,i,j,k)+U(PRES_VAR,i-1,j,k)+U(PRES_VAR,i,j,k-1)+U(PRES_VAR,i-1,j,k-1))*0.25
                 dens = (U(DENS_VAR,i,j,k)+U(DENS_VAR,i-1,j,k)+U(DENS_VAR,i,j,k-1)+U(DENS_VAR,i-1,j,k-1))*0.25
                 gamc = (U(GAMC_VAR,i,j,k)+U(GAMC_VAR,i-1,j,k)+U(GAMC_VAR,i,j,k-1)+U(GAMC_VAR,i-1,j,k-1))*0.25
                 soundspeed = sqrt(gamc*pres/dens)

                 totVel = sqrt(velx**2 + velz**2)

                 if (totVel/soundSpeed .le. eps0) then
                    factor = 0.25
                    velxP = 1.
                    velxN = 1.
                    velzP = 1.
                    velzN = 1.
                 else
                    if (abs(velx)/max(totVel,eps1) .le. eps) velx=0.
                    if (abs(velz)/max(totVel,eps1) .le. eps) velz=0.

                    velxP=0.5*(1.+signum(velx))*abs(signum(velx))
                    velxN=0.5*(1.-signum(velx))*abs(signum(velx))
                    velzP=0.5*(1.+signum(velz))*abs(signum(velz))
                    velzN=0.5*(1.-signum(velz))*abs(signum(velz))

                    if (velx*velz == 0.) then
                       factor = 1.
                    else
                       factor = 0.5
                    endif
                 endif

                 E(EY_SCRATCH_GRID_VAR,i,j,k) = &
                      factor* (  velzP*( flx(F08MAGZ_FLUX,i,  j,  k-1)) &
                                +velzN*( flx(F08MAGZ_FLUX,i,  j,  k  )) &
                                +velxP*(-flz(F06MAGX_FLUX,i-1,j,  k  )) &
                                +velxN*(-flz(F06MAGX_FLUX,i,  j,  k  )))

                 ! check local velocities for proper upwindings
                 vely = (U(VELY_VAR,i,j,k)+U(VELY_VAR,i,j-1,k)+U(VELY_VAR,i,j,k-1)+U(VELY_VAR,i,j-1,k-1))*0.25
                 velz = (U(VELZ_VAR,i,j,k)+U(VELZ_VAR,i,j-1,k)+U(VELZ_VAR,i,j,k-1)+U(VELZ_VAR,i,j-1,k-1))*0.25

                 ! check the local sound speed
                 pres = (U(PRES_VAR,i,j,k)+U(PRES_VAR,i,j-1,k)+U(PRES_VAR,i,j,k-1)+U(PRES_VAR,i,j-1,k-1))*0.25
                 dens = (U(DENS_VAR,i,j,k)+U(DENS_VAR,i,j-1,k)+U(DENS_VAR,i,j,k-1)+U(DENS_VAR,i,j-1,k-1))*0.25
                 gamc = (U(GAMC_VAR,i,j,k)+U(GAMC_VAR,i,j-1,k)+U(GAMC_VAR,i,j,k-1)+U(GAMC_VAR,i,j-1,k-1))*0.25
                 soundspeed = sqrt(gamc*pres/dens)

                 totVel = sqrt(vely**2 + velz**2)

                 if (totVel/soundSpeed .le. eps0) then
                    factor = 0.25
                    velyP = 1.
                    velyN = 1.
                    velzP = 1.
                    velzN = 1.
                 else
                    if (abs(vely)/max(totVel,eps1) .le. eps) vely=0.
                    if (abs(velz)/max(totVel,eps1) .le. eps) velz=0.

                    velyP=0.5*(1.+signum(vely))*abs(signum(vely))
                    velyN=0.5*(1.-signum(vely))*abs(signum(vely))
                    velzP=0.5*(1.+signum(velz))*abs(signum(velz))
                    velzN=0.5*(1.-signum(velz))*abs(signum(velz))

                    if (vely*velz == 0.) then
                       factor = 1.
                    else
                       factor = 0.5
                    endif
                 endif

                 E(EX_SCRATCH_GRID_VAR,i,j,k) = &
                      factor* (  velyP*( flz(F07MAGY_FLUX,i,  j-1,k  )) &
                                +velyN*( flz(F07MAGY_FLUX,i,  j,  k  )) &
                                +velzP*(-fly(F08MAGZ_FLUX,i,  j,  k-1)) &
                                +velzN*(-fly(F08MAGZ_FLUX,i,  j,  k  )))
#endif
              
              endif
           enddo
        enddo
     enddo
  endif


  !! Release Block pointer
  call Grid_releaseBlkPtr(blockID,E,SCRATCH)
  call Grid_releaseBlkPtr(blockID,U,CENTER)

#endif
! end of #if NDIM > 0
#endif
! end of #if FLASH_NEDGE_VAR > 0

End Subroutine hy_uhd_getElectricFields
