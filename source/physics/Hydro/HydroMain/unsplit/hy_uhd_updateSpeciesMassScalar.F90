!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_updateSpeciesMassScalar
!!
!! NAME
!!
!!  hy_uhd_updateSpeciesMassScalar
!!
!! SYNOPSIS
!!
!!  hy_uhd_updateSpeciesMassScalar( integer(IN) :: order,
!!                                  real(IN)    :: densNew,
!!                                  real(IN)    :: Sp(NSPECIES+NMASS_SCALARS,-3:3,*),
!!                                  real(IN)    :: U(6,-3:3,...*),
!!                                  real(IN)    :: FL,
!!                                  real(IN)    :: FR,
!!                                  real(IN)    :: GL,
!!                                  real(IN)    :: GR,
!!                                  real(IN)    :: HL,
!!                                  real(IN)    :: HR,
!!                                  real(IN)    :: dx,
!!                                  real(IN)    :: dy,
!!                                  real(IN)    :: dz,
!!                                  real(IN)    :: dt,
!!                                  real(OUT)   :: SpNew(NSPECIES+NMASS_SCALARS))
!!
!! ARGUMENTS
!!
!!  order    - order for spatial discretization
!!  densNew  - density at n+1 step. If densNew is passed in as 0.0, then
!!             then the updated species and mass scalars are left in multiplied-by-density
!!             form (thus the species variables will contain partial densities), and
!!             species mass fractions will not have been renormalized to sum to 1.
!!  Sp       - Species and mass scalars to be advected at n step
!!  U        - solution vectors at n step containing density, pressure, velx vely, velz, and game in order
!!  FL,FR,GL,GR,HL,HR - Godunov fluxes of density
!!  dx,dy,dz - grid deltas
!!  dt       - timestep
!!  SpNew    - updated solution vectors of species and mass scalars at n+1 step
!!
!! DESCRIPTION
!!
!!  This routine advances species by locally an Eulerian algorithm in an unspit way.
!!
!!*** 
#define DENS 1
#define PRES 2
#define VELX 3
#define VELY 4
#define VELZ 5
#define GAME 6

!!REORDER(4): Sp, U

Subroutine hy_uhd_updateSpeciesMassScalar&
     (order,densNew,Sp,U,FL,FR,GL,GR,HL,HR,dx,dy,dz,dt,SpNew)

  use Hydro_data, ONLY : hy_ContactSteepening, hy_tiny, hy_transOrder
  use hy_uhd_slopeLimiters, ONLY : mc,signum

  implicit none

!#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  integer, intent(IN) :: order
  real, intent(IN) :: densNew
#if NDIM == 1
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,1,1)   :: Sp
  real, intent(IN), dimension(6,-3:3,1,1) :: U
#elif NDIM == 2
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,-3:3,1)   :: Sp
  real, intent(IN), dimension(6,-3:3,-3:3,1) :: U
#elif NDIM == 3
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,-3:3,-3:3)   :: Sp
  real, intent(IN), dimension(6,-3:3,-3:3,-3:3)   :: U
#endif

  real, intent(IN)  :: FL,FR,GL,GR,HL,HR,dx,dy,dz,dt
  real, intent(OUT), dimension(NSPECIES+NMASS_SCALARS) :: SpNew


  real :: delbarX0,delbarXL,delbarXR
  real :: delbarY0,delbarYL,delbarYR
  real :: delbarZ0,delbarZL,delbarZR

  real :: delbarX, delbarY,delbarZ

  real :: XL,XR
  real :: YL,YR
  real :: ZL,ZR

  real,dimension(NSPECIES+NMASS_SCALARS) :: XLh,XRh,YLh,YRh,ZLh,ZRh
  real,dimension(-1:1) :: X0p,X0m
  real :: XLp,XRm
  real,dimension(-1:1) :: Y0p,Y0m
  real :: YLp,YRm
  real,dimension(-1:1) :: Z0p,Z0m
  real :: ZLp,ZRm

  real :: delX,delY,delZ,X6,Y6,Z6
  real :: sig,absig,transFlux
  integer :: spi, spn,sps
  integer :: ix,iy,iz
  integer :: iShift
  real :: sumSpeciesXL, sumSpeciesXR
  real :: sumSpeciesYL, sumSpeciesYR
  real :: sumSpeciesZL, sumSpeciesZR
  real :: sumSpecies
  real :: hdtx,hdty,hdtz

  real :: corr3D

  !! For contact steepening
  real    :: temp1, temp2, temp3, eta_steep, del2rhoR, del2rhoL
  real, PARAMETER :: eta1=20.E0, eta2=0.05E0,epsln=0.01E0,K0=0.1E0



  !!           |           |
  !!        XLp|X0m  X0 X0p|XRm 
  !!           |           |
  !!-----*-----|-----*-----|-----*-----|
  !!         i-1/2   i   i+1/2  i+1  
  !!
  !!
  

  ix = 0
#if NDIM == 1
  iy = 1
  iz = 1
#elif NDIM ==2
  iy = 0
  iz = 1
#elif NDIM == 3
  iy = 0
  iz = 0
#endif


  ! Initialize arrays
  X0m = 0.
  X0p = 0.
  Y0m = 0.
  Y0p = 0.
  Z0m = 0.
  Z0p = 0.

  ! Note:
  ! spn = NSPECIES+NMASS_SCALARS
  ! spi = spn-NPROP_VARS
  ! sps = NSPECIES

  sumSpeciesXL = 0.
  sumSpeciesXR = 0.
  sumSpeciesYL = 0.
  sumSpeciesYR = 0.
  sumSpeciesZL = 0.
  sumSpeciesZR = 0.

  hdtx = 0.5*dt/dx
  if (NDIM > 1) then
     hdty = 0.5*dt/dy
     if (NDIM > 2) then
        hdtz = 0.5*dt/dz
     endif
  endif

  do spn = SPECIES_BEGIN, MASS_SCALARS_END
     spi = spn-NPROP_VARS

     if (order < 3) then !! MUSCL-HANCOCK TYPE 2ND ORDER METHOD================================

        !! Use MC slope limiter
        delbarXL=mc(Sp(spi,ix-1,iy,iz)-Sp(spi,ix-2,iy,iz),Sp(spi,ix,  iy,iz)-Sp(spi,ix-1,iy,iz))
        delbarX0=mc(Sp(spi,ix,  iy,iz)-Sp(spi,ix-1,iy,iz),Sp(spi,ix+1,iy,iz)-Sp(spi,ix,  iy,iz))
        delbarXR=mc(Sp(spi,ix+1,iy,iz)-Sp(spi,ix,  iy,iz),Sp(spi,ix+2,iy,iz)-Sp(spi,ix+1,iy,iz))

        ! (1) Left state of center node at i-1/2  ---------------------------
        X0m(1)=Sp(spi,ix,iy,iz) -(0.5+min(U(VELX,ix,iy,iz),0.)*hdtx)*delbarX0

        ! (2) Right state of center node at i+1/2  --------------------------
        X0p(1)=Sp(spi,ix,iy,iz) +(0.5-max(U(VELX,ix,iy,iz),0.)*hdtx)*delbarX0

        if (NDIM >= 2) then
           ! y-flux
           sig = signum(U(VELY,ix,iy,iz))
           absig = abs(sig)

           call computeTransFlux(Sp(spi,ix,iy-2:iy+2,iz),hy_transOrder,hdty,sig,transFlux)
           X0m(1) = X0m(1) + transFlux*absig*U(VELY,ix,iy,iz)
           X0p(1) = X0p(1) + transFlux*absig*U(VELY,ix,iy,iz)

           ! z-flux
           if (NDIM == 3) then
              sig = signum(U(VELZ,ix,iy,iz))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix,iy,iz-2:iz+2),hy_transOrder,hdtz,sig,transFlux)
              X0m(1) = X0m(1) + transFlux*absig*U(VELZ,ix,iy,iz)
              X0p(1) = X0p(1) + transFlux*absig*U(VELZ,ix,iy,iz)

              ! y-z cross deriv at ix
              call compute2dCrossDerivative(U(VELY,ix,iy,iz),U(VELZ,ix,iy,iz),Sp(spi,ix,iy-1:iy+1,iz-1:iz+1),corr3D)
              corr3D = corr3D*4./3.*hdty*hdtz

              X0m(1) = X0m(1) + corr3D
              X0p(1) = X0p(1) + corr3D

           endif
        endif


        ! (3) Right state of left node ----------------------------------
        XLp=Sp(spi,ix-1,iy,iz)+(0.5-max(U(VELX,ix-1,iy,iz),0.)*hdtx)*delbarXL
        if (NDIM >= 2) then
           ! y-flux
           sig = signum(U(VELY,ix-1,iy,iz))
           absig = abs(sig)

           call computeTransFlux(Sp(spi,ix-1,iy-2:iy+2,iz),hy_transOrder,hdty,sig,transFlux)
           XLp = XLp + transFlux*absig*U(VELY,ix-1,iy,iz)

           if (NDIM == 3) then
              ! z-flux
              sig = signum(U(VELZ,ix-1,iy,iz))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix-1,iy,iz-2:iz+2),hy_transOrder,hdtz,sig,transFlux)
              XLp = XLp + transFlux*absig*U(VELZ,ix-1,iy,iz)

              ! y-z cross deriv at ix-1
              call compute2dCrossDerivative(U(VELY,ix-1,iy,iz),U(VELZ,ix-1,iy,iz),Sp(spi,ix-1,iy-1:iy+1,iz-1:iz+1),corr3D)
              corr3D =  corr3D*4./3.*hdty*hdtz

              XLp = XLp + corr3D

           endif
        endif

        !! Upwinding
        sig = signum(FL)
        XLh(spi) = 0.5*(1.+sig)*XLp+0.5*(1.-sig)*X0m(1)


        ! (4) Left state of right node ---------------------------------
        XRm=Sp(spi,ix+1,iy,iz)-(0.5+min(U(VELX,ix+1,iy,iz),0.)*hdtx)*delbarXR

        if (NDIM >= 2) then
           ! y-flux
           sig = signum(U(VELY,ix+1,iy,iz))
           absig = abs(sig)

           call computeTransFlux(Sp(spi,ix+1,iy-2:iy+2,iz),hy_transOrder,hdty,sig,transFlux)
           XRm = XRm + transFlux*absig*U(VELY,ix+1,iy,iz)

           if (NDIM == 3) then
              ! z-flux
              sig = signum(U(VELZ,ix+1,iy,iz))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix+1,iy,iz-2:iz+2),hy_transOrder,hdtz,sig,transFlux)
              XRm = XRm + transFlux*absig*U(VELZ,ix+1,iy,iz)

              ! y-z cross deriv at ix+1
              call compute2dCrossDerivative(U(VELY,ix+1,iy,iz),U(VELZ,ix+1,iy,iz),Sp(spi,ix+1,iy-1:iy+1,iz-1:iz+1),corr3D)
              corr3D =  corr3D*4./3.*hdty*hdtz

              XRm = XRm + corr3D

           endif
        endif


        !! Upwinding
        sig = signum(FR)
        XRh(spi) = 0.5*(1.+sig)*X0p(1)+0.5*(1.-sig)*XRm

     else !! PPM TYPE 3RD ORDER METHOD==================================================

        do iShift = -1,1
           !! Use MC slope limiter
           delbarXL=mc(Sp(spi,ix+iShift-1,iy,iz)-Sp(spi,ix+iShift-2,iy,iz),&
                       Sp(spi,ix+iShift,  iy,iz)-Sp(spi,ix+iShift-1,iy,iz))

           delbarX0=mc(Sp(spi,ix+iShift,  iy,iz)-Sp(spi,ix+iShift-1,iy,iz),&
                       Sp(spi,ix+iShift+1,iy,iz)-Sp(spi,ix+iShift,  iy,iz))

           delbarXR=mc(Sp(spi,ix+iShift+1,iy,iz)-Sp(spi,ix+iShift,  iy,iz),&
                       Sp(spi,ix+iShift+2,iy,iz)-Sp(spi,ix+iShift+1,iy,iz))


           !! PPM polynomial interpolation 
           ! (1) at i-1/2: Left state of center node (i,j,k)
           XL=0.5*(Sp(spi,ix+iShift,iy,iz)+Sp(spi,ix+iShift-1,iy,iz))-(delbarX0-delbarXL)/6.

           ! (2) at i+1/2: Right state of center node (i,j,k)
           XR=0.5*(Sp(spi,ix+iShift,iy,iz)+Sp(spi,ix+iShift+1,iy,iz))-(delbarXR-delbarX0)/6.

           !! Apply the same logic for steepening mass scalar/species 
           !! as for the contact steeping
           !! NOTE: THE FOLLOWING LINES ARE COMMENTED OUT, BUT CAN BE USED WHEN APPLYING STEEPENING.
!!$           temp1 = U(DENS,ix+iShift+1,iy,iz) - U(DENS,ix+iShift-1,iy,iz)
!!$           if (hy_ContactSteepening .and. abs(temp1) > hy_tiny) then
!!$              ! Eqn 1.17 : Second derivatives
!!$              del2rhoR = (    U(DENS,ix+iShift+2,iy,iz)&
!!$                          -2.*U(DENS,ix+iShift+1,iy,iz)&
!!$                            + U(DENS,ix+iShift,  iy,iz))/(6.*dx**2)
!!$
!!$              del2rhoL = (    U(DENS,ix+iShift,  iy,iz)&
!!$                          -2.*U(DENS,ix+iShift-1,iy,iz)&
!!$                             +U(DENS,ix+iShift-2,iy,iz))/(6.*dx**2)
!!$
!!$              ! Third derivative
!!$              eta_steep = (del2rhoL-del2rhoR)*dx**2/temp1
!!$              if (del2rhoR*del2rhoL >= 0.) then
!!$                 eta_steep = 0.
!!$              endif
!!$              if (epsln*min(U(DENS,ix+iShift+1,iy,iz),U(DENS,ix+iShift-1,iy,iz))&
!!$                   -abs(U(DENS,ix+iShift+1,iy,iz) - U(DENS,ix+iShift-1,iy,iz)) >= 0.) then
!!$                 eta_steep = 0.
!!$              endif
!!$
!!$              ! Eqn 1.16
!!$              eta_steep = max(0., min(1., eta1*(eta_steep - eta2)))
!!$
!!$              ! Eqn 3.2
!!$              temp2 = abs(U(PRES,ix+iShift+1,iy,iz)-U(PRES,ix+iShift-1,iy,iz))&
!!$              /min(U(PRES,ix+iShift+1,iy,iz),U(PRES,ix+iShift-1,iy,iz))
!!$
!!$              temp3 = abs(U(DENS,ix+iShift+1,iy,iz)-U(DENS,ix+iShift-1,iy,iz))&
!!$                   /min(U(DENS,ix+iShift+1,iy,iz),U(DENS,ix+iShift-1,iy,iz))
!!$
!!$              if (U(GAME,ix+iShift,iy,iz)*K0*temp3-temp2 < 0.0) then
!!$                 eta_steep = 0.
!!$              endif
!!$
!!$              ! Eqn 1.15
!!$              XL = XL*(1.-eta_steep) + (Sp(spi,ix+iShift-1,iy,iz)+0.5*delbarXL)*eta_steep
!!$              XR = XR*(1.-eta_steep) + (Sp(spi,ix+iShift+1,iy,iz)-0.5*delbarXR)*eta_steep
!!$           endif

           !! Apply limiting
           XL=max(min(Sp(spi,ix+iShift,iy,iz),Sp(spi,ix+iShift-1,iy,iz)),&
                min(max(Sp(spi,ix+iShift,iy,iz),Sp(spi,ix+iShift-1,iy,iz)),XL))

           !! Apply limiting
           XR=max(min(Sp(spi,ix+iShift,iy,iz),Sp(spi,ix+iShift+1,iy,iz)),&
                min(max(Sp(spi,ix+iShift,iy,iz),Sp(spi,ix+iShift+1,iy,iz)),XR))

           !! Check monotonicity
           if ((XR-Sp(spi,ix+iShift,iy,iz))*(Sp(spi,ix+iShift,iy,iz)-XL)<=0.) then
              XR=Sp(spi,ix+iShift,iy,iz)
              XL=Sp(spi,ix+iShift,iy,iz)
           endif
           if (6.*(XR-XL)*(Sp(spi,ix+iShift,iy,iz)-0.5*(XL+XR)) > (XR-XL)**2) then
              XL = 3.*Sp(spi,ix+iShift,iy,iz)-2.*XR
           endif
           if (6.*(XR-XL)*(Sp(spi,ix+iShift,iy,iz)-0.5*(XL+XR)) <-(XR-XL)**2) then
              XR = 3.*Sp(spi,ix+iShift,iy,iz)-2.*XL
           endif

           !! PPM reconstruction to advance 1/2 time step (i.e., n+1/2) at cell interfaces
           delX=XR-XL
           X6=6.*(Sp(spi,ix+iShift,iy,iz)-0.5*(XR+XL))

           X0m(iShift)=XL-min(U(VELX,ix+iShift,iy,iz),0.)*hdtx*&
                (delX+(1.0+min(U(VELX,ix+iShift,iy,iz),0.)*4./3.*hdtx)*X6)

           X0p(iShift)=XR-max(U(VELX,ix+iShift,iy,iz),0.)*hdtx*&
                (delX-(1.0-max(U(VELX,ix+iShift,iy,iz),0.)*4./3.*hdtx)*X6)


           !! Transverse flux components correction-------------------------------------
           if (NDIM >= 2) then
              ! y-flux
              sig = signum(U(VELY,ix+iShift,iy,iz))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix+iShift,iy-2:iy+2,iz),hy_transOrder,hdty,sig,transFlux)
              X0m(iShift) = X0m(iShift) + transFlux*absig*U(VELY,ix+iShift,iy,iz)
              X0p(iShift) = X0p(iShift) + transFlux*absig*U(VELY,ix+iShift,iy,iz)

              ! z-flux
              if (NDIM == 3) then
                 sig = signum(U(VELZ,ix+iShift,iy,iz))
                 absig = abs(sig)

                 call computeTransFlux(Sp(spi,ix+iShift,iy,iz-2:iz+2),hy_transOrder,hdtz,sig,transFlux)
                 X0m(iShift) = X0m(iShift) + transFlux*absig*U(VELZ,ix+iShift,iy,iz)
                 X0p(iShift) = X0p(iShift) + transFlux*absig*U(VELZ,ix+iShift,iy,iz)

                 ! y-z cross deriv at ix+iShift
                 call compute2dCrossDerivative(U(VELY,ix+iShift,iy,iz),U(VELZ,ix+iShift,iy,iz),&
                                               Sp(spi,ix+iShift,iy-1:iy+1,iz-1:iz+1),corr3D)
                 corr3D = corr3D*4./3.*hdty*hdtz

                 X0m(iShift) = X0m(iShift) + corr3D
                 X0p(iShift) = X0p(iShift) + corr3D

              endif
           endif

        enddo

        !! Updates ----------------------------------------------------------------
        !! Upwinding at i-1/2
        sig = signum(FL)
        XLh(spi) = 0.5*(1.+sig)*X0p(-1)+0.5*(1.-sig)*X0m(0)

        !! Upwinding at i+1/2
        sig = signum(FR)
        XRh(spi) = 0.5*(1.+sig)*X0p(0)+0.5*(1.-sig)*X0m(1)

     endif

     !! Take a summation for total species
     if (spn <= SPECIES_END) then
        sumSpeciesXL = sumSpeciesXL + XLh(spi)
        sumSpeciesXR = sumSpeciesXR + XRh(spi)
     endif

     if (NDIM >= 2) then

        if (order < 3) then !! MUSCL-HANCOCK TYPE 2ND ORDER METHOD================================

           delbarYL=mc(Sp(spi,ix,iy-1,iz)-Sp(spi,ix,iy-2,iz),Sp(spi,ix,iy,  iz)-Sp(spi,ix,iy-1,iz))
           delbarY0=mc(Sp(spi,ix,iy,  iz)-Sp(spi,ix,iy-1,iz),Sp(spi,ix,iy+1,iz)-Sp(spi,ix,iy,  iz))
           delbarYR=mc(Sp(spi,ix,iy+1,iz)-Sp(spi,ix,iy,  iz),Sp(spi,ix,iy+2,iz)-Sp(spi,ix,iy+1,iz))

           ! (1) Left state of center node at j-1/2 ----------------------------
           Y0m(1)=Sp(spi,ix,iy,iz) -(0.5+min(U(VELY,ix,iy,iz),0.)*hdty)*delbarY0

           ! (3) Right state of center node at j+1/2 ---------------------------
           Y0p(1)=Sp(spi,ix,iy,iz) +(0.5-max(U(VELY,ix,iy,iz),0.)*hdty)*delbarY0

           ! x-flux
           sig = signum(U(VELX,ix,iy,iz))
           absig = abs(sig)

           call computeTransFlux(Sp(spi,ix-2:ix+2,iy,iz),hy_transOrder,hdtx,sig,transFlux)
           Y0m(1) = Y0m(1) + transFlux*absig*U(VELX,ix,iy,iz)
           Y0p(1) = Y0p(1) + transFlux*absig*U(VELX,ix,iy,iz)

           if (NDIM == 3) then
              ! z-flux
              sig = signum(U(VELZ,ix,iy,iz))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix,iy,iz-2:iz+2),hy_transOrder,hdtz,sig,transFlux)
              Y0m(1) = Y0m(1) + transFlux*absig*U(VELZ,ix,iy,iz)
              Y0p(1) = Y0p(1) + transFlux*absig*U(VELZ,ix,iy,iz)

              ! x-z cross deriv at iy
              call compute2dCrossDerivative(U(VELX,ix,iy,iz),U(VELZ,ix,iy,iz),Sp(spi,ix-1:ix+1,iy,iz-1:iz+1),corr3D)
              corr3D =  corr3D*4./3.*hdtx*hdtz

              Y0m(1) = Y0m(1) + corr3D
              Y0p(1) = Y0p(1) + corr3D

           endif


           ! (2) Right state of left node ---------------------------------
           YLp=Sp(spi,ix,iy-1,iz)+(0.5-max(U(VELY,ix,iy-1,iz),0.)*hdty)*delbarYL

           ! x-flux
           sig = signum(U(VELX,ix,iy-1,iz))
           absig = abs(sig)

           call computeTransFlux(Sp(spi,ix-2:ix+2,iy-1,iz),hy_transOrder,hdtx,sig,transFlux)
           YLp = YLp + transFlux*absig*U(VELX,ix,iy-1,iz)

           if (NDIM == 3) then
              ! z-flux
              sig = signum(U(VELZ,ix,iy-1,iz))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix,iy-1,iz-2:iz+2),hy_transOrder,hdtz,sig,transFlux)
              YLp = YLp + transFlux*absig*U(VELZ,ix,iy-1,iz)

              ! x-z cross deriv at iy-1
              call compute2dCrossDerivative(U(VELX,ix,iy-1,iz),U(VELZ,ix,iy-1,iz),Sp(spi,ix-1:ix+1,iy-1,iz-1:iz+1),corr3D)
              corr3D =  corr3D*4./3.*hdtx*hdtz

              YLp = YLp + corr3D

           endif

           !! Upwinding
           sig = signum(GL)
           YLh(spi) = 0.5*(1.+sig)*YLp+0.5*(1.-sig)*Y0m(1)


           ! (4) Left state of right node ---------------------------------
           YRm=Sp(spi,ix,iy+1,iz)-(0.5+min(U(VELY,ix,iy+1,iz),0.)*hdty)*delbarYR

           ! x-flux
           sig = signum(U(VELX,ix,iy+1,iz))
           absig = abs(sig)

           call computeTransFlux(Sp(spi,ix-2:ix+2,iy+1,iz),hy_transOrder,hdtx,sig,transFlux)
           YRm = YRm + transFlux*absig*U(VELX,ix,iy+1,iz)

           if (NDIM == 3) then
              ! z-flux
              sig = signum(U(VELZ,ix,iy+1,iz))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix,iy+1,iz-2:iz+2),hy_transOrder,hdtz,sig,transFlux)
              YRm = YRm + transFlux*absig*U(VELZ,ix,iy+1,iz)

              ! x-z cross deriv at iy+1
              call compute2dCrossDerivative(U(VELX,ix,iy+1,iz),U(VELZ,ix,iy+1,iz),Sp(spi,ix-1:ix+1,iy+1,iz-1:iz+1),corr3D)
              corr3D =  corr3D*4./3.*hdtx*hdtz

              YRm = YRm + corr3D



           endif

           !! Upwinding
           sig = signum(GR)
           YRh(spi) = 0.5*(1.+sig)*Y0p(1)+0.5*(1.-sig)*YRm

        else  !! PPM TYPE 3RD ORDER METHOD==================================================

           do iShift = -1,1
              !! Use MC slope limiter
              delbarYL=mc(Sp(spi,ix,iy+iShift-1,iz)-Sp(spi,ix,iy+iShift-2,iz),&
                          Sp(spi,ix,iy+iShift,  iz)-Sp(spi,ix,iy+iShift-1,iz))

              delbarY0=mc(Sp(spi,ix,iy+iShift,iz  )-Sp(spi,ix,iy+iShift-1,iz),&
                          Sp(spi,ix,iy+iShift+1,iz)-Sp(spi,ix,iy+iShift,  iz))

              delbarYR=mc(Sp(spi,ix,iy+iShift+1,iz)-Sp(spi,ix,iy+iShift,iz  ),&
                          Sp(spi,ix,iy+iShift+2,iz)-Sp(spi,ix,iy+iShift+1,iz))


              !! PPM polynomial interpolation 
              ! (1) at i-1/2: Left state of center node (i,j,k)
              YL=0.5*(Sp(spi,ix,iy+iShift,iz)+Sp(spi,ix,iy+iShift-1,iz))-(delbarY0-delbarYL)/6.

              ! (2) at i+1/2: Right state of center node (i,j,k)
              YR=0.5*(Sp(spi,ix,iy+iShift,iz)+Sp(spi,ix,iy+iShift+1,iz))-(delbarYR-delbarY0)/6.

              !! Apply the same logic for steepening mass scalar/species 
              !! as for the contact steeping
              !! NOTE: THE FOLLOWING LINES ARE COMMENTED OUT, BUT CAN BE USED WHEN APPLYING STEEPENING.
              temp1 = U(DENS,ix,iy+iShift+1,iz) - U(DENS,ix,iy+iShift-1,iz)
!!$              if (hy_ContactSteepening .and. abs(temp1) > hy_tiny) then
!!$                 ! Eqn 1.17 : Second derivatives
!!$                 del2rhoR = (    U(DENS,ix,iy+iShift+2,iz)&
!!$                             -2.*U(DENS,ix,iy+iShift+1,iz)&
!!$                               + U(DENS,ix,iy+iShift,  iz))/(6.*dy**2)
!!$
!!$                 del2rhoL = (    U(DENS,ix,iy+iShift,  iz)&
!!$                             -2.*U(DENS,ix,iy+iShift-1,iz)&
!!$                                +U(DENS,ix,iy+iShift-2,iz))/(6.*dy**2)
!!$
!!$                 ! Third derivative
!!$                 eta_steep = (del2rhoL-del2rhoR)*dy**2/temp1
!!$                 if (del2rhoR*del2rhoL >= 0.) then
!!$                    eta_steep = 0.
!!$                 endif
!!$                 if (epsln*min(U(DENS,ix,iy+iShift+1,iz),U(DENS,ix,iy+iShift-1,iz))&
!!$                      -abs(U(DENS,ix,iy+iShift+1,iz) - U(DENS,ix,iy+iShift-1,iz)) >= 0.) then
!!$                    eta_steep = 0.
!!$                 endif
!!$
!!$                 ! Eqn 1.16
!!$                 eta_steep = max(0., min(1., eta1*(eta_steep - eta2)))
!!$
!!$                 ! Eqn 3.2
!!$                 temp2 = abs(U(PRES,ix,iy+iShift+1,iz)-U(PRES,ix,iy+iShift-1,iz))&
!!$                 /min(U(PRES,ix,iy+iShift+1,iz),U(PRES,ix,iy+iShift-1,iz))
!!$
!!$                 temp3 = abs(U(DENS,ix,iy+iShift+1,iz)-U(DENS,ix,iy+iShift-1,iz))&
!!$                      /min(U(DENS,ix,iy+iShift+1,iz),U(DENS,ix,iy+iShift-1,iz))
!!$
!!$                 if (U(GAME,ix,iy+iShift,iz)*K0*temp3-temp2 < 0.0) then
!!$                    eta_steep = 0.
!!$                 endif
!!$
!!$                 ! Eqn 1.15
!!$                 YL = YL*(1.-eta_steep) + (Sp(spi,ix,iy+iShift-1,iz)+0.5*delbarYL)*eta_steep
!!$                 YR = YR*(1.-eta_steep) + (Sp(spi,ix,iy+iShift+1,iz)-0.5*delbarYR)*eta_steep
!!$              endif

              !! Apply limiting
              YL=max(min(Sp(spi,ix,iy+iShift,iz),Sp(spi,ix,iy+iShift-1,iz)),&
                   min(max(Sp(spi,ix,iy+iShift,iz),Sp(spi,ix,iy+iShift-1,iz)),YL))

              !! Apply limiting
              YR=max(min(Sp(spi,ix,iy+iShift,iz),Sp(spi,ix,iy+iShift+1,iz)),&
                   min(max(Sp(spi,ix,iy+iShift,iz),Sp(spi,ix,iy+iShift+1,iz)),YR))

              !! Check monotonicity
              if ((YR-Sp(spi,ix,iy+iShift,iz))*(Sp(spi,ix,iy+iShift,iz)-YL)<=0.) then
                 YR=Sp(spi,ix,iy+iShift,iz)
                 YL=Sp(spi,ix,iy+iShift,iz)
              endif
              if (6.*(YR-YL)*(Sp(spi,ix,iy+iShift,iz)-0.5*(YL+YR)) > (YR-YL)**2) then
                 YL = 3.*Sp(spi,ix,iy+iShift,iz)-2.*YR
              endif
              if (6.*(YR-YL)*(Sp(spi,ix,iy+iShift,iz)-0.5*(YL+YR)) <-(YR-YL)**2) then
                 YR = 3.*Sp(spi,ix,iy+iShift,iz)-2.*YL
              endif

              !! PPM reconstruction to advance 1/2 time step (i.e., n+1/2) at cell interfaces
              delY=YR-YL
              Y6=6.*(Sp(spi,ix,iy+iShift,iz)-0.5*(YR+YL))

              Y0m(iShift)=YL-min(U(VELY,ix,iy+iShift,iz),0.)*hdty*&
                   (delY+(1.0+min(U(VELY,ix,iy+iShift,iz),0.)*4./3.*hdty)*Y6)

              Y0p(iShift)=YR-max(U(VELY,ix,iy+iShift,iz),0.)*hdty*&
                   (delY-(1.0-max(U(VELY,ix,iy+iShift,iz),0.)*4./3.*hdty)*Y6)


              !! Transverse flux components correction-------------------------------------
              ! x-flux
              sig = signum(U(VELX,ix,iy+iShift,iz)) ! it was wrong in the original version, missing iShift
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix-2:ix+2,iy+iShift,iz),hy_transOrder,hdtx,sig,transFlux)
              Y0m(iShift) = Y0m(iShift) + transFlux*absig*U(VELX,ix,iy+iShift,iz)
              Y0p(iShift) = Y0p(iShift) + transFlux*absig*U(VELX,ix,iy+iShift,iz)

              if (NDIM == 3) then
                 ! z-flux
                 sig = signum(U(VELZ,ix,iy+iShift,iz))
                 absig = abs(sig)

                 call computeTransFlux(Sp(spi,ix,iy+iShift,iz-2:iz+2),hy_transOrder,hdtz,sig,transFlux)
                 Y0m(iShift) = Y0m(iShift) + transFlux*absig*U(VELZ,ix,iy+iShift,iz)
                 Y0p(iShift) = Y0p(iShift) + transFlux*absig*U(VELZ,ix,iy+iShift,iz)


                 ! x-z cross deriv at iy+iShift
                 call compute2dCrossDerivative(U(VELX,ix,iy+iShift,iz),U(VELZ,ix,iy+iShift,iz),&
                                               Sp(spi,ix-1:ix+1,iy+iShift,iz-1:iz+1),corr3D)
                 corr3D = corr3D*4./3.*hdtx*hdtz

                 Y0m(iShift) = Y0m(iShift) + corr3D
                 Y0p(iShift) = Y0p(iShift) + corr3D

              endif

           enddo


           !! Updates ---------------------------------------------------------------------
           !! Upwinding at i-1/2
           sig = signum(GL)
           YLh(spi) = 0.5*(1.+sig)*Y0p(-1)+0.5*(1.-sig)*Y0m(0)

           !! Upwinding at i+1/2
           sig = signum(GR)
           YRh(spi) = 0.5*(1.+sig)*Y0p(0)+0.5*(1.-sig)*Y0m(1)

        endif

        if (spn <= SPECIES_END) then
           sumSpeciesYL = sumSpeciesYL + YLh(spi)
           sumSpeciesYR = sumSpeciesYR + YRh(spi)
        endif



        if (NDIM == 3) then

           if (order < 3) then !! MUSCL-HANCOCK TYPE 2ND ORDR METHOD================================

              delbarZL=mc(Sp(spi,ix,iy,iz-1)-Sp(spi,ix,iy,iz-2),Sp(spi,ix,iy,iz  )-Sp(spi,ix,iy,iz-1))
              delbarZ0=mc(Sp(spi,ix,iy,iz  )-Sp(spi,ix,iy,iz-1),Sp(spi,ix,iy,iz+1)-Sp(spi,ix,iy,iz  ))
              delbarZR=mc(Sp(spi,ix,iy,iz+1)-Sp(spi,ix,iy,iz  ),Sp(spi,ix,iy,iz+2)-Sp(spi,ix,iy,iz+1))

              ! (1) Left state of center node at i-1/2 ---------------------------
              Z0m(1)=Sp(spi,ix,iy,iz) -(0.5+min(U(VELZ,ix,iy,iz),0.)*hdtz)*delbarZ0

              ! (2) Right state of center node at i+1/2 ---------------------------
              Z0p(1)=Sp(spi,ix,iy,iz) +(0.5-max(U(VELZ,ix,iy,iz),0.)*hdtz)*delbarZ0

              ! x-flux
              sig = signum(U(VELX,ix,iy,iz))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix-2:ix+2,iy,iz),hy_transOrder,hdtx,sig,transFlux)
              Z0m(1) = Z0m(1) + transFlux*absig*U(VELX,ix,iy,iz)
              Z0p(1) = Z0p(1) + transFlux*absig*U(VELX,ix,iy,iz)

              ! y-flux
              sig = signum(U(VELY,ix,iy,iz))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix,iy-2:iy+2,iz),hy_transOrder,hdty,sig,transFlux)
              Z0m(1) = Z0m(1) + transFlux*absig*U(VELY,ix,iy,iz)
              Z0p(1) = Z0p(1) + transFlux*absig*U(VELY,ix,iy,iz)

              ! x-y cross deriv at iz
              call compute2dCrossDerivative(U(VELX,ix,iy,iz),U(VELY,ix,iy,iz),Sp(spi,ix-1:ix+1,iy-1:iy+1,iz),corr3D)
              corr3D =  corr3D*4./3.*hdtx*hdty

              Z0m(1) = Z0m(1) + corr3D
              Z0p(1) = Z0p(1) + corr3D


              ! (2) Right state of left node ---------------------------------
              ZLp=Sp(spi,ix,iy,iz-1)+(0.5-max(U(VELZ,ix,iy,iz-1),0.)*hdtz)*delbarZL

              ! x-flux
              sig = signum(U(VELX,ix,iy,iz-1))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix-2:ix+2,iy,iz-1),hy_transOrder,hdtx,sig,transFlux)
              ZLp = ZLp + transFlux*absig*U(VELX,ix,iy,iz-1)

              ! y-flux
              sig = signum(U(VELY,ix,iy,iz-1))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix,iy-2:iy+2,iz-1),hy_transOrder,hdty,sig,transFlux)
              ZLp = ZLp + transFlux*absig*U(VELY,ix,iy,iz-1)

              ! x-y cross deriv at iz-1
              call compute2dCrossDerivative(U(VELX,ix,iy,iz-1),U(VELY,ix,iy,iz-1),Sp(spi,ix-1:ix+1,iy-1:iy+1,iz-1),corr3D)
              corr3D =  corr3D*4./3.*hdtx*hdty

              ZLp = ZLp + corr3D


              !! Upwinding
              sig = signum(HL)
              ZLh(spi) = 0.5*(1.+sig)*ZLp+0.5*(1.-sig)*Z0m(1)


              ! (4) Left state of right node ---------------------------------
              ZRm=Sp(spi,ix,iy,iz+1)-(0.5+min(U(VELZ,ix,iy,iz+1),0.)*hdtz)*delbarZR

              ! x-flux
              sig = signum(U(VELX,ix,iy,iz+1))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix-2:ix+2,iy,iz+1),hy_transOrder,hdtx,sig,transFlux)
              ZRm = ZRm + transFlux*absig*U(VELX,ix,iy,iz+1)

              ! y-flux
              sig = signum(U(VELY,ix,iy,iz+1))
              absig = abs(sig)

              call computeTransFlux(Sp(spi,ix,iy-2:iy+2,iz+1),hy_transOrder,hdty,sig,transFlux)
              ZRm = ZRm + transFlux*absig*U(VELY,ix,iy,iz+1)

              ! x-y cross deriv at iz+1
              call compute2dCrossDerivative(U(VELX,ix,iy,iz+1),U(VELY,ix,iy,iz+1),Sp(spi,ix-1:ix+1,iy-1:iy+1,iz+1),corr3D)
              corr3D =  corr3D*4./3.*hdtx*hdty

              ZRm = ZRm + corr3D

              !! Upwinding
              sig = signum(HR)
              ZRh(spi) = 0.5*(1.+sig)*Z0p(1)+0.5*(1.-sig)*ZRm

           else  !! PPM TYPE 3RD ORDER METHOD==================================================

              do iShift = -1,1
                 !! Use MC slope limiter
                 delbarZL=mc(Sp(spi,ix,iy,iz+iShift-1)-Sp(spi,ix,iy,iz+iShift-2),&
                             Sp(spi,ix,iy,iz+iShift  )-Sp(spi,ix,iy,iz+iShift-1))

                 delbarZ0=mc(Sp(spi,ix,iy,iz+iShift  )-Sp(spi,ix,iy,iz+iShift-1),&
                             Sp(spi,ix,iy,iz+iShift+1)-Sp(spi,ix,iy,  iz+iShift))

                 delbarZR=mc(Sp(spi,ix,iy,iz+iShift+1)-Sp(spi,ix,iy,iz+iShift  ),&
                             Sp(spi,ix,iy,iz+iShift+2)-Sp(spi,ix,iy,iz+iShift+1))


                 !! PPM polynomial interpolation 
                 ! (1) at i-1/2: Left state of center node (i,j,k)
                 ZL=0.5*(Sp(spi,ix,iy,iz+iShift)+Sp(spi,ix,iy,iz+iShift-1))-(delbarZ0-delbarZL)/6.

                 ! (2) at i+1/2: Right state of center node (i,j,k)
                 ZR=0.5*(Sp(spi,ix,iy,iz+iShift)+Sp(spi,ix,iy,iz+iShift+1))-(delbarZR-delbarZ0)/6.

                 !! Apply the same logic for steepening mass scalar/species 
                 !! as for the contact steeping
                 !! NOTE: THE FOLLOWING LINES ARE COMMENTED OUT, BUT CAN BE USED WHEN APPLYING STEEPENING.
                 temp1 = U(DENS,ix,iy,iz+iShift+1) - U(DENS,ix,iy,iz+iShift+1)
!!$                 if (hy_ContactSteepening .and. abs(temp1) > hy_tiny) then
!!$                    ! Eqn 1.17 : Second derivatives
!!$                    del2rhoR = (    U(DENS,ix,iy,iz+iShift+2)&
!!$                                -2.*U(DENS,ix,iy,iz+iShift+1)&
!!$                                  + U(DENS,ix,iy,iz+iShift  ))/(6.*dz**2)
!!$
!!$                    del2rhoL = (    U(DENS,ix,iy,iz+iShift  )&
!!$                                -2.*U(DENS,ix,iy,iz+iShift-1)&
!!$                                   +U(DENS,ix,iy,iz+iShift-2))/(6.*dz**2)
!!$
!!$                    ! Third derivative
!!$                    eta_steep = (del2rhoL-del2rhoR)*dz**2/temp1
!!$                    if (del2rhoR*del2rhoL >= 0.) then
!!$                       eta_steep = 0.
!!$                    endif
!!$                    if (epsln*min(U(DENS,ix,iy,iz+iShift+1),U(DENS,ix,iy,iz+iShift-1))&
!!$                         -abs(U(DENS,ix,iy,iz+iShift+1) - U(DENS,ix,iy,iz+iShift-1)) >= 0.) then
!!$                       eta_steep = 0.
!!$                    endif
!!$
!!$                    ! Eqn 1.16
!!$                    eta_steep = max(0., min(1., eta1*(eta_steep - eta2)))
!!$
!!$                    ! Eqn 3.2
!!$                    temp2 = abs(U(PRES,ix,iy,iz+iShift+1)-U(PRES,ix,iy,iz+iShift-1))&
!!$                    /min(U(PRES,ix,iy,iz+iShift+1),U(PRES,ix,iy,iz+iShift-1))
!!$
!!$                    temp3 = abs(U(DENS,ix,iy,iz+iShift+1)-U(DENS,ix,iy,iz+iShift-1))&
!!$                         /min(U(DENS,ix,iy,iz+iShift+1),U(DENS,ix,iy,iz+iShift-1))
!!$
!!$                    if (U(GAME,ix,iy,iz+iShift)*K0*temp3-temp2 < 0.0) then
!!$                       eta_steep = 0.
!!$                    endif
!!$
!!$                    ! Eqn 1.15
!!$                    ZL = ZL*(1.-eta_steep) + (Sp(spi,ix,iy,iz+iShift-1)+0.5*delbarZL)*eta_steep
!!$                    ZR = ZR*(1.-eta_steep) + (Sp(spi,ix,iy,iz+iShift+1)-0.5*delbarZR)*eta_steep
!!$                 endif

                 !! Apply limiting
                 ZL=max(min(Sp(spi,ix,iy,iz+iShift),Sp(spi,ix,iy,iz+iShift-1)),&
                      min(max(Sp(spi,ix,iy,iz+iShift),Sp(spi,ix,iy,iz+iShift-1)),ZL))

                 !! Apply limiting
                 ZR=max(min(Sp(spi,ix,iy,iz+iShift),Sp(spi,ix,iy,iz+iShift+1)),&
                      min(max(Sp(spi,ix,iy,iz+iShift),Sp(spi,ix,iy,iz+iShift+1)),ZR))

                 !! Check monotonicity
                 if ((ZR-Sp(spi,ix,iy,iz+iShift))*(Sp(spi,ix,iy,iz+iShift)-ZL)<=0.) then
                    ZR=Sp(spi,ix,iy,iz+iShift)
                    ZL=Sp(spi,ix,iy,iz+iShift)
                 endif
                 if (6.*(ZR-ZL)*(Sp(spi,ix,iy,iz+iShift)-0.5*(ZL+ZR)) > (ZR-ZL)**2) then
                    ZL = 3.*Sp(spi,ix,iy,iz+iShift)-2.*ZR
                 endif
                 if (6.*(ZR-ZL)*(Sp(spi,ix,iy,iz+iShift)-0.5*(ZL+ZR)) <-(ZR-ZL)**2) then
                    ZR = 3.*Sp(spi,ix,iy,iz+iShift)-2.*ZL
                 endif

                 !! PPM reconstruction to advance 1/2 time step (i.e., n+1/2) at cell interfaces
                 delZ=ZR-ZL
                 Z6=6.*(Sp(spi,ix,iy,iz+iShift)-0.5*(ZR+ZL))

                 Z0m(iShift)=ZL-min(U(VELZ,ix,iy,iz+iShift),0.)*hdtz*&
                      (delZ+(1.0+min(U(VELZ,ix,iy,iz+iShift),0.)*4./3.*hdtz)*Z6)

                 Z0p(iShift)=ZR-max(U(VELZ,ix,iy,iz+iShift),0.)*hdtz*&
                      (delZ-(1.0-max(U(VELZ,ix,iy,iz+iShift),0.)*4./3.*hdtz)*Z6)



                 !! Transverse flux components correction-------------------------------------
                 ! x-flux
                 sig = signum(U(VELX,ix,iy,iz+iShift))
                 absig = abs(sig)

                 call computeTransFlux(Sp(spi,ix-2:ix+2,iy,iz+iShift),hy_transOrder,hdtx,sig,transFlux)
                 Z0m(iShift) = Z0m(iShift) + transFlux*absig*U(VELX,ix,iy,iz+iShift)
                 Z0p(iShift) = Z0p(iShift) + transFlux*absig*U(VELX,ix,iy,iz+iShift)

                 ! y-flux
                 sig = signum(U(VELY,ix,iy,iz+iShift))
                 absig = abs(sig)

                 call computeTransFlux(Sp(spi,ix,iy-2:iy+2,iz+iShift),hy_transOrder,hdty,sig,transFlux)
                 Z0m(iShift) = Z0m(iShift) + transFlux*absig*U(VELY,ix,iy,iz+iShift)
                 Z0p(iShift) = Z0p(iShift) + transFlux*absig*U(VELY,ix,iy,iz+iShift)


                 ! x-y cross deriv at iz+iShift
                 call compute2dCrossDerivative(U(VELX,ix,iy,iz+iShift),U(VELY,ix,iy,iz+iShift),&
                                               Sp(spi,ix-1:ix+1,iy-1:iy+1,iz+iShift),corr3D)
                 corr3D = corr3D*4./3.*hdtx*hdty

                 Z0m(iShift) = Z0m(iShift) + corr3D
                 Z0p(iShift) = Z0p(iShift) + corr3D


              enddo

              !! Updates -------------------------------------------------------------------
              !! Upwinding at i-1/2
              sig = signum(HL)
              ZLh(spi) = 0.5*(1.+sig)*Z0p(-1)+0.5*(1.-sig)*Z0m(0)

              !! Upwinding at i+1/2
              sig = signum(HR)
              ZRh(spi)= 0.5*(1.+sig)*Z0p(0)+0.5*(1.-sig)*Z0m(1)

           endif

           if (spn <= SPECIES_END) then
              sumSpeciesZL = sumSpeciesZL + ZLh(spi)
              sumSpeciesZR = sumSpeciesZR + ZRh(spi)
           endif

        endif ! end of NDIM == 3
     endif ! end of NDIM >= 2
  enddo


  !! Initialize summation for mass fractions
  sumSpecies = 0.

  !! Now update to n+1
  do spn = SPECIES_BEGIN, MASS_SCALARS_END
     spi = spn-NPROP_VARS

     !! Convserving fluxes of mass fractions
     if (spn <= SPECIES_END) then
        XRh(spi) = XRh(spi)/sumSpeciesXR
        XLh(spi) = XLh(spi)/sumSpeciesXL
        if (NDIM >= 2) then
           YRh(spi) = YRh(spi)/sumSpeciesYR
           YLh(spi) = YLh(spi)/sumSpeciesYL
           if (NDIM == 3) then
              ZRh(spi) = ZRh(spi)/sumSpeciesZR
              ZLh(spi) = ZLh(spi)/sumSpeciesZL
           endif
        endif
     endif

#if NDIM == 1
        SpNew(spi) =  Sp(spi,ix,iy,iz)*U(DENS,ix,iy,iz) &
             - dt*(XRh(spi)*FR-XLh(spi)*FL)/dx

#elif NDIM == 2
        SpNew(spi) =  Sp(spi,ix,iy,iz)*U(DENS,ix,iy,iz) &
             - dt*((XRh(spi)*FR-XLh(spi)*FL)/dx &
                  +(YRh(spi)*GR-YLh(spi)*GL)/dy)

#elif NDIM == 3
        SpNew(spi) =  Sp(spi,ix,iy,iz)*U(DENS,ix,iy,iz) &
             - dt*((XRh(spi)*FR-XLh(spi)*FL)/dx &
                  +(YRh(spi)*GR-YLh(spi)*GL)/dy &
                  +(ZRh(spi)*HR-ZLh(spi)*HL)/dz)
#endif

     if (densNew .NE. 0.0) then

        SpNew(spi) = SpNew(spi)/densNew
        
        !! Conserving mass fractions
        if (spn <= SPECIES_END) then
           sumSpecies = SpNew(spi) + sumSpecies
        endif

     end if

  enddo

  if (densNew .NE. 0.0) then
     !! Conserving mass fractions
     do spn = SPECIES_BEGIN, SPECIES_END
        spi = spn-NPROP_VARS
        if (spn <= SPECIES_END) then
           SpNew(spi) = SpNew(spi)/sumSpecies
        endif
     enddo
  end if

end Subroutine hy_uhd_updateSpeciesMassScalar




Subroutine computeTransFlux(FluxArray,transOrder,hdtx,sign,transFlux)

  implicit none
  real, dimension(-2:2), intent(IN) :: FluxArray
  integer, intent(IN) :: transOrder
  real, intent(IN)    :: hdtx,sign
  real, intent(OUT)   :: transFlux

  real   :: transPos, transNeg, secondDeriv, num, denom

  select case(transOrder)
  case(0)
     transPos=0.
     transNeg=0.
  case(1)
     transPos = FluxArray( 0)-FluxArray(-1)
     transNeg = FluxArray( 1)-FluxArray( 0)
  case(3)
     transPos = (2.*FluxArray( 1)+3.*FluxArray( 0)-6.*FluxArray(-1)   +FluxArray(-2))/6. 
     transNeg = (  -FluxArray( 2)+6.*FluxArray( 1)-3.*FluxArray( 0)-2.*FluxArray(-1))/6.
  end select

   transFlux =  - hdtx*(transPos*(1.+sign) + transNeg*(1.-sign))*0.5

end Subroutine computeTransFlux


Subroutine compute2dCrossDerivative(u1,u2,vec,corr3D)

  implicit none
  real, intent(IN) :: u1, u2
  real, intent(IN), dimension(-1:1,-1:1) :: vec
  real, intent(OUT) :: corr3D

  if     ((u1>0.) .and. (u2>0.)) then

     corr3D = (vec( 0, 0)-vec(-1, 0)) &
             -(vec( 0,-1)-vec(-1,-1))

  elseif ((u1>0.) .and. (u2<0.)) then

     corr3D = (vec( 0, 1)-vec(-1, 1)) &
             -(vec( 0, 0)-vec(-1, 0))


  elseif ((u1<0.) .and. (u2>0.)) then

     corr3D = (vec( 1, 0)-vec( 0, 0)) &
             -(vec( 1,-1)-vec( 0,-1))


  elseif ((u1<0.) .and. (u2<0.)) then

     corr3D = (vec( 1, 1)-vec( 0, 1)) &
             -(vec( 1, 0)-vec( 0, 0))


  elseif ((u1==0.) .and. (u2>0.)) then

     corr3D = (vec( 1, 0)-vec(-1, 0)) &
             -(vec( 1,-1)-vec(-1,-1))

     corr3D = 0.5*corr3D

  elseif ((u1==0.) .and. (u2<0.)) then

     corr3D = (vec( 1, 1)-vec(-1, 1)) &
             -(vec( 1, 0)-vec(-1, 0))

     corr3D = 0.5*corr3D

  elseif ((u1>0.) .and. (u2==0.)) then

     corr3D = (vec( 0, 1)-vec(-1, 1)) &
             -(vec( 0,-1)-vec(-1,-1))

     corr3D = 0.5*corr3D

  elseif ((u1<0.) .and. (u2==0.)) then

     corr3D = (vec( 1, 1)-vec( 0, 1)) &
             -(vec( 1,-1)-vec( 0,-1))

     corr3D = 0.5*corr3D

  elseif ((u1==0.) .and. (u2==0.)) then

     corr3D = (vec( 1, 1)-vec(-1, 1)) &
             -(vec( 1,-1)-vec(-1,-1))

     corr3D = 0.25*corr3D

  endif

  corr3D = corr3D*u1*u2

End Subroutine compute2dCrossDerivative
