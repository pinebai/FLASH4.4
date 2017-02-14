!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_getFluxDeriv
!!
!! NAME
!!
!!  hy_uhd_getFluxDeriv
!!
!! SYNOPSIS
!!
!!  hy_uhd_getFluxDeriv ( integer(IN)  :: ix,
!!                        integer(IN)  :: iy,
!!                        integer(IN)  :: iz,
!!                        integer(IN)  :: blkLimitsGC(2,MDIM)
!!                        integer(IN)  :: fluxType,
!!                        integer(IN)  :: DerivDir,
!!                        real(IN)     :: faceFlux(:,:,:,:),
!!                        real(OUT)    :: Flux1Deriv,
!!                        real(OUT)    :: Flux2Deriv)  
!!
!! ARGUMENTS
!!
!!  ix,iy,iz    - local indices where the flux derivatives are computed
!!  blkLimitsGC - the start and end indices for the block with gcells
!!  fluxType    - flux in x,y,z directions
!!  DerivDir    - a direction of derivative
!!  faceFlux    - face flux
!!  Flux1Deriv  - 1st derivative of the given faceFlux
!!  Flux2Deriv  - 2nd derivative of the given faceFlux
!!
!!
!! DESCRIPTION
!!
!!  This routine calculates 1st and 2nd derivatives of the faceFlux that are used to compute
!!  the modified electric fields construction.
!!
!! REFERENCES
!!
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!
!!*** 

!!REORDER(4): faceFlux

Subroutine hy_uhd_getFluxDeriv( ix,iy,iz,blkLimitsGC,&
                                fluxType,DerivDir,   &
                                faceFlux,            &
                                Flux1Deriv,Flux2Deriv)

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration ------------------------
  integer, intent(IN) :: ix,iy,iz
  integer, intent(IN) :: fluxType,DerivDir
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
#ifdef FIXEDBLOCKSIZE
  real, intent(IN) :: faceFlux(NFLUXES, &
                               GRID_ILO_GC:GRID_IHI_GC,  &
                               GRID_JLO_GC:GRID_JHI_GC,  &
                               GRID_KLO_GC:GRID_KHI_GC)
#else
  real, intent(IN) :: faceFlux(NFLUXES,&
                               blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) 
#endif
  real, intent(OUT):: Flux1Deriv,Flux2Deriv
  !! ---------------------------------------------------

  select case(fluxType)
  case(DIR_X)
     !! [Flux_x(F07MAGY_FLUX) = -Ez, Flux_x(F08MAGZ_FLUX) = Ey]
     select case(DerivDir)
     case(DIR_Y)
        !! dyFlux_x, dy2Flux_x (needed for Ez)
        Flux1Deriv=.5*(-faceFlux(F07MAGY_FLUX,ix,iy+1,iz)&
                       +faceFlux(F07MAGY_FLUX,ix,iy-1,iz))

        Flux2Deriv=    -faceFlux(F07MAGY_FLUX,ix,iy+1,iz)&
                    +2.*faceFlux(F07MAGY_FLUX,ix,iy,  iz)&
                       -faceFlux(F07MAGY_FLUX,ix,iy-1,iz)
     case(DIR_Z)
        !! dzFlux_x, dz2Flux_x (needed for Ey)
        Flux1Deriv=.5*(faceFlux(F08MAGZ_FLUX,ix,iy,iz+1)&
                      -faceFlux(F08MAGZ_FLUX,ix,iy,iz-1))

        Flux2Deriv=    faceFlux(F08MAGZ_FLUX,ix,iy,iz+1)&
                   -2.*faceFlux(F08MAGZ_FLUX,ix,iy,iz  )&
                      +faceFlux(F08MAGZ_FLUX,ix,iy,iz-1)
     end select

  case (DIR_Y)
     !! [Flux_y(F06MAGX_FLUX) = Ez, Flux_y(F08MAGZ_FLUX) = -Ex]
     select case(DerivDir)
     case(DIR_X)
        !! dxFlux_y, dx2Flux_y (needed for Ez)
        Flux1Deriv=.5*(faceFlux(F06MAGX_FLUX,ix+1,iy,iz)&
                      -faceFlux(F06MAGX_FLUX,ix-1,iy,iz))

        Flux2Deriv=    faceFlux(F06MAGX_FLUX,ix+1,iy,iz)&
                   -2.*faceFlux(F06MAGX_FLUX,ix,  iy,iz)&
                      +faceFlux(F06MAGX_FLUX,ix-1,iy,iz)
     case(DIR_Z)
        !! dzFlux_y, dz2Flux_y (needed for Ex)
        Flux1Deriv=.5*(-faceFlux(F08MAGZ_FLUX,ix,iy,iz+1)&
                       +faceFlux(F08MAGZ_FLUX,ix,iy,iz-1))

        Flux2Deriv=    -faceFlux(F08MAGZ_FLUX,ix,iy,iz+1)&
                    +2.*faceFlux(F08MAGZ_FLUX,ix,iy,iz  )&
                       -faceFlux(F08MAGZ_FLUX,ix,iy,iz-1)
     end select

  case(DIR_Z)
     !! [Flux_z(F06MAGX_FLUX) = -Ey, Flux_z(F07MAGY_FLUX) = Ex]
     select case(DerivDir)
     case(DIR_X)
        !! dxFlux_z, dx2Flux_z (needed for Ey)
        Flux1Deriv=.5*(-faceFlux(F06MAGX_FLUX,ix+1,iy,iz)&
                       +faceFlux(F06MAGX_FLUX,ix-1,iy,iz))
        Flux2Deriv=    -faceFlux(F06MAGX_FLUX,ix+1,iy,iz)&
                    +2.*faceFlux(F06MAGX_FLUX,ix,  iy,iz)&
                       -faceFlux(F06MAGX_FLUX,ix-1,iy,iz)
     case(DIR_Y)
        !! dyFlux_z, dy2Flux_z (needed for Ex)
        Flux1Deriv=.5*(faceFlux(F07MAGY_FLUX,ix,iy+1,iz)&
                      -faceFlux(F07MAGY_FLUX,ix,iy-1,iz))

        Flux2Deriv=    faceFlux(F07MAGY_FLUX,ix,iy+1,iz)&
                   -2.*faceFlux(F07MAGY_FLUX,ix,iy,  iz)&
                      +faceFlux(F07MAGY_FLUX,ix,iy-1,iz)
     end select

  end select

End Subroutine hy_uhd_getFluxDeriv
