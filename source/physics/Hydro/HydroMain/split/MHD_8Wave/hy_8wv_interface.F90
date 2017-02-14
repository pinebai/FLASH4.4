!!****ih* source/physics/Hydro/HydroMain/split/MHD_8Wave/hy_8wv_interface
!!
!! NAME
!!  hy_8wv_interface
!!
!! SYNOPSIS
!!  use hy_8wv_interface
!!
!! DESCRIPTION
!!  This is an interface specific for the 8Wave 
!!  MHD module that defines its public interfaces.
!!
!!***
Module hy_8wv_interface

  implicit none

#include "constants.h"
#include "Flash.h"

  interface
     subroutine hy_8wv_divb(blockCount,blockList,dt)
       implicit none
       integer, intent(IN) :: blockCount
       integer, intent(IN), dimension(blockCount) :: blockList
       real, intent(IN)    :: dt
     end subroutine hy_8wv_divb
  end interface


  interface
     subroutine hy_8wv_fluxes(Um,Up,Flux,speed,vint,n,dir)
       implicit none
       integer, INTENT(in) :: n
       real, DIMENSION(NUNK_VARS,n), INTENT(in)  :: Um,Up
       real, DIMENSION(NFLUXES,n), INTENT(out) :: Flux
       real, DIMENSION(n), INTENT(out) :: speed, vint
       integer, INTENT(in):: dir
     end subroutine hy_8wv_fluxes
  end interface


  interface
     subroutine hy_8wv_interpolate(U,Uc,Um,Up,grav,dt,xc,dx,n,dir)
       implicit none
       integer, INTENT(IN) :: n,dir
       real, DIMENSION(NUNK_VARS,n), INTENT(IN) :: U
       real, DIMENSION(NUNK_VARS,n), INTENT(OUT) :: Uc,Um,Up
       real, DIMENSION(n), INTENT(IN) :: grav,xc,dx
       real, INTENT(IN) :: dt
     end subroutine hy_8wv_interpolate
  end interface


  interface
     subroutine hy_8wv_setTimestep(hy_meshMe,dtmin,i,j,k,blockID)
       implicit none
       real, INTENT(in) :: dtmin
       integer, INTENT(in) :: i,j,k,blockID,hy_meshMe
     end subroutine hy_8wv_setTimestep
  end interface


  interface
     subroutine hy_8wv_sources(Uc,Um,Up,S,grav,xc,dx,n,dir)
       implicit none
       integer, INTENT(IN) :: n,dir
       real, DIMENSION(NUNK_VARS,n), INTENT(OUT) :: S
       real, DIMENSION(NUNK_VARS,n), INTENT(IN) :: Uc,Um,Up
       real, DIMENSION(n), INTENT(IN) :: grav,xc,dx
       integer :: i,VELN_VAR,MAGN_VAR
     end subroutine hy_8wv_sources
  end interface


  interface
     subroutine hy_8wv_sweep( blockCount, blockList, dt, sweepDir)
       implicit none
       integer, intent(IN) :: sweepDir, blockCount
       integer, intent(IN), dimension(blockCount) :: blockList
       real, intent(inout) :: dt
     end subroutine hy_8wv_sweep
  end interface


  interface
     subroutine hy_8wv_addResistiveFluxes(i1,i2,B,Flux,visc,x,y,z,nx,ny,nz,sweepDir)
       implicit none
       integer, INTENT(IN) :: i1,i2,nx,ny,nz,sweepDir
       real, DIMENSION(3,nx,ny,nz), INTENT(IN) :: B
       real, DIMENSION(NFLUXES,max(nx,ny,nz)), INTENT(OUT) :: Flux
       real, DIMENSION(max(nx,ny,nz)), INTENT(IN) :: x,y,z,visc
     end subroutine hy_8wv_addResistiveFluxes
  end interface


  interface
     subroutine hy_8wv_addThermalFluxes(i1,i2,T,Flux,cond,x,y,z,nx,ny,nz,sweepDir)
       implicit none
       integer, INTENT(IN) :: i1,i2,nx,ny,nz,sweepDir
       real, DIMENSION(nx,ny,nz), INTENT(IN) :: T
       real, DIMENSION(NFLUXES,max(nx,ny,nz)), INTENT(OUT) :: Flux
       real, DIMENSION(max(nx,ny,nz)), INTENT(IN) :: x,y,z,cond
     end subroutine hy_8wv_addThermalFluxes
  end interface


  interface
     subroutine hy_8wv_addViscousFluxes(i1,i2,V,Flux,visc,x,y,z,nx,ny,nz,sweepDir)
       implicit none
       integer, INTENT(IN) :: i1,i2,nx,ny,nz,sweepDir
       real, DIMENSION(3,nx,ny,nz), INTENT(IN) :: V
       real, DIMENSION(NFLUXES,max(nx,ny,nz)), INTENT(OUT) :: Flux
       real, DIMENSION(max(nx,ny,nz)), INTENT(IN) :: x,y,z,visc
     end subroutine hy_8wv_addViscousFluxes
  end interface

end Module hy_8wv_interface
