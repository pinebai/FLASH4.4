!!****if* source/physics/IncompNS/IncompNSMain/constdens/ins_vt
!!
!! NAME
!!
!!  ins_vt
!!
!! SYNOPSIS
!!
!!  call ins_vt(:: isgs,
!!               :: ng,
!!               :: nxc,
!!               :: nyc,
!!               :: nzc,
!!               :: ru1,
!!               :: dx,
!!               :: dy,
!!               :: dz,
!!               :: coord,
!!               :: bsize,
!!              real, pointer, dimension(:,:,:,:)  :: facexdata,
!!              real, pointer, dimension(:,:,:,:)  :: faceydata,
!!              real, pointer, dimension(:,:,:,:)  :: facezdata,
!!              real, pointer, dimension(:,:,:,:)  :: solndata)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   isgs : 
!!
!!   ng : 
!!
!!   nxc : 
!!
!!   nyc : 
!!
!!   nzc : 
!!
!!   ru1 : 
!!
!!   dx : 
!!
!!   dy : 
!!
!!   dz : 
!!
!!   coord : 
!!
!!   bsize : 
!!
!!   facexdata : 
!!
!!   faceydata : 
!!
!!   facezdata : 
!!
!!   solndata : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!!$c
!!$c-----------------------------------------------------------------------
!!$c                 ***************************                         
!!$c                 *         turvis.f        *                        
!!$c                 ***************************                       
!!$c----------------------------------------------------------------------- 
!!$c
!!$c    - Turvis:       calls the different eddy viscosity routines
!!$c
!!$c----------------------------------------------------------------------- 
!!$c
!!$c
!!$c-----SUBROUTINE-Turvis------------------------E. Balaras 7/12/98-------
!!$c
!!$c-----Adapted to AMR: M. Vanella, June 2007.----------------------------
!!$c
!!$c
  SUBROUTINE ins_vt(isgs,ng,nxc,nyc,nzc,RU1,dx,dy,dz,   &
                    coord,bsize,&
                    facexData,&
                    faceyData,&
                    facezData,&
                    solnData)             

  use ins_interface, only : ins_velgradtensor


  implicit none
 
#include "constants.h"
#include "Flash.h"


  integer isgs,ng,nxc,nyc,nzc
  real  RU1,dx,dy,dz       
  real coord(3),bsize(3)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  ! Local variables:
  real  dxdydz
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: UC,VC,WC
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: SXX,SYY,SZZ,&
        SXY,SYZ,SXZ,tpdvdxc,tpdwdxc,G

  integer i,j,k

  real ALPHA,RPLUS

  real ywall,ypos,dyaux
  real, parameter :: css = 0.16
  integer nxi,nyj,nzk

!!$c
!!$c----------------------------------------------------------------------
!!$c                      square of the ratio between test and grid filter
!!$c----------------------------------------------------------------------

  nxi = NXB + ng
  nyj = NYB + ng
  nzk = NZB + ng


  dxdydz = (dx*dy*dz)**(0.33333333333333333)
!!$c
!!$c----------------------------------------------------------------------
!!$c                       calculate Strain Rates and |S| at cell centers
!!$c----------------------------------------------------------------------

  call ins_velgradtensor(ng,facexData,faceyData,facezData, &              
                         dx,dy,dz,SXX,SXY,SXZ,             &
                         tpdvdxc,SYY,SYZ,tpdwdxc,G,SZZ)

  SXY = 0.5*(SXY + tpdvdxc)
  SXZ = 0.5*(SXZ + tpdwdxc)
  SYZ = 0.5*(SYZ + G)

  G = SQRT(2.*(SXX**2+SYY**2+SZZ**2)+4.*(SXY**2+SYZ**2+SXZ**2)) 


  ! Smagorinsky model
  IF(ISGS==1) THEN

     DO j= ng+1 , nyj

        solnData(TVIS_VAR,ng+1:nxi,j,ng+1:nzk) = &
        (CSS*DXDYDZ)**2*G(ng+1:nxi,j,ng+1:nzk)

     ENDDO

  endif



  RETURN

END SUBROUTINE ins_VT


