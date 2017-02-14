
! This routine fixes inconsistent fluxes across block interfaces 
! of different refinement levels as defined by the normal 
! velocity component.
!
! Adapted to Flash 3, Marcos Vanella, September 2007.


  SUBROUTINE ins_fluxfix(ng,nxc,nyc,nzc,nxi,nyj,nzk,blockCount,&
                         blockList)

  
  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes
  use IncompNS_data, ONLY : ins_meshMe

  implicit none
  
#include "constants.h"
#include "Flash.h"  
#include "IncompNS.h"

  integer, intent(IN) :: ng,nxc,nyc,nzc,nxi,nyj,nzk, &
                         blockCount     
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList


  integer, dimension(MDIM) :: datasize
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData

  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_u
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_v
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_w

  integer :: sx,sy,sz,ex,ey,ez
  integer :: lb,blockID,level
  real :: del(MDIM)
  
  flxint_u = 0.
  flxint_v = 0.
  flxint_w = 0.

  level = -1

  ! Dump fluxes into flux variables:   
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take cell depth to be 1:
     del(DIR_Z) = 1. 
#endif

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   

     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(1)-NGUARD
     ey = dataSize(2)-NGUARD*K2D
     ez = dataSize(3)-NGUARD*K3D

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif
     ! copy interface velocities to flux variables
     ! X direction:
     flxint_u(VELC_FLUX,sx,sy:ey,sz:ez) = &
     facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez)*del(DIR_Y)*del(DIR_Z)
     flxint_u(VELC_FLUX,ex+1,sy:ey,sz:ez) = &
     facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez)*del(DIR_Y)*del(DIR_Z)

     call Grid_putFluxData(blockID, IAXIS, &
                           flxint_u, dataSize)             ! dy*dz


     ! Y direction:
     flxint_v(VELC_FLUX,sx:ex,sy,sz:ez) = &
     faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez)*del(DIR_X)*del(DIR_Z)
     flxint_v(VELC_FLUX,sx:ex,ey+1,sz:ez) = &
     faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez)*del(DIR_X)*del(DIR_Z)

     call Grid_putFluxData(blockID, JAXIS, &
                           flxint_v, dataSize)             ! dx*dz     


#if NDIM == 3
     ! Z direction:
     flxint_w(VELC_FLUX,sx:ex,sy:ey,sz) = &
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz)*del(DIR_X)*del(DIR_Y)
     flxint_w(VELC_FLUX,sx:ex,sy:ey,ez+1) = &
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1)*del(DIR_X)*del(DIR_Y)

     call Grid_putFluxData(blockID, KAXIS, &
                           flxint_w, dataSize)             ! dx*dy 
#endif


     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo

!  write(*,*) 'Mype=',Ins_MeshMe, 'Before Grid Conserve fluxes'

  ! fix fluxes at refinement interfaces
  call Grid_conserveFluxes( ALLDIR, level)

!  write(*,*) 'Mype=',Ins_MeshMe, 'After Grid Conserve fluxes'

 
  ! Get fixed fluxes:
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take cell depth to be 1:
     del(DIR_Z) = 1. 
#endif

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   


     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(1)-NGUARD
     ey = dataSize(2)-NGUARD*K2D
     ez = dataSize(3)-NGUARD*K3D


     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     ! update interface velocities with flux variables
     ! X direction:
     call Grid_getFluxData(blockID, IAXIS, &
     flxint_u, dataSize)            

     facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez) = &
       flxint_u(VELC_FLUX,sx,sy:ey,sz:ez)/(del(DIR_Y)*del(DIR_Z))  ! dy*dz
     facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez) = &
       flxint_u(VELC_FLUX,ex+1,sy:ey,sz:ez)/(del(DIR_Y)*del(DIR_Z))


     ! Y direction:
     call Grid_getFluxData(blockID, JAXIS, &
     flxint_v, dataSize)                 

     faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez) = &
       flxint_v(VELC_FLUX,sx:ex,sy,sz:ez)/(del(DIR_X)*del(DIR_Z))  ! dx*dz
     faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez) = &
       flxint_v(VELC_FLUX,sx:ex,ey+1,sz:ez)/(del(DIR_X)*del(DIR_Z))

#if NDIM == 3
     ! Z direction:
     call Grid_getFluxData(blockID, KAXIS, &
     flxint_w, dataSize)              

     facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz) = &
       flxint_w(VELC_FLUX,sx:ex,sy:ey,sz)/(del(DIR_X)*del(DIR_Y))  ! dx*dy
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1) = &
       flxint_w(VELC_FLUX,sx:ex,sy:ey,ez+1)/(del(DIR_X)*del(DIR_Y))
#endif

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo


  return

END SUBROUTINE ins_fluxfix

!##########################################################

! This routine fixes inconsistent fluxes across block interfaces 
! of different refinement levels as defined by the pressure
! gradient.

  SUBROUTINE ins_fluxfix_p(ng,nxc,nyc,nzc,nxi,nyj,nzk,pvar,&
                         blockCount,blockList)

  use Grid_interface, ONLY : Grid_getDeltas,         &                             
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes

  implicit none
  
#include "constants.h"
#include "Flash.h"  
#include "IncompNS.h"

  integer, intent(in) :: ng,nxc,nyc,nzc,nxi,nyj,nzk,pvar,&
                         blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList


  integer, dimension(MDIM) :: datasize
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  
  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData

  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_u
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_v
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_w

  integer :: sx,sy,sz,ex,ey,ez
  integer :: lb,blockID,level
  real :: del(MDIM)  

  ! Dump fluxes into flux variables:   
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take cell depth to be 1:
     del(DIR_Z) = 1. 
#endif

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   

     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(1)-NGUARD
     ey = dataSize(2)-NGUARD*K2D
     ez = dataSize(3)-NGUARD*K3D

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     ! copy interface velocities to flux variables
     ! X direction:
     flxint_u(VELC_FLUX,sx,sy:ey,sz:ez) =    &
     (solnData(pvar,sx,sy:ey,sz:ez)- &
      solnData(pvar,sx-1,sy:ey,sz:ez))*del(DIR_Y)*del(DIR_Z)/del(DIR_X)

     flxint_u(VELC_FLUX,ex+1,sy:ey,sz:ez) =    &
     (solnData(pvar,ex+1,sy:ey,sz:ez)- &
      solnData(pvar,ex,sy:ey,sz:ez))*del(DIR_Y)*del(DIR_Z)/del(DIR_X)

     call Grid_putFluxData(blockID, IAXIS, &
                           flxint_u, dataSize)             ! dy*dz


     ! Y direction:
     flxint_v(VELC_FLUX,sx:ex,sy,sz:ez) =    &
     (solnData(pvar,sx:ex,sy,sz:ez)- &
      solnData(pvar,sx:ex,sy-1,sz:ez))*del(DIR_X)*del(DIR_Z)/del(DIR_Y)

     flxint_v(VELC_FLUX,sx:ex,ey+1,sz:ez) =    &
     (solnData(pvar,sx:ex,ey+1,sz:ez)- &
      solnData(pvar,sx:ex,ey,sz:ez))*del(DIR_X)*del(DIR_Z)/del(DIR_Y)

     call Grid_putFluxData(blockID, JAXIS, &
                           flxint_v, dataSize)             ! dx*dz     

#if NDIM == 3
     ! Z direction:
     flxint_w(VELC_FLUX,sx:ex,sy:ey,sz) =    &
     (solnData(pvar,sx:ex,sy:ey,sz)- &
      solnData(pvar,sx:ex,sy:ey,sz-1))*del(DIR_X)*del(DIR_Y)/del(DIR_Z)

     flxint_w(VELC_FLUX,sx:ex,sy:ey,ez+1) =    &
     (solnData(pvar,sx:ex,sy:ey,ez+1)- &
      solnData(pvar,sx:ex,sy:ey,ez))*del(DIR_X)*del(DIR_Y)/del(DIR_Z)


     call Grid_putFluxData(blockID, KAXIS, &
                           flxint_w, dataSize)             ! dx*dy 
#endif

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  enddo


  ! fix fluxes at refinement interfaces
  call Grid_conserveFluxes( ALLDIR, level)

  return

END SUBROUTINE ins_fluxfix_p

