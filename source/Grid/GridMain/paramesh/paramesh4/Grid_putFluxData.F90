!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_putFluxData
!!
!! NAME
!!  Grid_putFluxData
!!
!! SYNOPSIS
!!
!!
!!  call Grid_putFluxData(integer(IN) :: blockID,
!!                   integer(IN) :: axis,
!!                   real(IN)    :: fluxes(NFLUXES,dataSize(1),dataSize(2),dataSize(3)),
!!                   integer(IN) :: dataSize(3),
!!          OPTIONAL,integer(IN) :: pressureSlots(:),
!!          OPTIONAL,real(IN)    :: areaLeft(:,:,:))
!!
!! DESCRIPTION 
!!
!!
!!  Put the fluxes in a direction specified by axis for boundary cells
!!  for block blockID. This routine needs to be used with adaptive mesh
!!  since fluxes calculated by the two blocks that are at fine/coarse boundary have 
!!  different accuracy. The fluxes calculated by individual blocks are reported to 
!!  the Grid through this call. Once that is done, a call to Grid_conserveFluxes 
!!  applies the flux conservation algorithm to make it consistent across the fine/coarse 
!!  boundaries.
!!
!! ARGUMENTS
!!
!!  blockID : The local blockid
!!
!!
!!  axis : integer value specifying on which cell faces to put fluxes. 
!!         The options are IAXIS, JAXIS, or KAXIS defined in constants.h
!!
!!
!!  fluxes :  real array with space for fluxes, through one axis, 
!!            for all cells of a block and for all flux variables.
!!            fluxes(VAR, i, j, k) is VAR's flux through 
!!            the left cell face for cell i, j, k.
!!
!!
!!  dataSize : integer array specifying the dimensions for the array, fluxes
!!
!!             dataSize (1) holds the number of cells provided in the i direction
!!
!!             dataSize (2) holds the number of cells provided in the j direction
!!                          if 1 d problem, set datasize(2) = 1
!!
!!             dataSize (3) holds the number of cells provided in the k direction
!!                          if 1 or 2 d problem, set datasize(3) = 1
!!
!!             fluxes should contain space for fluxes of all cells in the block, 
!!             including guardcells, and the  fluxes must be correct for 
!!             the interior cells of the block, as this interface does not know which 
!!             cell fluxes the Grid will need to store.
!!
!!  pressureSlots: If present and greater than zero, this indicates one or more flux
!!                 variables in the fluxes array that may need special handling because
!!                 they really scale like flux densities when other flux variables scale
!!                 like fluxes.  For the split PPM Hydro implementation, for example,
!!                 this is normally used for the pressure "flux" variable, but
!!                 it can be applied to other flux variables that the caller keeps in
!!                 flux density form.
!!                 The special handling consists in multiplying the flux variables with
!!                 the appropriate face areas, which are taken from the areaLeft array
!!                 argument, before storing them in the arrays on which the Grid unit
!!                 acts.
!!                 Special handling should only be requested if the Grid units handles
!!                 flux variable "as fluxes" (or else if it does nto matter anyway,
!!                 as is the case in 1D).
!!                 The pressureSlots argument in the corresponding Grid_getFluxData
!!                 call should generally match the one in the Grid_putFluxData call.
!!
!!  areaLeft :     areas of left and right faces, only used if special scaling is
!!                 requested with the pressureSlot argument.
!!                 The areaLeft argument in the corresponding Grid_getFluxData
!!                 call should generally match the one in the Grid_putFluxData call.
!!
!! NOTES 
!!
!!   Any code calling this subroutine needs to know the explicit interface,
!!   since this interface contains optional dummy arguments and assumed-shape
!!   dummy arrays. Calling FORTRAN units should therefore contain a line like
!!       use Grid_interface, ONLY: Grid_putFluxData
!!
!!   This implementation is specific to Paramesh 4.
!!
!! SEE ALSO
!!
!!   Grid_getFluxData
!!   Grid_conserveFluxes
!!   hy_sweep
!!***

!!REORDER(5): flux_[xyz],gr_[xyz]flx
!!REORDER(4): fluxes
!!REORDER(5): gr_xflx_[yz]face, gr_yflx_[xz]face, gr_zflx_[xy]face
#include "Flash.h"
subroutine Grid_putFluxData(blockID, axis, fluxes, dataSize, pressureSlots, areaLeft)

  use physicaldata, ONLY : flux_x, flux_y, flux_z, nfluxes
  use tree, ONLY : surr_blks, nodetype
  use Grid_data, ONLY : gr_xflx, gr_yflx, gr_zflx
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >=2
  use Grid_data, ONLY : gr_xflx_yface, gr_yflx_xface
#if NDIM == 3
  use Grid_data, ONLY : gr_xflx_zface, gr_yflx_zface, gr_zflx_xface, gr_zflx_yface
#endif
#endif
#endif

  implicit none

#include "constants.h"

  integer, intent(IN) :: blockID
  integer, intent(IN) :: axis
  integer, intent(IN), dimension(3) :: dataSize
  real, intent(IN), dimension(NFLUXES,dataSize(1),dataSize(2),dataSize(3)) :: fluxes
  integer, intent(IN), OPTIONAL,target :: pressureSlots(:)
  real, intent(IN), OPTIONAL :: areaLeft(:,:,:)

#if NFLUXES > 0
  integer :: presVar, np
  integer,save,dimension(1),target :: presDefault = (/-1/)
  integer,pointer,dimension(:) :: presP
  integer :: sx,ex,sy,ey,sz,ez

  if (present(pressureSlots)) then
     presP => pressureSlots
  else
     presP => presDefault
  end if


  sx = NGUARD+1
  sy = NGUARD*K2D+1
  sz = NGUARD*K3D+1
  ex = dataSize(1)-NGUARD
  ey = dataSize(2)-NGUARD*K2D
  ez = dataSize(3)-NGUARD*K3D

  select case(axis)
  case(IAXIS)
     flux_x(:nfluxes,1,:,:,blockID) = fluxes(:,sx,sy:ey,sz:ez) 
     flux_x(:nfluxes,2,:,:,blockID) = fluxes(:,ex+1,sy:ey,sz:ez) 
#ifdef CHOMBO_COMPATIBLE_HYDRO
     gr_xflx(:,1,:,:,blockID) = fluxes(:,sx,sy:ey,sz:ez)
     gr_xflx(:,2,:,:,blockID) = fluxes(:,ex+1,sy:ey,sz:ez)
#else
     gr_xflx(:,1,:,:,blockID) = fluxes(:,sx+1,sy:ey,sz:ez)
     gr_xflx(:,2,:,:,blockID) = fluxes(:,ex,sy:ey,sz:ez)
#ifdef FLASH_HYDRO_UNSPLIT
     !! Store transverse components for the faces in global scratch arrays.
#if NDIM >= 2
     gr_xflx_yface(:,:,1,:,blockID) = fluxes(:,sx+1:ex,sy,sz:ez)
     gr_xflx_yface(:,:,2,:,blockID) = fluxes(:,sx+1:ex,ey,sz:ez)
#if NDIM == 3
     gr_xflx_zface(:,:,:,1,blockID) = fluxes(:,sx+1:ex,sy:ey,sz)
     gr_xflx_zface(:,:,:,2,blockID) = fluxes(:,sx+1:ex,sy:ey,ez)
#endif
#endif
#endif
#endif
     do np = 1,size(presP,1)
        presVar = presP(np)
        if (presVar > 0) then
           if (.NOT.(surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,1,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
              where (areaLeft(sx,sy:ey,sz:ez).NE.0.0)
                 flux_x(presVar,1,:,:,blockID) = flux_x(presVar,1,:,:,blockID) * areaLeft(sx,sy:ey,sz:ez)
              end where
           end if
           if (.NOT.(surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,3,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
              flux_x(presVar,2,:,:,blockID) = flux_x(presVar,2,:,:,blockID) * areaLeft(ex+1,sy:ey,sz:ez)
           end if
        end if
     end do

  case(JAXIS)
     flux_y(:nfluxes,:,1,:,blockID)  = fluxes(:,sx:ex,sy,sz:ez)
     flux_y(:nfluxes,:,2,:,blockID)  = fluxes(:,sx:ex,ey+1,sz:ez) 
#ifdef CHOMBO_COMPATIBLE_HYDRO
     gr_yflx(:,:,1,:,blockID) =  fluxes(:,sx:ex,sy,sz:ez)
     gr_yflx(:,:,2,:,blockID) =  fluxes(:,sx:ex,ey+1,sz:ez)
#else
     gr_yflx(:,:,1,:,blockID) =  fluxes(:,sx:ex,sy+1,sz:ez)
     gr_yflx(:,:,2,:,blockID) =  fluxes(:,sx:ex,ey,sz:ez)
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >= 2
     !! Store transverse components for the faces in global scratch arrays.
     gr_yflx_xface(:,1,:,:,blockID) = fluxes(:,sx,sy+1:ey,sz:ez)
     gr_yflx_xface(:,2,:,:,blockID) = fluxes(:,ex,sy+1:ey,sz:ez)
#if NDIM == 3
     gr_yflx_zface(:,:,:,1,blockID) = fluxes(:,sx:ex,sy+1:ey,sz)
     gr_yflx_zface(:,:,:,2,blockID) = fluxes(:,sx:ex,sy+1:ey,ez)
#endif
#endif
#endif
#endif
#if NDIM > 1
     do np = 1,size(presP,1)
        presVar = presP(np)
        if (presVar > 0) then
           if (.NOT.(surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,2,1,1+K3D,blockID) == nodetype(blockID))) then
              where (areaLeft(sx:ex,sy,sz:ez).NE.0.0)
                 flux_y(presVar,:,1,:,blockID) = flux_y(presVar,:,1,:,blockID) * areaLeft(sx:ex,sy,sz:ez)
              end where
           end if
           if (.NOT.(surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,2,3,1+K3D,blockID) == nodetype(blockID))) then
              where (areaLeft(sx:ex,ey+1,sz:ez).NE.0.0)
                 flux_y(presVar,:,2,:,blockID) = flux_y(presVar,:,2,:,blockID) * areaLeft(sx:ex,ey+1,sz:ez)

              end where
           end if
        end if
     end do
#endif

  case(KAXIS)
     flux_z(:nfluxes,:,:,1,blockID) = fluxes(:,sx:ex,sy:ey,sz) 
     flux_z(:nfluxes,:,:,2,blockID) = fluxes(:,sx:ex,sy:ey,ez+1) 
#ifdef CHOMBO_COMPATIBLE_HYDRO
     gr_zflx(:,:,:,1,blockID) = fluxes(:,sx:ex,sy:ey,sz)
     gr_zflx(:,:,:,2,blockID) = fluxes(:,sx:ex,sy:ey,ez+1)
#else
     gr_zflx(:,:,:,1,blockID) = fluxes(:,sx:ex,sy:ey,sz+1)
     gr_zflx(:,:,:,2,blockID) = fluxes(:,sx:ex,sy:ey,ez)
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM == 3
     !! Store transverse components for the faces in global scratch arrays.
     gr_zflx_xface(:,1,:,:,blockID) = fluxes(:,sx,sy:ey,sz+1:ez)
     gr_zflx_xface(:,2,:,:,blockID) = fluxes(:,ex,sy:ey,sz+1:ez)

     gr_zflx_yface(:,:,1,:,blockID) = fluxes(:,sx:ex,sy,sz+1:ez)
     gr_zflx_yface(:,:,2,:,blockID) = fluxes(:,sx:ex,ey,sz+1:ez)
#endif
#endif
#endif
#if NDIM > 2
     do np = 1,size(presP,1)
        presVar = presP(np)
        if (presVar > 0) then
           if (.NOT.(surr_blks(1,2,2,1,blockID) > 0 .AND. &
             surr_blks(3,2,2,1,blockID) == nodetype(blockID))) then
              flux_z(presVar,:,:,1,blockID) = flux_z(presVar,:,:,1,blockID) * areaLeft(sx:ex,sy:ey,sz)
           end if
           if (.NOT.(surr_blks(1,2,2,3,blockID) > 0 .AND. &
             surr_blks(3,2,2,3,blockID) == nodetype(blockID))) then
              flux_z(presVar,:,:,2,blockID) = flux_z(presVar,:,:,2,blockID) * areaLeft(sx:ex,sy:ey,ez+1)
           end if
        end if
     end do
#endif
  end select
#endif

  return
end subroutine Grid_putFluxData
