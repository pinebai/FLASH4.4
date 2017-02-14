!!****if* source/Grid/GridMain/paramesh/Grid_conserveFluxes
!!
!! NAME
!!  Grid_conserveFluxes
!!
!! SYNOPSIS
!!
!!  Grid_conserveFluxes(integer(IN) :: axis,
!!                      integer(IN) :: level)
!!  
!! DESCRIPTION 
!!  
!!  Flux conservation is necessary when 2 blocks of differing
!!  levels (meaning having different grid spacings) border 
!!  one another. 
!!  
!!  This routine can perform flux conservation on the finest
!!  blocks, the most typical usage for the Paramesh Grid or on
!!  blocks of a certain level.
!!  
!!  The routine overwrites the flux arrays maintained by the Grid
!!  
!! ARGUMENTS 
!!
!!
!!  axis - conserve fluxes in just one direction if 
!!         IAXIS, JAXIS, KAXIS, or in all directions if ALLDIR.
!!         These constants are defined in constants.h.
!!
!!  level - refinement level. Ignored.
!!
!!***
subroutine Grid_conserveFluxes( axis, level)
  use paramesh_interfaces, ONLY : amr_flux_conserve
#include "Flash.h"
#include "constants.h"
#ifndef FLASH_GRID_PARAMESH2
  use physicaldata, ONLY: no_permanent_guardcells
#endif
  use Grid_data, ONLY : gr_meshMe
implicit none
  integer, intent(in) ::  axis, level
  integer :: gridDataStruct

  if (axis == ALLDIR) then
     call amr_flux_conserve(gr_meshMe, 0, 0)          
  else
     call amr_flux_conserve(gr_meshMe, 0, axis)     
  endif


#ifndef FLASH_GRID_PARAMESH2
  gridDataStruct = CENTER
#if NFACE_VARS > 0
  gridDataStruct = CENTER_FACES
#endif
  if (no_permanent_guardcells) then
     call gr_commSetup(gridDataStruct)
  else
     call gr_freeCommRecvBuffer
  end if
#endif
end subroutine Grid_conserveFluxes
