!!****if* source/Grid/GridMain/paramesh/Paramesh2/Grid_putFluxData
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
!!  pressureSlot : If present and greater than zero, this indicates one flux variable
!!                 in the fluxes array that may need special handling because it
!!                 really scales like a flux; normally this would be pressure,
!!                 but it could be another flux variable that the caller keeps in
!!                 flux density form. Ignored in this implementation since with
!!                 Paramesh2, all "fluxes" are assumed to be give as flux densities
!!                 anyway as far as the Grid unit is concerned.
!!
!!  areaLeft :     areas of left and right faces, only used if special scaling is
!!                 requested with the pressureSlot argument.
!!                 Ignored in this implementation since with
!!                 Paramesh2, all "fluxes" are assumed to be give as flux densities
!!                 anyway as far as the Grid unit is concerned.
!!
!! NOTES 
!!
!!   This implementation is specific to Paramesh 2.
!!
!!   Any code calling this subroutine needs to know the explicit interface,
!!   since this interface contains optional dummy arguments and assumed-shape
!!   dummy arrays. Calling FORTRAN units should therefore contain a line like
!!       use Grid_interface, ONLY: Grid_putFluxData
!!
!!***


subroutine Grid_putFluxData(blockID, axis, fluxes, dataSize, pressureSlots, areaLeft)

  use physicaldata, ONLY : flux_x, tflux_x, flux_y, tflux_y, flux_z, tflux_z, nfluxes

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: blockID
  integer, intent(IN) :: axis
  integer, intent(IN), dimension(3) :: dataSize
  real, intent(IN), dimension(NFLUXES,dataSize(1),dataSize(2),dataSize(3)) :: fluxes
  integer, intent(IN), OPTIONAL,target :: pressureSlots(:)
  real, intent(IN), OPTIONAL :: areaLeft(:,:,:)

  integer :: leftIndex = NGUARD+1
  integer :: rightIndex 

  select case(axis)

  case(IAXIS)
     rightIndex=leftIndex+NXB-1
     flux_x(:nfluxes,1,:,:,blockID) = fluxes(:,leftIndex,:,:) 
     flux_x(:nfluxes,2,:,:,blockID) = fluxes(:,rightIndex+1,:,:) 
     tflux_x(:nfluxes,1,:,:,blockID) = fluxes(:,leftIndex+1,:,:) 
     tflux_x(:nfluxes,2,:,:,blockID) = fluxes(:,rightIndex,:,:) 

  case(JAXIS)
     rightIndex=leftIndex+NYB-1
     flux_y(:nfluxes,:,1,:,blockID) = fluxes(:,:,leftIndex,:) 
     flux_y(:nfluxes,:,2,:,blockID) = fluxes(:,:,rightIndex+1,:) 
     tflux_y(:nfluxes,:,1,:,blockID) = fluxes(:,:,leftIndex+1,:) 
     tflux_y(:nfluxes,:,2,:,blockID) = fluxes(:,:,rightIndex,:) 

  case(KAXIS)
     rightIndex=leftIndex+NZB-1
     flux_z(:nfluxes,:,:,1,blockID) = fluxes(:,:,:,leftIndex) 
     flux_z(:nfluxes,:,:,2,blockID) = fluxes(:,:,:,rightIndex+1) 
     tflux_z(:nfluxes,:,:,1,blockID) = fluxes(:,:,:,leftIndex+1) 
     tflux_z(:nfluxes,:,:,2,blockID) = fluxes(:,:,:,rightIndex) 

  end select

  return
end subroutine Grid_putFluxData





