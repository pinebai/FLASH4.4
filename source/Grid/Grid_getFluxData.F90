!!****f* source/Grid/Grid_getFluxData
!!
!! NAME
!!  Grid_getFluxData
!!
!! SYNOPSIS
!!
!!
!!  call Grid_getFluxData(integer(IN) :: blockID,
!!                   integer(IN) :: axis,
!!                   real(INOUT) :: fluxes(NFLUXES,dataSize(1),dataSize(2),dataSize(3)),
!!                   integer(IN) :: dataSize(3),
!!          OPTIONAL,integer(IN) :: pressureSlots,
!!          OPTIONAL,real(IN)    :: areaLeft(:,:,:))
!!
!! DESCRIPTION 
!!
!!  Get the fluxes in a direction specified by axis for boundary cells
!!  for block blockID. This routine needs to be used when using adaptive mesh
!!  since fluxes calculated by the two blocks that are at a fine/coarse
!!  boundary have different accuracy.
!!  
!!  This should be called after Grid_conserveFluxes, which gets the 
!!  fluxes updated for cells at fine/coarse boundaries and makes them consistent.
!!
!! ARGUMENTS
!!
!!  blockID : The local blockid
!!
!!
!!  axis : integer value specifying on which cell faces to get fluxes. 
!!         The options are IAXIS, JAXIS, or KAXIS defined in constants.h
!!
!!
!!  fluxes :  real array with space for fluxes, through one axis, 
!!            for all cells of a block and for all flux variables.
!!            fluxes(VAR, i, j, k) is VAR's flux through 
!!            the left cell face for cell i, j, k.  The 
!!            fluxes of the boundary cells of coarse blocks
!!            will have been appropriately changed, if the flux matching
!!            between fine and coarse boundaries (normally by a call to
!!            Grid_conserveFluxes), has occured since the last call
!!            of Grid_putFluxData.
!!
!!
!!  dataSize : integer array specifying the dimensions for fluxes
!!
!!             dataSize (1) holds the number of cells returned in the i direction
!!
!!             dataSize (2) holds the number of cells returned in the j direction
!!                          if 1 d problem, set datasize(2) = 1
!!
!!             dataSize (3) holds the number of cells returned in the k direction
!!                          if 1 or 2 d problem, set datasize(3) = 1
!!
!!             fluxes should contain space for fluxes of all cells in the block, 
!!             including guardcells.
!!
!!  pressureSlots: an optional array of integer values indicating fluxes that are to be
!!                 handled like pressures, rather than genuine fluxes.  If present, each
!!                 element greater than zero indicates one flux variable in the fluxes
!!                 array that may need special handling because it really does not scale
!!                 like a flux. Normally this would be pressure, or (e.g., with multiTemp
!!                 physics) several partial pressures given at cell interfaces. The
!!                 pressureSlots array could also list other flux variables that the
!!                 caller keeps in flux density form.
!!
!!  areaLeft :     areas of left and right faces, only used if special scaling is
!!                 requested with the pressureSlots argument.
!!
!! NOTES 
!!
!!   Specific implementations may ignore the optional dummy arguments,
!!   pressureSlots and areaLeft, if the corresponding Grid_conserveFluxes
!!   only implements handling of fluxes as flux densities (i.e., if
!!   Grid_conserveFluxes assumes that the fluxes stored by Grid_putFluxData
!!   and retrieved by Grid_getFluxData contain flux densities).
!!
!!   For the Unform Grid, no implementation is provided, since there is never
!!   a need for flux correction.
!!
!!   Any code calling this subroutine needs to know the explicit interface,
!!   since this interface contains optional dummy arguments and assumed-shape
!!   dummy arrays. Calling FORTRAN units should therefore contain a line like
!!       use Grid_interface, ONLY: Grid_getFluxData
!!
!! SEE ALSO
!!
!!   Grid_putFluxData
!!   Grid_conserveFluxes
!!   hy_sweep
!!***

subroutine Grid_getFluxData(blockID, axis, fluxes, dataSize, pressureSlots, areaLeft)

  implicit none

#include "Flash.h"

  integer, intent(IN) :: blockID
  integer, intent(IN) :: axis
  integer, intent(IN), dimension(3) :: dataSize
  real, intent(INOUT), dimension(NFLUXES,dataSize(1),dataSize(2),dataSize(3)) :: fluxes
  integer, intent(IN), OPTIONAL,target :: pressureSlots(:)
  real, intent(IN), OPTIONAL :: areaLeft(:,:,:)

  !dummy values for stub - not required for INOUT dummy argument
  !!fluxes = 0

  return
end subroutine Grid_getFluxData





