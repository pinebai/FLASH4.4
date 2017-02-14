!!****f* source/Grid/Grid_putFluxData
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
!!   Specific implementations may ignore the optional dummy arguments,
!!   pressureSlot and areaLeft, if the corresponding Grid_conserveFluxes
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
!!       use Grid_interface, ONLY: Grid_putFluxData
!!
!! SEE ALSO
!!
!!   Grid_getFluxData
!!   Grid_conserveFluxes
!!   hy_sweep
!!***

subroutine Grid_putFluxData(blockID, axis, fluxes, dataSize, pressureSlots, areaLeft)

  implicit none

#include "Flash.h"

  integer, intent(IN) :: blockID
  integer, intent(IN) :: axis
  integer, intent(IN), dimension(3) :: dataSize
  real, intent(IN), dimension(NFLUXES,dataSize(1),dataSize(2),dataSize(3)) :: fluxes
  integer, intent(IN), OPTIONAL,target :: pressureSlots(:)
  real, intent(IN), OPTIONAL :: areaLeft(:,:,:)

  return
end subroutine Grid_putFluxData





