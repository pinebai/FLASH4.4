!!****f* source/physics/Eos/Eos_guardCells
!!
!! NAME
!!
!!  Eos_guardCells
!!
!! SYNOPSIS
!!
!!  call Eos_guardCells(integer(IN)  :: eosmode,
!!                      integer(IN)  :: blockid,
!!                      logical(IN)  :: corners,
!!             optional,integer(IN)  :: layers(MDIM),
!!             optional,logical(IN)  :: skipSrl)
!!
!! DESCRIPTION
!!
!!  Another layer of wrapping around Eos_wrapped calls, provided
!!  as a convenience to make it easy to apply the EOS to guard cells only.
!!  
!! ARGUMENTS
!!
!!  eosmode : determines which variables are used as Eos input variables.
!!            The valid values are MODE_DENS_EI (where density and internal
!!            energy are inputs), MODE_DENS_PRES (density and pressure as inputs)
!!            MODE_DENS_TEMP (density and temperature are inputs).
!!            These quantities are defined in constants.h.
!!            The argument is passed unchanged and unexamined to Eos_wrapped calls.
!!
!!  blockid : ID of block in current processor
!!
!!  corners : indicates whether Eos should be called on corner
!!            guard cells (i.e., diagonal guard cells)
!!  layers  : the number of guard cells to be included along each dimension
!!  skipSrl : whether to skip guard cell regions that are coming from
!!            neighboring blocks at the same refinement (or from boundary
!!            conditions) and thus have not undergone interpolation or
!!            restrictions.
!!
!! SEE ALSO
!!  Eos_wrapped
!!
!!***

subroutine Eos_guardCells(eosMode, blockID,corners,layers,skipSrl)

#include "constants.h"
  implicit none

  integer,intent(IN) :: eosMode,blockID
  logical,intent(IN) :: corners
  integer,dimension(MDIM),optional, intent(IN) :: layers
  logical,optional, intent(IN) :: skipSrl

end subroutine Eos_guardCells
