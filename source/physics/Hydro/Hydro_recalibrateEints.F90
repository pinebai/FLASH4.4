!!****f* source/physics/Hydro/Hydro_recalibrateEints
!! NAME
!!
!!  Hydro_recalibrateEints
!! 
!! SYNOPSIS
!!
!!  call Hydro_recalibrateEints(
!!                     integer(IN) :: range(HIGH, MDIM),
!!                     integer(IN) :: blockID)
!!
!!  call Hydro_recalibrateEintsForCell(
!!                     integer(in)    :: eint,
!!                     integer(INOUT) :: eion,
!!                     integer(INOUT) :: eele,
!!            OPTIONAL,integer(INOUT) :: erad,
!!            OPTIONAL,integer(INOUT) :: e1,e2,e3)
!!  call Hydro_recalibrateEints(
!!                     integer(in)    :: eint,
!!                     integer(INOUT) :: eion,
!!                     integer(INOUT) :: eele,
!!            OPTIONAL,integer(INOUT) :: erad,
!!            OPTIONAL,integer(INOUT) :: e1,e2,e3)
!!
!!
!! DESCRIPTION
!!
!! This function recalibrates multiTemp component internal energies and energies
!! so that their sums agree with the overall values for eint.
!!
!! In applications not configured for multiTemp Hydro, the call will act like
!! stubs and leave all arguments unchanged.
!!
!! The assumption is that
!!
!!    eint = eion + eele + erad    (in general)
!!    eint = eion + eele           (in the one-cell form if erad is not present;
!!                                  in the per-range/block form if runtime parameter
!!                                  hy_3Ttry_D==3.0).
!!
!! In the per-range/block form of the routine, those variables represent
!! the permanent solution variables (in UNK) EINT_VAR, EION_VAR, EELE_VAR,
!! and ERAD_VAR.
!!
!! If the equation is not true, or is only approximately try, the call will make
!! it true if possible by multiplying the provided two or three component energies
!! with a common scaling factor.
!!
!! Hydro_recalibrateEints is available as a generic function.
!! There is a specific implementation that acts on all cells within a range for
!! a given block.
!! Another specific implementation, available with the specific name
!! Hydro_recalibrateEintsForCell, can be used to adjust energies fo the
!! components of one cell.
!!
!!  ARGUMENTS 
!!
!!   
!!   range: an array that holds the lower and upper indices of the section
!!          of block on which recalibration is to be applied. The example
!!          shows how the array describes the block section.
!!
!!   blockID: current block number
!!
!!   eint     : combined internal energy (input)
!!   eion,eele: "ion" and "electron" component energies; to be adjusted.
!!   erad     : radiation component energy; to be adjusted, if present.
!!   e1,e2,e3 : alternative energy variables to be co-adjusted with eion,eele,
!!              erad, respectively,  if present (may not be useful)
!!
!!  EXAMPLE 
!!      if range(LOW,IAXIS)=1,range(HIGH,IAXIS)=iguard,
!!         range(LOW,JAXIS)=1,range(HIGH,JAXIS)=jguard,
!!         range(LOW,KAXIS)=1,range(HIGH,KAXIS)=kguard,
!!      then recalibration is applied to the lower left hand corner of the guard
!!      cells in the block. 
!!
!!      However if the value were
!!         range(LOW,IAXIS)=iguard+1,range(HIGH,IAXIS)=iguard+nxb,
!!         range(LOW,JAXIS)=jguard+1,range(HIGH,JAXIS)=jguard+nyb,
!!         range(LOW,KAXIS)=kguard+1,range(HIGH,KAXIS)=kguard+nzb,
!!      then recalibration is applied to all the interior cells in the block.
!!
!!  NOTES
!!      This interface is defined in Fortran Module 
!!      Hydro_interface. All functions calling this routine should include
!!      a statement like
!!      use Hydro_interface, ONLY : Hydro_recalibrateEints
!!
!!      This routine cannot use "INTERIOR" mode of indexing the range.  In the
!!      second example given above, although only the interior cells are being
!!      calculated with EOS, the range indices still must include the guard cells.
!!      See, for example, IsentropicVortex/Simulation_initBlock where the data is
!!      generated on INTERIOR cells with Grid_putRowData, but the same indices cannot
!!      be used for the EOS call.
!!
!!  SEE ALSO
!!
!!     Eos_wrapped
!!     Eos
!!     Eos.h
!!
!!***

subroutine Hydro_recalibrateEints(range,blockID)

  implicit none

#include "constants.h"

  integer, dimension(2,MDIM), intent(in) :: range
  integer,intent(in) :: blockID

  return
end subroutine Hydro_recalibrateEints

subroutine Hydro_recalibrateEintsForCell(eint,eion,eele,erad,e1,e2,e3)
  implicit none
  real,intent(in)    :: eint
  real,intent(INOUT) :: eion,eele
  real,intent(INOUT),OPTIONAL :: erad
  real,intent(INOUT),OPTIONAL :: e1,e2,e3

end subroutine Hydro_recalibrateEintsForCell
