!!****f* source/physics/RadTrans/RadTrans_mgdEFromT
!!
!!  NAME 
!!
!!  RadTrans_mgdEFromT
!!
!!  SYNOPSIS
!!
!!  call RadTrans_mgdEFromT( integer(IN) :: blockId
!!                           integer(IN) :: axis(MDIM)
!!                           real(IN)    :: trad
!!                           real(OUT), optional :: tradActual )
!!
!!  DESCRIPTION 
!!
!!  This routine uses a radiation temperature to set the group
!!  specific radiation energies for all of the groups owned by this
!!  process in a given cell and block. It also sets the total specific
!!  radiation energy.
!!
!!  This routine was created to make it easy to initialize the group
!!  energies from a temperature in Simulation_initBlock.
!!
!!  The final, optional argument is tradActual. This is the radiation
!!  temperature computed from the specific energy. This is necessary
!!  to handle the case where the highest radiation energy group
!!  boundary is too small. In this case part of the black body
!!  spectrum will be cutoff and ERAD_VAR will not exactly equal
!!  a*trad**4. Thus, this routine will optionally return the correct
!!  value of trad which accounts for this possibility. Users can use
!!  this value to set TRAD_VAR in Simulation_initBlock.
!!
!!  ARGUMENTS
!!
!!    blockId    : The blockId of the cell
!!    axis       : An array storing the i,j,k coordinate of the cell
!!    trad       : The radiation temperature (K)
!!    tradActual : The actual radiation temperature computed from erad
!!***
subroutine RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
  implicit none

#include "constants.h"

  ! Arguments:
  integer, intent(in) :: blockId
  integer, intent(in) :: axis(MDIM)
  real,    intent(in) :: trad
  real,    intent(out), optional :: tradActual

  if (present(tradActual)) tradActual = trad

end subroutine RadTrans_mgdEFromT
