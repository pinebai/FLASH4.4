!!****f* source/physics/sourceTerms/Flame/Flame_rhJump
!!
!! NAME
!!
!!  Flame_rhJump
!!
!! SYNOPSIS
!!
!!  call Flame_rhJump (
!!                     real(inout) :: eosData_u(EOS_NUM),
!!                     real(inout) :: eosData_b(EOS_NUM),
!!                     real(in)    :: q,
!!                     real(in)    :: s,
!!                     integer(in) :: mode,
!!           optional, real(in)    :: mfrac_u(NSPECIES),
!!           optional, real(in)    :: mfrac_b(NSPECIES) )
!!
!! DESCRIPTION
!!
!!  Calculate the thermodynamic state of the burned material
!!  by applying Rankine-Hugoniot jump condition.
!!  Unburned state is calculated frome *_u with mode as the eos mode.
!!  Burned state is calculated from energy release q (in erg/gram) and
!!  the flame speed, s, (wrt the fuel) with the composition information
!!  of the fuel specified using either eosData_?(EOS_ABAR/ZBAR) or mfrac_?
!!  depending on whether EOS_YEYI is being used
!!
!! ARGUMENTS
!!
!!   eosData_u - thermodynamic state info about unburned material
!!   eosData_b - thermodynamic state info about burned material
!!      q      - energy release (erg/g)
!!      s      - flame speed (cm/s)
!!      mode   - EOS mode used to calculate full unburned state
!!    mfrac_u  - composition of unburned state (if not present abar and zbar from eosData_u is used)
!!    mfrac_b  - composition of burned state (if not present abar and zbar from eosData_b is used)
!!
!! SEE ALSO
!!
!!  See Flame_interface.F90 for possible updates
!!
!!***
!
! Dean Townsley 2007,2008
!

#include "Flash.h"
#include "Eos.h"
subroutine Flame_rhJump(eosData_u, eosData_b, q, s, mode, mfrac_u, mfrac_b)

  implicit none

  real, dimension(EOS_NUM),  intent(inout) :: eosData_u
  real, dimension(EOS_NUM),  intent(inout) :: eosData_b
  real,    intent(in)     :: q, s
  integer, intent(in)     :: mode  !! This is the Eos mode
  real, optional, dimension(NSPECIES), intent(in)    :: mfrac_u
  real, optional, dimension(NSPECIES), intent(in)    :: mfrac_b

end subroutine Flame_rhJump
