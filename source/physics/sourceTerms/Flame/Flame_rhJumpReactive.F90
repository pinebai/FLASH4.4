!!****f* source/physics/sourceTerms/Flame/Flame_rhJumpReactive
!!
!! NAME
!!
!!  Flame_rhJumpReactive
!!
!! SYNOPSIS
!!
!!  call Flame_rhJumpReactive ( real(inout) :: eosData_u(:),
!!                              real(in) :: qbar_u,
!!                              real(out) :: eosData_b(:),
!!                              real(out) :: qbar_b,
!!                              integer(in) :: eos_mode )
!!
!! DESCRIPTION
!!
!!  Calculate state of burned (and unburned) material by applying
!!  Rankine-Hugoniot jump condition.
!!  This version is for reactive ash (NSE), so that all information
!!  including the energy release (for a constant pressure burn) is
!!  calculated.
!!  Unburned state is found from calling Eos on eosData_u with mode
!!  eos_mode
!!  qbar values are in MeV/Baryon
!!
!! ARGUMENTS
!!
!!   eosData_u - thermodynamic state data for unburned material
!!      qbar_u - average binding energy per baryon (MeV) of unburned material
!!   eosData_b - thermodynamic state data for burned material
!!      qbar_u - average binding energy per baryon (MeV) of burned material
!!    eos_mode - equation of state mode for calculating complete unburned state.
!!
!! SEE ALSO
!!
!!  See Flame_interface.F90 for possible updates
!!
!!***

! This is a stub for when the flame unit is not compiled in
!
! Dean Townsley 2007,2008
!

#include "Eos.h"
subroutine Flame_rhJumpReactive(eosData_u, qbar_u, eosData_b, qbar_b, eos_mode)

   implicit none

   real, dimension(EOS_NUM), intent(inout) :: eosData_u
   real,    intent(in)                     :: qbar_u
   real, dimension(EOS_NUM), intent(out)   :: eosData_b
   real,    intent(out)                    :: qbar_b
   integer, intent(in)                     :: eos_mode

   eosData_b(:) = eosData_u
   qbar_b = qbar_u
end subroutine Flame_rhJumpReactive
