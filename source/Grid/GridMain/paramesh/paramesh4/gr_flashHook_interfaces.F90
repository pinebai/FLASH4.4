!!****ih* source/Grid/GridMain/paramesh/paramesh4/gr_flashHook_interfaces
!!
!! NAME
!!
!!  gr_flashHook_interfaces
!!
!! SYNOPSIS
!!
!!  use gr_flashHook_interfaces
!!
!! DESCRIPTION
!!
!!  Interfaces for some hooks that are called back from within PARAMESH 4,
!!  so that FLASH can provide the conversion to (and from) conserved form
!!  in mesh interpolation.
!!
!! NOTES
!!
!!  Names and organization of these subprograms could change,
!!  if calling them becomes part of official PARAMESH code.
!!  Currently the calls to these routines are FLASH modifications
!!  of the PARAMESH code.
!!
!!***

! Modification history:
!     Created   October 2006  KW

      module gr_flashHook_interfaces

      implicit none
      
      integer,parameter :: &
      &  gr_cells_INTERIOR=1, gr_cells_GUARD=2, gr_cells_ALL=3

      integer,parameter :: &
      gr_callReason_PROLONG=1, gr_callReason_RESTRICT=2


      interface
         subroutine flash_convert_cc_hook(datainout, nvars, i1,i2,j1,j2,k1,k2, why)
         integer, intent(IN) :: nvars, i1,i2,j1,j2,k1,k2
         real,dimension(nvars, i1:i2,j1:j2,k1:k2),intent(INOUT) :: datainout
         integer, intent(IN) :: why
         end subroutine flash_convert_cc_hook
      end interface

      interface
         subroutine flash_unconvert_cc_hook(datainout, nvars, i1,i2,j1,j2,k1,k2, where, why, &
     &                                      nlayers_in_data)
         integer, intent(IN) :: nvars, i1,i2,j1,j2,k1,k2
         real,dimension(nvars, i1:i2,j1:j2,k1:k2),intent(INOUT) :: datainout
         integer, intent(IN) :: where, why
         integer, intent(IN), optional :: nlayers_in_data
         end subroutine flash_unconvert_cc_hook
      end interface

      end module gr_flashHook_interfaces
