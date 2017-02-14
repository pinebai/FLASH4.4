!!****h* source/physics/sourceTerms/Flame/Flame_interface
!!
!! NAME
!!
!!  Flame_interface
!!
!! SYNOPSIS
!!
!!  use Flame_interface
!!
!! DESCRIPTION
!!
!!  Interface module for the Flame code unit.
!!
!!***
! declarations of public interface for Flame tracking Unit
!
! Dean Townsley 2007,2008
!
! The basic function of each routine is also described here.


Module Flame_interface
#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  implicit none

  interface Flame_step
     subroutine Flame_step( num_blocks, blockList, dt  )    
       
       integer, INTENT(in)                        :: num_blocks
       integer, INTENT(in), DIMENSION(num_blocks) :: blockList
       real,    INTENT(in)                        :: dt
       !  Pricipal public function.
       !  Evolve flame forward for all blocks in blockList by one step
       !  of size dt.
       !  Flame speed and flame effects sub-units are called within this
       !  subroutine as they are required
       !  for ADR Applies unsplit reaction-diffusion operater to FLAM_MSCALAR
       !  May or may not deposit energy, depending on which
       !  Flame_Effects module had been included
     end subroutine Flame_step
  end interface

  interface
     subroutine Flame_computeAbarZbar(solnScalars, abarData, zbarData)
       implicit none
       real, dimension(:,:), intent(in)  :: solnScalars
       real, dimension(:), intent(inout) :: abarData, zbarData
       ! A callback, typically called by Eos unit implementations to get
       ! values for EOS_ABAR and EOS_ZBAR input elements of eosData
       ! before the EOS computation proper.
     end subroutine Flame_computeAbarZbar
  end interface

  interface Flame_getProfile
     subroutine Flame_getProfile(x, f)
       
       real, intent(in)  :: x
       real, intent(out) :: f
       !  get value of progress variable a distance x from center of
       !  flame front in the steady state propagating flame (used to
       !  initialize data on mesh).
       !  x is defined such that positive x is in the direction of
       !  propagation and f = 0.5 at x = 0
     end subroutine Flame_getProfile
  end interface

  interface Flame_getWidth
     subroutine Flame_getWidth(laminarWidth)

       real, intent(OUT) :: laminarWidth
       !  approximate total width of the flame front
       !  more than about twice this far away progress
       !  variable can be initialized  to 0 or 1 for
       !  unburned and burned respectively
     end subroutine Flame_getWidth
  end interface

  interface Flame_init
     subroutine Flame_init()

     end subroutine Flame_init
  end interface

  interface Flame_finalize
     subroutine Flame_finalize()

     end subroutine Flame_finalize
  end interface

  interface Flame_laminarSpeed
     subroutine Flame_laminarSpeed(dens, s, ds, info)

       real, intent(in)   :: dens
       real, intent(out)  :: s
       real, optional, intent(out) :: ds
       real, dimension(:),optional,intent(in)   :: info
       ! return the physical laminar flame speed s
       ! generally only used in initialization
       ! behavior is strongly dependent on what FlameSpeed subunit is included
       ! the info array is to provide imlementation-dependent information such
       ! as composition if applicable
     end subroutine Flame_laminarSpeed
  end interface

  interface
     subroutine Flame_rhJump(eosData_u, eosData_b, q, s, mode, mfrac_u, mfrac_b)

        real, dimension(EOS_NUM),  intent(inout) :: eosData_u
        real, dimension(EOS_NUM),  intent(inout) :: eosData_b
        real,    intent(in)     :: q, s
        integer, intent(in)     :: mode  !! This is the Eos mode
        real, optional, dimension(NSPECIES), intent(in)    :: mfrac_u
        real, optional, dimension(NSPECIES), intent(in)    :: mfrac_b
       ! Calculate the thermodynamic state of the burned material
       ! (and unburned) material by applying Rankine-Hugoniot jump condition.
       ! Unburned state is calculated frome *_u with mode as the eos mode.
       ! Burned state is calculated from energy release q (in erg/gram) and
       ! the flame speed, s, (wrt the fuel) with the composition information
       ! of the fuel specified using either eosData_?(EOS_ABAR/ZBAR) or mfrac_?
       ! depending on whether EOS_YEYI is being used
     end subroutine Flame_rhJump
  end interface

  interface Flame_rhJumpReactive
     subroutine Flame_rhJumpReactive(eosData_u, qbar_u, eosData_b, qbar_b, eos_mode)

        real, dimension(EOS_NUM), intent(inout) :: eosData_u
        real,    intent(in)                     :: qbar_u
        real, dimension(EOS_NUM), intent(out)   :: eosData_b
        real,    intent(out)                    :: qbar_b
        integer, intent(in)                     :: eos_mode
        ! Calculate state of burned (and unburned) material by applying
        ! Rankine-Hugoniot jump condition.
        ! This version is for reactive ash (NSE), so that all information
        ! including the energy release (for a constant pressure burn) is
        ! calculated.
        ! Unburned state is found from calling Eos on eosData_u with mode
        ! eos_mode
        ! qbar values are in MeV/Baryon
     end subroutine Flame_rhJumpReactive
  end interface

end Module Flame_interface
