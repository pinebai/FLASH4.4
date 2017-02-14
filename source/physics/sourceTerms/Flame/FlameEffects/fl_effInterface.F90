!
! Dean Townsley 2008
!

#include "FortranLangFeatures.fh"

module fl_effInterface
  interface fl_effInit
     subroutine fl_effInit()
        implicit none
     end subroutine fl_effInit
  end interface

  interface fl_effFinalize
     subroutine fl_effFinalize()
        implicit none
     end subroutine
  end interface

  interface fl_effects
     subroutine fl_effects( solnData, flamdot, dt, blockID)
        implicit none
        real, dimension(:,:,:,:), POINTER_INTENT_IN :: solnData
        real, dimension(:,:,:), intent(in)  :: flamdot
        real, intent(in) :: dt
        integer, intent(in) :: blockID
                ! Applies ancillary effects of the flame (e.g. energy release,
                ! composition change) by updating variables in the block passed
                ! via solnData
                ! Depending on the implementation, energy may be deposited
                ! in a unit called later (like Burn)
     end subroutine
  end interface

end module fl_effInterface
