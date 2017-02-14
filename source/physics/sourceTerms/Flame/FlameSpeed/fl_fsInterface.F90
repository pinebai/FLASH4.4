! Dean Townsley 2008

#include "FortranLangFeatures.fh"

module fl_fsInterface

  interface fl_fsInit
    subroutine fl_fsInit()
       implicit none
    end subroutine fl_fsInit
  end interface

  interface fl_fsFinalize
     subroutine fl_fsFinalize()
        implicit none
     end subroutine
  end interface

  interface fl_flameSpeed
     subroutine fl_flameSpeed( solnData, flamespeed, blockID, nlayers)
        implicit none
        real, dimension(:,:,:,:),POINTER_INTENT_IN :: solnData
        real, dimension(:,:,:),intent(out) :: flamespeed
        integer, intent(in) :: blockID, nlayers
     end subroutine
  end interface

  interface fl_fsGcMask
     subroutine fl_fsGcMask(fl_gcMask,fl_gcDoEos)
        implicit none
        logical, dimension(:), intent(inout) ::  fl_gcMask
        logical, intent(inout) :: fl_gcDoEos
        ! This function should add to the mask the variables
        ! needed to calculate the flame speed.  Also setting
        ! the doeos flag as necessary
        ! This intended to be able to be used in an "accumulation"
        ! thus the flags should only be changed if .true. is required,
        ! otherwise the values on input should be preserved
     end subroutine
  end interface

end module
