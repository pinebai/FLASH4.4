!
! Aaron Jackson 2010

#include "constants.h"
#include "Flash.h"

module fl_fsTFIInterface

  ! no real reason to declare explicit interfaces for
  ! Init and Finalize subroutines because they take no arguments

  interface fl_fsTFIFlameSpeedBlock
     subroutine fl_fsTFIFlameSpeedBlock(solnData, s, ds, dx, compLimits, &
                                        quench_limit)
       implicit none

       real, dimension(:,:,:,:), pointer :: solnData
       real, dimension(:,:,:), intent(inout) :: s
       real, dimension(:,:,:), intent(in) :: ds
       real, intent(in) :: dx
       integer, dimension(LOW:HIGH,MDIM), intent(in) :: compLimits
       real, dimension(:,:,:), intent(out), optional :: quench_limit

       ! the laminar flame speed is supplied and the turbulent
       ! flame speed is returned.
    end subroutine fl_fsTFIFlameSpeedBlock
  end interface

  interface fl_fsTFIEnhance
     subroutine fl_fsTFIEnhance(E, up, s, de, dl0, de_over_dl1, E_lim)
       implicit none

       real, intent(out) :: E
       real, intent(in) :: up, s, de, dl0, de_over_dl1
       real, intent(out), optional :: E_lim

       ! calculates the enhancement factor to the laminar flame speed
     end subroutine fl_fsTFIEnhance
  end interface

  interface fl_fsTFIGamma
     subroutine fl_fsTFIGamma(G, up_over_s, de_over_dl)
       implicit none

       real, intent(out) :: G
       real, intent(in) :: up_over_s, de_over_dl

       ! calculates Gamma
     end subroutine fl_fsTFIGamma
  end interface

  interface fl_fsTFIFunc
     real function fl_fsTFIFunc(x)
        implicit none
        real, intent(in) :: x
     end function fl_fsTFIFunc
  end interface

  interface fl_fsTFIAlpha
     subroutine fl_fsTFIAlpha(alpha, up, de, dl0)
       implicit none

       real, intent(out) :: alpha
       real, intent(in)  :: up, de, dl0

       ! calculates model constant alpha from Colin et al. (2000) TFI
       ! implementation. fl_fsBeta is a runtime parameter and is used to
       ! calculate alpha. Any scaling of alpha should use fl_fsBeta.
     end subroutine
  end interface

  interface fl_fsTFIGcMask
     subroutine fl_fsTFIGcMask(fl_gcMask, fl_gcDoEos)
       implicit none

       logical, dimension(:), intent(inout) :: fl_gcMask
       logical, intent(inout) :: fl_gcDoEos

       ! sets gc mask for any variable the TFI model needs
     end subroutine
  end interface

end module
