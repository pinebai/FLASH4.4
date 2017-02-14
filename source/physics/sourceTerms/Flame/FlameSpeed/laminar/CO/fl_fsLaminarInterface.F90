!
! Dean Towsley 2008

#include "constants.h"

module fl_fsLaminarInterface

  interface fl_fsLaminarFlameSpeedBlock
     subroutine fl_fsLaminarFlameSpeedBlock(solnData, s, dens, calcLimits, ds)
     implicit none

     real,dimension(:,:,:,:),pointer :: solnData
     real,dimension(:,:,:),intent(OUT) :: s
     real,dimension(:,:,:),intent(OUT) :: dens
     integer,dimension(LOW:HIGH,MDIM), intent(IN) :: calcLimits
     real,OPTIONAL,dimension(:,:,:),intent(OUT) :: ds
     ! Calculates the unburned density and the speed of the laminar nuclear
     ! flame, s,  over the block provided.  The structure calcLimits contains
     ! the index ranges over which the calculation is to be performed.  It is
     ! assumed that the output arrays have been correctly sized by the calling
     ! routine to be able to hold data within the specified index ranges.
     end subroutine
  end interface

end module

