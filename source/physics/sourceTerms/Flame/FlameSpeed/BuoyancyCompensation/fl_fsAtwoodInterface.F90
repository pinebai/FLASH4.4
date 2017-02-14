!
! Dean Townsley 2008

#include "constants.h"

module fl_fsAtwoodInterface

  ! no real reason to declare explicit interfaces for
  ! Init and Finalize subroutines because they take no arguments

  interface fl_flAtwood
     subroutine fl_fsAtwood(solnData, atwood, dens_u, calcLimits)

        implicit none

        real, dimension(:,:,:,:), pointer :: solnData
        real, dimension(:,:,:), intent(out) :: atwood
        real, dimension(:,:,:), intent(in)  :: dens_u
        integer, dimension(LOW:HIGH,MDIM), intent(in)    :: calcLimits

        ! fills the atwood array with an estimate of the atwood number
        ! at the fuel density dens_u.  The solnData structure is used for
        ! ancillary data (such as initial composition).  Only the
        ! cells in the ranges indicated by calcLimits are filled.

     end subroutine
  end interface

end module
