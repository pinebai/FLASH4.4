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

  interface fl_fsUnburnDensBlock
     subroutine fl_fsUnburnDensBlock(solnData, dens, calcLimits)
     implicit none

     real,dimension(:,:,:,:),pointer :: solnData
     real,dimension(:,:,:),intent(OUT) :: dens
     integer,dimension(LOW:HIGH,MDIM),intent(IN) :: calcLimits
     ! Calculates the unburned density over the block provided.  The structure
     ! calcLimits contains the index ranges over which the calculation is to be
     ! performed. It is assumed that the output arrays have been correctly
     ! sized by the calling routine to be able to hold data within the
     ! specified index ranges.
     end subroutine
  end interface

! private functions for Laminar flamespeed sub-sub Unit

  interface
     subroutine fl_fsConeInterp(c12i,ne22i,den,s,ds)
       implicit none

       real, intent(IN) :: c12i, ne22i, den
       real, intent(OUT) :: s, ds
       ! performs a single-point interpolation of the laminar nuclear flame
       ! speed at the specified mass fractions of carbon 12 and Neon 22 and
       ! density
     end subroutine
  end interface

  interface
     subroutine fl_fsUnburnDens(pres, ye, dens)
       implicit none

       real, intent(IN)  :: pres, ye
       real, intent(OUT) :: dens
       ! estimates the unburned density from the pressure and composition
       ! comes from Hansen and Kawaler and fit for White Dwarf
     end subroutine
  end interface

end module

