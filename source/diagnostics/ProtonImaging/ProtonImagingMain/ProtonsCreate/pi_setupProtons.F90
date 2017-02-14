!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonsCreate/pi_setupProtons
!!
!! NAME
!!
!!  pi_setupProtons
!!
!! SYNOPSIS
!!
!!  call pi_setupProtons ()
!!
!! DESCRIPTION
!!
!!  Allocates the needed proton array and initializes the counter for the proton
!!  tags.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pi_setupProtons ()

  use ProtonImaging_data,    ONLY : pi_maxProtonCount, &
                                    pi_tagMax,         &
                                    pi_protons

  implicit none

#include "ProtonImaging.h"
!
!
!     ...Allocate the proton array.
!
!
  allocate (pi_protons (1:PROTON_ATTRCOUNT,1:pi_maxProtonCount))

  pi_tagMax = 0
!
!
!     ...Ready!
!
!
  return
end subroutine pi_setupProtons
