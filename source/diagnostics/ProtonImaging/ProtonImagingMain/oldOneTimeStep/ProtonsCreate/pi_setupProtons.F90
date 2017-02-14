!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonsCreate/pi_setupProtons
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
!!  Sets up the protons, which means to allocate the needed proton array. The proton array
!!  is used to track the protons through the domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pi_setupProtons ()

  use ProtonImaging_data,    ONLY : pi_maxProtonCount, &
                                    pi_protons

  implicit none

#include "ProtonImaging.h"
!
!
!     ...Allocate the proton array.
!
!
  allocate (pi_protons (1:PROTON_ATTRCOUNT,1:pi_maxProtonCount))
!
!
!     ...Ready!
!
!
  return
end subroutine pi_setupProtons
