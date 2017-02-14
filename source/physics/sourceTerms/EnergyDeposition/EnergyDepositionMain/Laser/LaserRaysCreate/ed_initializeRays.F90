!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRaysCreate/ed_initializeRays
!!
!! NAME
!!
!!  ed_initializeRays
!!
!! SYNOPSIS
!!
!!  call ed_initializeRays ()
!!
!! DESCRIPTION
!!
!!  Initializes the complete ray array(s) and sets the counter(s) to zero.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!  Besides inititalizing always the main rays array, the routine also inititalizes a
!!  rays saving array if requested. Currently this is activated by only one keyword, but
!!  this might change in the future.
!!
!!***

subroutine ed_initializeRays ()

  use EnergyDeposition_data,   ONLY : ed_maxRayCount,         &
                                      ed_notSetInteger,       &
                                      ed_notSetReal,          &
                                      ed_numberOfSavedRays,   &
                                      ed_rayCount,            &
                                      ed_rays,                &
                                      ed_raysSaved,           &
                                      ed_saveOutOfDomainRays

  implicit none

#include "EnergyDeposition.h"
!
!
!     ...Initialize the complete main ray array and (optionally) the rays saving array.
!
!
  ed_rayCount          = 0
  ed_numberOfSavedRays = 0

  ed_rays (1:RAY_ATTR_COUNT , 1:ed_maxRayCount)  = ed_notSetReal

  if (ed_saveOutOfDomainRays) then
      ed_raysSaved (1:ed_maxRayCount) % rayX     = ed_notSetReal
      ed_raysSaved (1:ed_maxRayCount) % rayY     = ed_notSetReal
      ed_raysSaved (1:ed_maxRayCount) % rayZ     = ed_notSetReal
      ed_raysSaved (1:ed_maxRayCount) % rayPower = ed_notSetReal
      ed_raysSaved (1:ed_maxRayCount) % rayTag   = ed_notSetInteger
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_initializeRays
