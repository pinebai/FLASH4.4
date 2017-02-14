!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRaysCreate/ed_setupRays
!!
!! NAME
!!
!!  ed_setupRays
!!
!! SYNOPSIS
!!
!!  call ed_setupRays ()
!!
!! DESCRIPTION
!!
!!  Sets up the rays array(s). This simply means to allocate the needed rays array(s).
!!  Since the rays will be created afresh at each time step, there is no point in
!!  initializing the rays array(s) at this stage. Initialization of the rays is done
!!  right before they are going to be created.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!  Besides setting up always the main rays array, the routine also sets up a (opional)
!!  rays saving array. Currently this is activated by only one keyword, but this might
!!  change in the future.
!!
!!***

subroutine ed_setupRays ()

  use EnergyDeposition_data,    ONLY : ed_maxRayCount,         &
                                       ed_rays,                &
                                       ed_raysSaved,           &
                                       ed_saveOutOfDomainRays

  implicit none

#include "EnergyDeposition.h"
!
!
!     ...Allocate the main rays array and (optionally) the rays saving array.
!
!
  allocate (ed_rays (1:RAY_ATTR_COUNT,1:ed_maxRayCount))

  if (ed_saveOutOfDomainRays) then
      allocate (ed_raysSaved (1:ed_maxRayCount))
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_setupRays
