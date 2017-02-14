!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beamsCheck
!!
!! NAME
!!
!!  ed_beamsCheck
!!
!! SYNOPSIS
!!
!!  call ed_beamsCheck ()
!!
!! DESCRIPTION
!!
!!  Checks the collected beams data for specific grid geometries. This routine calls the appropriate
!!  subroutines according to the domain geometry specified. If the geometry does not figure inside
!!  this routine it simply does nothing and exits.
!!
!! ARGUMENTS
!!
!!***

subroutine ed_beamsCheck ()

  use ed_interface,           ONLY : ed_beamsCheck1DRec,   &
                                     ed_beamsCheck2DCyl3D, &
                                     ed_beamsCheck2DRec,   &
                                     ed_beamsCheck3DRec

  use EnergyDeposition_data,  ONLY : ed_gridGeometry, &
                                     ed_laser3Din2D

  implicit none

#include "EnergyDeposition.h"
!
!
!     ...Select the appropriate subroutine.
!
!
  select case (ed_gridGeometry)

    case (GRID_1DCARTESIAN)

      call ed_beamsCheck1DRec ()

    case (GRID_2DCARTESIAN)

      call ed_beamsCheck2DRec ()

    case (GRID_2DCYLINDRICAL)

      if (ed_laser3Din2D) then
          call ed_beamsCheck2DCyl3D ()
      else
          call ed_beamsCheck2DRec ()
      end if

    case (GRID_3DCARTESIAN)

      call ed_beamsCheck3DRec ()

    case default

      return

  end select
!
!
!    ...Ready!
!
!
  return
end subroutine ed_beamsCheck
