!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beamsInfo
!!
!! NAME
!!
!!  ed_beamsInfo
!!
!! SYNOPSIS
!!
!!  call ed_beamsInfo ()
!!
!! DESCRIPTION
!!
!!  Generates information about the beams for specific grid geometries. In here all beam
!!  information is generated that can be generated at initialization. This routine calls
!!  the appropriate subroutines according to the domain geometry specified. If the geometry
!!  does not figure inside this routine it simply does nothing and exits.
!!
!! ARGUMENTS
!!
!!***

subroutine ed_beamsInfo ()

  use ed_interface,           ONLY : ed_beamsInfo1DRec, &
                                     ed_beamsInfo2DRec, &
                                     ed_beamsInfo3DRec

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

      call ed_beamsInfo1DRec ()

    case (GRID_2DCARTESIAN)

      call ed_beamsInfo2DRec ()

    case (GRID_2DCYLINDRICAL)

      if (ed_laser3Din2D) then
          call ed_beamsInfo3DRec ()
      else
          call ed_beamsInfo2DRec ()
      end if

    case (GRID_3DCARTESIAN)

      call ed_beamsInfo3DRec ()

    case default

      return

  end select
!
!
!    ...Ready!
!
!
  return
end subroutine ed_beamsInfo
