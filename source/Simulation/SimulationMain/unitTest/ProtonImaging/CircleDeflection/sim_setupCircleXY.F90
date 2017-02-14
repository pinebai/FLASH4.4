!!****if* source/Simulation/SimulationMain/unitTest/ProtonImaging/CircleDeflection/sim_setupCircleXY
!!
!! NAME 
!!
!!  sim_setupCircleXY
!!
!! SYNOPSIS
!!
!!  sim_setupCircleXY ()
!!
!! DESCRIPTION
!!
!!  This routine sets up beam proton (x,y) coordinates, such that a circle will result
!!  on the beam's crossectional plane. The radius of the circle created is equal to 1
!!  circular beam radius, hence the circle will lay on the surface of the beam.
!!
!! ARGUMENTS
!!
!!***

subroutine sim_setupCircleXY ()

  use Simulation_Data,  ONLY : sim_nBeamPoints, &
                               sim_beamPointsX, &
                               sim_beamPointsY

  implicit none

#include "constants.h"

  integer :: i
  real    :: angle, angularPart
!
!
!     ...Divide the beam crossectional area into # of cricular points angular parts.
!        The first point has an angle of 0 radians.
!
!
  angle = 0.0
  angularPart = (PI + PI) / real (sim_nBeamPoints)

  do i = 1,sim_nBeamPoints
     sim_beamPointsX (i) = sin (angle)
     sim_beamPointsY (i) = cos (angle)
     angle = angle + angularPart
  end do
!
!
!     ...Ready!
!
!
  return
end subroutine sim_setupCircleXY
