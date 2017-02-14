!!****ih* source/Simulation/SimulationMain/unitTest/RungeKutta/2Dellipse/sim_interface
!!
!! NAME
!!
!!  sim_interface
!!
!! SYNOPSIS
!!
!!  use sim_interface
!!
!! DESCRIPTION
!!
!!  The interface for the current simulation unit.
!!
!!***

Module sim_interface

  interface
    subroutine sim_calculateInitialData ()
    end subroutine sim_calculateInitialData
  end interface

  interface
     function sim_distancePoint2Ellipse2Dxy (a, b, Rot, Cx, Cy, Px, Py)
       real, intent (in)  :: a, b
       real, intent (in)  :: Rot
       real, intent (in)  :: Cx, Cy
       real, intent (in)  :: Px, Py
       real               :: sim_distancePoint2Ellipse2Dxy (1:2)
     end function sim_distancePoint2Ellipse2Dxy
  end interface

  interface
    function sim_ODEfunction (t,y)
      real, intent (in) :: t
      real, intent (in) :: y (:)
      real              :: sim_ODEfunction (1:size (y))
    end function sim_ODEfunction
  end interface

end Module sim_interface
