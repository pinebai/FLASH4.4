!!****ih* source/Simulation/SimulationMain/unitTest/RungeKutta/3Dcircle/sim_interface
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
    function sim_ODEfunction (t,y)
      real, intent (in) :: t
      real, intent (in) :: y (:)
      real              :: sim_ODEfunction (1:size (y))
    end function sim_ODEfunction
  end interface

end Module sim_interface
