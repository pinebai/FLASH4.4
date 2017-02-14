!!****ih* source/Simulation/SimulationMain/radflaHD/BondiAccretion/sim_interface
!!
!! This module includes interface blocks for some subroutines used
!! internally in the Bubble Lab simulation.
!! public interfaces.
!!***
Module sim_interface

  implicit none

  interface
     subroutine sim_computeAnaBondiScaled(x,u,alpha)
       implicit none
       real,intent(IN) :: x
       real,intent(OUT) :: u,alpha
     end subroutine sim_computeAnaBondiScaled

     subroutine sim_computeAnaBondi(r,velR,rho)
       use Simulation_data, ONLY: sim_rho_vac
       implicit none
       real,intent(IN) :: r
       real,intent(OUT) :: velR,rho
     end subroutine sim_computeAnaBondi
  end interface

end Module sim_interface
