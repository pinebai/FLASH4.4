!!****if* source/Simulation/SimulationMain/ConductionDeltaSaDiff/sim_xToXi
!!
!! NAME
!!
!!  sim_xToXi
!!
!! SYNOPSIS
!!
!!  call sim_xToXi(real(in) :: x,
!!                 real(in) :: t,
!!                 real(in) :: n,
!!                 real(out) :: xi)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   x : 
!!
!!   t : 
!!
!!   n : 
!!
!!   xi : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine sim_xToXi(x, t, n, xi)
  use Simulation_data, ONLY: sim_xi0, sim_alpha, sim_Q, sim_toffset, &
       sim_condTemperatureExponent
  implicit none
  real, intent(in) :: x, t, n
  real, intent(out) :: xi

!!$  n = sim_condTemperatureExponent

  xi = x / ( (sim_alpha * sim_Q**n * (t+sim_toffset))**(1.0/(n+2)) )
end subroutine sim_xToXi
