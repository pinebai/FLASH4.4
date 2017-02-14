!!****if* source/Simulation/SimulationMain/unitTest/RungeKutta/2Dellipse/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash ()
!!
!! DESCRIPTION
!!
!!  Simple stripped down version for testing single units.
!!
!! NOTES
!!
!!***

subroutine Driver_evolveFlash ()

  implicit none

  real :: told, tnew

  call cpu_time (told)

  call sim_RungeKuttaTest ()

  call cpu_time (tnew)

  write (*,*) ' Cpu time = ',tnew - told

  return
end subroutine Driver_evolveFlash
