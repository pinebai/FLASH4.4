!!****if* source/Driver/DriverMain/Driver_initNumericalTools
!!
!! NAME
!!   
!!  Driver_initNumericalTools
!!
!! SYNOPSIS
!!
!!  Driver_initNumericalTools ()
!!
!! DESCRIPTION
!!
!!  Initializes all numerical tool Units by calling their respective
!!  initialization routines viz. Roots_init, RungeKutta_init, etc.
!!  
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine Driver_initNumericalTools ()

  use Roots_interface,      ONLY:  Roots_init
  use RungeKutta_interface, ONLY:  RungeKutta_init

  implicit none

  call Roots_init      ()
  call RungeKutta_init ()

end subroutine Driver_initNumericalTools
