!!****if* source/numericalTools/RungeKutta/RungeKuttaMain/rk_orderRKmethod
!!
!! NAME
!!
!!  rk_orderRKmethod
!!
!! SYNOPSIS
!!
!!  rk_orderRKmethod (character (len=*), intent (in) :: method)
!!
!! DESCRIPTION
!!
!!  This integer function returns the order of the specified Runge Kutta method.
!!
!! ARGUMENTS
!!
!!  method : the Runge Kutta method
!!
!! NOTES
!!
!!***

integer function rk_orderRKmethod (method)

  use Driver_interface,  ONLY : Driver_abortFlash

  implicit none

  character (len=*), intent (in) :: method
!
!
!   ...Retrieve the order.
!
!
  select case (method)

  case ('EulerHeu12'); rk_orderRKmethod = 1
  case ('BogShamp23'); rk_orderRKmethod = 2
  case ('Fehlberg34'); rk_orderRKmethod = 3
  case ('Fehlberg45'); rk_orderRKmethod = 4
  case ('CashKarp45'); rk_orderRKmethod = 4

  case default
        call Driver_abortFlash ('[rk_orderRKmethod] ERROR: unknown RK method')
  end select
!
!
!   ...Ready! 
!
!
  return
end function rk_orderRKmethod
