!!****if* source/numericalTools/RungeKutta/localAPI/rk_orderRKmethod
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

  implicit none

  character (len=*), intent (in) :: method

  rk_orderRKmethod = 0

  return
end function rk_orderRKmethod
