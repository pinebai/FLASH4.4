!!****if* source/numericalTools/RungeKutta/RungeKuttaMain/RungeKutta_finalize
!!
!! NAME
!!  
!!  RungeKutta_finalize
!!
!! SYNOPSIS
!! 
!!  call RungeKutta_finalize ()
!!
!! DESCRIPTION
!!
!!  Finalizes the RungeKutta unit.
!!
!! ARGUMENTS
!!
!!***

subroutine RungeKutta_finalize ()

  use RungeKutta_data, ONLY: rk_a, rk_b, rk_c

  implicit none
!
!
!   ...Deassociate all pointers.
!
!
  if (associated (rk_a)) nullify (rk_a)
  if (associated (rk_b)) nullify (rk_b)
  if (associated (rk_c)) nullify (rk_c)
!
!
!    ...Ready!
!
!
  return
end subroutine RungeKutta_finalize
