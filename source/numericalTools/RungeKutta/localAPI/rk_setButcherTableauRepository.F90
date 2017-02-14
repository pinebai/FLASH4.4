!!****if* source/numericalTools/RungeKutta/localAPI/rk_setButcherTableauRepository
!!
!! NAME
!!
!!  rk_setButcherTableauRepository
!!
!! SYNOPSIS
!!
!!  call rk_setButcherTableauRepository ()
!!
!! DESCRIPTION
!!
!!  This routine sets up the whole collection of Butcher tableaus the Runge Kutta
!!  integrator is able to handle. It is the repository from which the Runge Kutta
!!  stepper will choose its appropriate tableaus.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine rk_setButcherTableauRepository ()

  implicit none

  return
end subroutine rk_setButcherTableauRepository
