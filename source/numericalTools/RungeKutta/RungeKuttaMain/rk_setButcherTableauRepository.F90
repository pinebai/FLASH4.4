!!****if* source/numericalTools/RungeKutta/RungeKuttaMain/rk_setButcherTableauRepository
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

  use RungeKutta_data,   ONLY : rk_aTableauCashKarp45, &
                                rk_bTableauCashKarp45, &
                                rk_cTableauCashKarp45, &
                                rk_aTableauBogShamp23, &
                                rk_bTableauBogShamp23, &
                                rk_cTableauBogShamp23, &
                                rk_aTableauFehlberg34, &
                                rk_bTableauFehlberg34, &
                                rk_cTableauFehlberg34, &
                                rk_aTableauFehlberg45, &
                                rk_bTableauFehlberg45, &
                                rk_cTableauFehlberg45, &
                                rk_aTableauEulerHeu12, &
                                rk_bTableauEulerHeu12, &
                                rk_cTableauEulerHeu12

  implicit none
!
!
!     ...Add the Euler-Heun 1(2) tableaus to the repository.
!
!
  rk_cTableauEulerHeu12 (1) = 0.0
  rk_cTableauEulerHeu12 (2) = 1.0

  rk_bTableauEulerHeu12 (1,1) = 1.0 / 2.0
  rk_bTableauEulerHeu12 (2,1) = 1.0 / 2.0

  rk_bTableauEulerHeu12 (1,2) = 1.0
  rk_bTableauEulerHeu12 (2,2) = 0.0

  rk_aTableauEulerHeu12 (2,1) = 1.0
!
!
!     ...Add the Bogacki-Shampine 2(3) tableaus to the repository.
!
!
  rk_cTableauBogShamp23 (1) = 0.0
  rk_cTableauBogShamp23 (2) = 1.0 / 2.0
  rk_cTableauBogShamp23 (3) = 3.0 / 4.0
  rk_cTableauBogShamp23 (4) = 1.0

  rk_bTableauBogShamp23 (1,1) = 2.0 / 9.0
  rk_bTableauBogShamp23 (2,1) = 1.0 / 3.0
  rk_bTableauBogShamp23 (3,1) = 4.0 / 9.0
  rk_bTableauBogShamp23 (4,1) = 0.0

  rk_bTableauBogShamp23 (1,2) = 7.0 / 24.0
  rk_bTableauBogShamp23 (2,2) = 1.0 / 4.0
  rk_bTableauBogShamp23 (3,2) = 1.0 / 3.0
  rk_bTableauBogShamp23 (4,2) = 1.0 / 8.0

  rk_aTableauBogShamp23 (2,1) = 1.0 / 2.0
  rk_aTableauBogShamp23 (3,1) = 0.0
  rk_aTableauBogShamp23 (4,1) = 2.0 / 9.0
  rk_aTableauBogShamp23 (3,2) = 3.0 / 4.0
  rk_aTableauBogShamp23 (4,2) = 1.0 / 3.0
  rk_aTableauBogShamp23 (4,3) = 4.0 / 9.0
!
!
!     ...Add the Fehlberg 3(4) tableaus to the repository.
!
!
  rk_cTableauFehlberg34 (1) = 0.0
  rk_cTableauFehlberg34 (2) = 1.0 / 4.0
  rk_cTableauFehlberg34 (3) = 4.0 / 9.0
  rk_cTableauFehlberg34 (4) = 6.0 / 7.0
  rk_cTableauFehlberg34 (5) = 1.0

  rk_bTableauFehlberg34 (1,1) = 1.0 / 6.0
  rk_bTableauFehlberg34 (2,1) = 0.0
  rk_bTableauFehlberg34 (3,1) = 27.0 / 52.0
  rk_bTableauFehlberg34 (4,1) = 49.0 / 156.0
  rk_bTableauFehlberg34 (5,1) = 0.0

  rk_bTableauFehlberg34 (1,2) = 43.0 / 288.0
  rk_bTableauFehlberg34 (2,2) = 0.0
  rk_bTableauFehlberg34 (3,2) = 243.0 / 416.0
  rk_bTableauFehlberg34 (4,2) = 343.0 / 1872.0
  rk_bTableauFehlberg34 (5,2) = 1.0 / 12.0

  rk_aTableauFehlberg34 (2,1) = 1.0 / 4.0
  rk_aTableauFehlberg34 (3,1) = 4.0 / 81.0
  rk_aTableauFehlberg34 (4,1) = 57.0 / 98.0
  rk_aTableauFehlberg34 (5,1) = 1.0 / 6.0
  rk_aTableauFehlberg34 (3,2) = 32.0 / 81.0
  rk_aTableauFehlberg34 (4,2) = - 432.0 / 343.0
  rk_aTableauFehlberg34 (5,2) = 0.0
  rk_aTableauFehlberg34 (4,3) = 1053.0 / 686.0
  rk_aTableauFehlberg34 (5,3) = 27.0 / 52.0
  rk_aTableauFehlberg34 (5,4) = 49.0 / 156.0
!
!
!     ...Add the Fehlberg 4(5) tableaus to the repository.
!
!
  rk_cTableauFehlberg45 (1) = 0.0
  rk_cTableauFehlberg45 (2) = 1.0 / 4.0
  rk_cTableauFehlberg45 (3) = 3.0 / 8.0
  rk_cTableauFehlberg45 (4) = 12.0 / 13.0
  rk_cTableauFehlberg45 (5) = 1.0
  rk_cTableauFehlberg45 (6) = 1.0 / 2.0

  rk_bTableauFehlberg45 (1,1) = 16.0 / 135.0
  rk_bTableauFehlberg45 (2,1) = 0.0
  rk_bTableauFehlberg45 (3,1) = 6656.0 / 12825.0
  rk_bTableauFehlberg45 (4,1) = 28561.0 / 56430.0
  rk_bTableauFehlberg45 (5,1) = - 9.0 / 50.0
  rk_bTableauFehlberg45 (6,1) = 2.0 / 55.0

  rk_bTableauFehlberg45 (1,2) = 25.0 / 216.0
  rk_bTableauFehlberg45 (2,2) = 0.0
  rk_bTableauFehlberg45 (3,2) = 1408.0 / 2565.0
  rk_bTableauFehlberg45 (4,2) = 2197.0 /4104.0
  rk_bTableauFehlberg45 (5,2) = - 1.0 / 5.0
  rk_bTableauFehlberg45 (6,2) = 0.0

  rk_aTableauFehlberg45 (2,1) = 1.0 / 4.0
  rk_aTableauFehlberg45 (3,1) = 3.0 / 32.0
  rk_aTableauFehlberg45 (4,1) = 1932.0 / 2197.0
  rk_aTableauFehlberg45 (5,1) = 439.0 / 216.0
  rk_aTableauFehlberg45 (6,1) =  - 8.0 / 27.0
  rk_aTableauFehlberg45 (3,2) = 9.0 / 32.0
  rk_aTableauFehlberg45 (4,2) = - 7200.0 / 2197.0
  rk_aTableauFehlberg45 (5,2) = - 8.0
  rk_aTableauFehlberg45 (6,2) = 2.0
  rk_aTableauFehlberg45 (4,3) = 7296.0 / 2197.0
  rk_aTableauFehlberg45 (5,3) = 3680.0 / 513.0
  rk_aTableauFehlberg45 (6,3) = - 3544.0 / 2565.0
  rk_aTableauFehlberg45 (5,4) = - 845.0 / 4104.0
  rk_aTableauFehlberg45 (6,4) = 1859.0 / 4104.0
  rk_aTableauFehlberg45 (6,5) = - 11.0 / 40.0
!
!
!     ...Add the Cash Karp 4(5) tableaus to the repository.
!
!
  rk_cTableauCashKarp45 (1) = 0.0
  rk_cTableauCashKarp45 (2) = 1.0 / 5.0
  rk_cTableauCashKarp45 (3) = 3.0 / 10.0
  rk_cTableauCashKarp45 (4) = 3.0 / 5.0
  rk_cTableauCashKarp45 (5) = 1.0
  rk_cTableauCashKarp45 (6) = 7.0 / 8.0

  rk_bTableauCashKarp45 (1,1) = 37.0 / 378.0
  rk_bTableauCashKarp45 (2,1) = 0.0
  rk_bTableauCashKarp45 (3,1) = 250.0 / 621.0
  rk_bTableauCashKarp45 (4,1) = 125.0 / 594.0
  rk_bTableauCashKarp45 (5,1) = 0.0
  rk_bTableauCashKarp45 (6,1) = 512.0 / 1771.0

  rk_bTableauCashKarp45 (1,2) = 2825.0 / 27648.0
  rk_bTableauCashKarp45 (2,2) = 0.0
  rk_bTableauCashKarp45 (3,2) = 18575.0 / 48384.0
  rk_bTableauCashKarp45 (4,2) = 13525.0 / 55296.0
  rk_bTableauCashKarp45 (5,2) = 277.0 / 14336.0
  rk_bTableauCashKarp45 (6,2) = 1.0 / 4.0

  rk_aTableauCashKarp45 (2,1) = 1.0 / 5.0
  rk_aTableauCashKarp45 (3,1) = 3.0 / 40.0
  rk_aTableauCashKarp45 (4,1) = 3.0 / 10.0
  rk_aTableauCashKarp45 (5,1) = - 11.0 / 54.0
  rk_aTableauCashKarp45 (6,1) = 1631.0 / 55296.0
  rk_aTableauCashKarp45 (3,2) = 9.0 / 40.0
  rk_aTableauCashKarp45 (4,2) = - 9.0 / 10.0
  rk_aTableauCashKarp45 (5,2) = 5.0 / 2.0
  rk_aTableauCashKarp45 (6,2) = 175.0 / 512.0
  rk_aTableauCashKarp45 (4,3) = 6.0 / 5.0
  rk_aTableauCashKarp45 (5,3) = - 70.0 / 27.0
  rk_aTableauCashKarp45 (6,3) = 575.0 / 13824.0
  rk_aTableauCashKarp45 (5,4) = 35.0 / 27.0
  rk_aTableauCashKarp45 (6,4) = 44275.0 / 110592.0
  rk_aTableauCashKarp45 (6,5) = 253.0 / 4096.0
!
!
!   ...Ready! 
!
!
  return
end subroutine rk_setButcherTableauRepository
