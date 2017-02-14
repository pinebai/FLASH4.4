!!****if* source/Simulation/SimulationMain/unitTest/Roots/x3Polynomials/sim_setx3PolynomialCoeffs
!!
!! NAME
!!
!!  sim_setx3PolynomialCoeffs
!! 
!! SYNOPSIS
!!
!!  call sim_setx3PolynomialCoeffs
!!
!! DESCRIPTION
!!
!!  Sets the x3 polynomial coefficients A,B,C for the x^2, x^1 and x^0 terms:
!!
!!              x3 polynomial = x^3 + Ax^2 + Bx + C
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_setx3PolynomialCoeffs ()

  use  Driver_interface, ONLY: Driver_abortFlash

  use  Simulation_data,  ONLY: sim_numberOfx3Polynomials, &
                               sim_x3Polynomialx0Coeff,   &
                               sim_x3Polynomialx1Coeff,   &
                               sim_x3Polynomialx2Coeff

  implicit none
!
!
!   ...Check array storage size.
!
!
  if (sim_numberOfx3Polynomials > 9) then
      call Driver_abortFlash ('[sim_setx3PolynomialCoeffs] ERROR: Not enough storage for coefficients!')
  end if
!
!
!   ...Set the coefficients.
!
!
  sim_x3Polynomialx2Coeff (1:9) = (/ -1.00000010000001e+14, &
                                     -3.000003,             &
                                     +8.999e+80,            &
                                     +1.e+24,               &
                                     -1.99999999999999e+14, &
                                     -3.e+5,                &
                                     -3.0,                  &
                                     -2.0000001e+7,         &
                                     +0.99999999999998e+14  /)

  sim_x3Polynomialx1Coeff (1:9) = (/ +1.00000010000001e+21, &
                                     +3.000006000002,       &
                                     -1.0009e+161,          &
                                     -1.0,                  &
                                     +0.99999999999998e+28, &
                                     +3.0000000001e+10,     &
                                     +1.00000000000003e+14, &
                                     +1.00000020000001e+14, &
                                     -1.99999999999998e+14  /)

  sim_x3Polynomialx0Coeff (1:9) = (/ -1.e+21,               &
                                     -1.000003000002,       &
                                     +1.e+238,              &
                                     -1.e+24,               &
                                     +1.e+28,               &
                                     -1.0000000001e+15,     &
                                     -1.00000000000001e+14, &
                                     -1.00000000000001e+14, &
                                     +2.e+14                /)
!
!
!   ...Ready!
!
!
  return
end subroutine sim_setx3PolynomialCoeffs
