!!****if* source/Simulation/SimulationMain/unitTest/Roots/x4Polynomials/sim_setx4PolynomialCoeffs
!!
!! NAME
!!
!!  sim_setx4PolynomialCoeffs
!! 
!! SYNOPSIS
!!
!!  call sim_setx4PolynomialCoeffs
!!
!! DESCRIPTION
!!
!!  Sets the x4 polynomial coefficients A,B,C,D for the x^3, x^2, x^1 and x^0 terms:
!!
!!              x4 polynomial = x^4 + Ax^3 + Bx^2 + Cx + D
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_setx4PolynomialCoeffs ()

  use  Driver_interface, ONLY: Driver_abortFlash

  use  Simulation_data,  ONLY: sim_numberOfx4Polynomials, &
                               sim_x4Polynomialx0Coeff,   &
                               sim_x4Polynomialx1Coeff,   &
                               sim_x4Polynomialx2Coeff,   &
                               sim_x4Polynomialx3Coeff

  implicit none
!
!
!   ...Check array storage size.
!
!
  if (sim_numberOfx4Polynomials > 13) then
      call Driver_abortFlash ('[sim_setx4PolynomialCoeffs] ERROR: Not enough storage for coefficients!')
  end if
!
!
!   ...Set the coefficients.
!
!
  sim_x4Polynomialx3Coeff (1:13) = (/ -1.001001001e+9,       &
                                      -4.006,                &
                                      -9.98899e+79,          &
                                      -1.00000000000002e+14, &
                                      +1.e+7,                &
                                      -9.000002e+6,          &
                                      +2.000011e+6,          &
                                      -1.00002011e+8,        &
                                      -2.0000002e+7,         &
                                      -1.9986e+4,            &
                                      -4.006,                &
                                      -4.0e+3,               &
                                      -6.0                   /)

  sim_x4Polynomialx2Coeff (1:13) = (/ +1.001002001001e+15,   &
                                      +6.018011,             &
                                      -1.1008989e+157,       &
                                      +1.99999999999999e+14, &
                                      -2.00000000000001e+14, &
                                      -0.9999981999998e+13,  &
                                      +1.010022000028e+12,   &
                                      +2.01101022001e+11,    &
                                      +1.01000040000005e+14, &
                                      +1.00720058e+8,        &
                                      +5.6008018e+1,         &
                                      +6.00001e+6,           &
                                      +1.01000013e+8         /)

  sim_x4Polynomialx1Coeff (1:13) = (/ -1.001001001e+18,      &
                                      -4.018022006,          &
                                      -1.010999e+233,        &
                                      +1.00000000000002e+14, &
                                      -1.e+7,                &
                                      +1.9999982e+13,        &
                                      +1.1110056e+13,        &
                                      -1.02200111000011e+14, &
                                      -2.020001e+14,         &
                                      -1.8600979874e+10,     &
                                      -1.04148036024e+2,     &
                                      -4.00002e+9,           &
                                      -2.04000012e+8         /)

  sim_x4Polynomialx0Coeff (1:13) = (/ +1.e+18,               &
                                      +1.006011006,          &
                                      -1.e+307,              &
                                      -2.e+14,               &
                                      +2.e+14,               &
                                      -2.e+13,               &
                                      +2.828e+13,            &
                                      +1.1000011e+15,        &
                                      +5.05e+14,             &
                                      +1.00004909000441e+14, &
                                      +6.75896068064016e+2,  &
                                      +1.000010000009e+12,   &
                                      +1.00000104000004e+14  /)
!
!
!   ...Ready!
!
!
  return
end subroutine sim_setx4PolynomialCoeffs
