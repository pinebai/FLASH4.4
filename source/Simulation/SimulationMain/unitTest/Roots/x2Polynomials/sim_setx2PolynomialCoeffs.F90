!!****if* source/Simulation/SimulationMain/unitTest/Roots/x2Polynomials/sim_setx2PolynomialCoeffs
!!
!! NAME
!!
!!  sim_setx2PolynomialCoeffs
!! 
!! SYNOPSIS
!!
!!  call sim_setx2PolynomialCoeffs
!!
!! DESCRIPTION
!!
!!  Sets the x2 polynomial coefficients A,B for the x^1 and x^0 terms:
!!
!!                   x2 polynomial = x^2 + Ax + B
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_setx2PolynomialCoeffs ()

  use  Driver_interface, ONLY: Driver_abortFlash

  use  Simulation_data,  ONLY: sim_LPN, sim_SPN,          &
                               sim_numberOfx2Polynomials, &
                               sim_x2Polynomialx0Coeff,   &
                               sim_x2Polynomialx1Coeff

  implicit none
!
!
!   ...Check array storage size.
!
!
  if (sim_numberOfx2Polynomials > 16) then
      call Driver_abortFlash ('[sim_setx2PolynomialCoeffs] ERROR: Not enough storage for coefficients!')
  end if
!
!
!   ...Set the coefficients for: x^2 + ax + b.
!
!
  sim_x2Polynomialx1Coeff (1:16) = (/ + sim_LPN, &
                                      + sim_LPN, &
                                      - sim_LPN, &
                                      - sim_LPN, &
                                      + sim_LPN, &
                                      + sim_LPN, &
                                      - sim_LPN, &
                                      - sim_LPN, &
                                      + sim_SPN, &
                                      + sim_SPN, &
                                      - sim_SPN, &
                                      - sim_SPN, &
                                      + sim_SPN, &
                                      + sim_SPN, &
                                      - sim_SPN, &
                                      - sim_SPN  /)

  sim_x2Polynomialx0Coeff (1:16) = (/ + sim_LPN, &
                                      - sim_LPN, &
                                      + sim_LPN, &
                                      - sim_LPN, &
                                      + sim_SPN, &
                                      - sim_SPN, &
                                      + sim_SPN, &
                                      - sim_SPN, &
                                      + sim_LPN, &
                                      - sim_LPN, &
                                      + sim_LPN, &
                                      - sim_LPN, &
                                      + sim_SPN, &
                                      - sim_SPN, &
                                      + sim_SPN, &
                                      - sim_SPN  /)
!
!
!   ...Ready!
!
!
  return
end subroutine sim_setx2PolynomialCoeffs
