!!****if* source/numericalTools/RungeKutta/RungeKuttaMain/RungeKutta_data
!!
!! NAME
!!
!!  RungeKutta_data
!!
!! SYNOPSIS
!!
!!  use RungeKutta_data
!!  
!! DESCRIPTION
!!
!!  Data module for RungeKutta ODE integration
!!  ------------------------------------------
!!   
!!   Legend: (P) means data that is set as (P)arameters
!!           (G) means data that is (G)et from other units (driver, physical constants, etc ...)
!!           (R) means data that is supplied as input through (R)untime parameters
!!           (C) means data that is (C)alculated internally by the Runge Kutta code
!!           (I) means data that is (I)nitialized internally by the Runge Kutta code
!!
!!   rk_a                            (C) : Pointers for a(i,j)
!!   rk_aTableauBogShamp23           (I) : Stores a(i,j) for the Bogacki-Shampine 2(3) method
!!   rk_aTableauCashKarp45           (I) : Stores a(i,j) for the Cash Karp 4(5) method
!!   rk_aTableauEulerHeu12           (I) : Stores a(i,j) for the Euler-Heun 1(2) method
!!   rk_aTableauFehlberg34           (I) : Stores a(i,j) for the Fehlberg 3(4) method
!!   rk_aTableauFehlberg45           (I) : Stores a(i,j) for the Fehlberg 4(5) method
!!   rk_b                            (C) : Pointers for b(i), b*(i)
!!   rk_bTableauBogShamp23           (I) : Stores b(i), b*(i) for the Bogacki-Shampine 2(3) method
!!   rk_bTableauCashKarp45           (I) : Stores b(i), b*(i) for the Cash Karp 4(5) method
!!   rk_bTableauEulerHeu12           (I) : Stores b(i), b*(i) for the Euler-Heun 1(2) method
!!   rk_bTableauFehlberg34           (I) : Stores b(i), b*(i) for the Fehlberg 3(4) method
!!   rk_bTableauFehlberg45           (I) : Stores b(i), b*(i) for the Fehlberg 4(5) method
!!   rk_c                            (C) : Pointers for c(i)
!!   rk_cTableauBogShamp23           (I) : Stores c(i) for the Bogacki-Shampine 2(3) method
!!   rk_cTableauCashKarp45           (I) : Stores c(i) for the Cash Karp 4(5) method
!!   rk_cTableauEulerHeu12           (I) : Stores c(i) for the Euler-Heun 1(2) method
!!   rk_cTableauFehlberg34           (I) : Stores c(i) for the Fehlberg 3(4) method
!!   rk_cTableauFehlberg45           (I) : Stores c(i) for the Fehlberg 4(5) method
!!   rk_cubeRootMacheps            (C,I) : Cube root of the machine epsilon value
!!   rk_kVectors                     (C) : Stores the intermediate k vectors
!!   rk_maxButcherTableauDimension   (P) : Maximum size of the (n x n) Butcher tableau
!!   rk_maxNumberDependentVariables  (P) : Maximum number of dependent variables
!!   rk_stepSizeConfinementFactor    (R) : Reduction factor for step size reduction for confined RK runs
!!   rk_stepSizeSafetyFactor         (R) : Build in safety factor for new step size estimate
!!
!! NOTES
!!
!!  1) The maximum number of variables and Butcher tableaus that the Runge-Kutta can
!!     handle is 'hardwired' as parameters. This is done in order to create threadsave
!!     working arrays and pointers for possible multithreaded applications. The attribute
!!     'threadprivate' is only allowed on arrays that have NOT been allocated right
!!     before a multithreaded region opens up. Hence these arrays must not be allocated
!!     during initialization of the Runge-Kutta unit. An allocation/deallocation pair
!!     for each Runge-Kutta step is prohibitive due to performance reasons.
!!
!!***

Module RungeKutta_data

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, parameter :: rk_maxButcherTableauDimension  = 6
  integer, parameter :: rk_maxNumberDependentVariables = 20

  real, save :: rk_cubeRootMacheps
  real, save :: rk_stepSizeConfinementFactor
  real, save :: rk_stepSizeSafetyFactor

  real, target, save :: rk_aTableauBogShamp23 (1:4,1:4)
  real, target, save :: rk_aTableauCashKarp45 (1:6,1:6)
  real, target, save :: rk_aTableauEulerHeu12 (1:2,1:2)
  real, target, save :: rk_aTableauFehlberg34 (1:5,1:5)
  real, target, save :: rk_aTableauFehlberg45 (1:6,1:6)

  real, target, save :: rk_bTableauBogShamp23 (1:4,1:2)
  real, target, save :: rk_bTableauCashKarp45 (1:6,1:2)
  real, target, save :: rk_bTableauEulerHeu12 (1:2,1:2)
  real, target, save :: rk_bTableauFehlberg34 (1:5,1:2)
  real, target, save :: rk_bTableauFehlberg45 (1:6,1:2)

  real, target, save :: rk_cTableauBogShamp23 (1:4)
  real, target, save :: rk_cTableauCashKarp45 (1:6)
  real, target, save :: rk_cTableauEulerHeu12 (1:2)
  real, target, save :: rk_cTableauFehlberg34 (1:5)
  real, target, save :: rk_cTableauFehlberg45 (1:6)

  real, pointer :: rk_c (:)
  real, pointer :: rk_b (:,:)
  real, pointer :: rk_a (:,:)

  real, save :: rk_initialFunction (1:rk_maxNumberDependentVariables)
  real, save :: rk_kVectors        (1:rk_maxNumberDependentVariables , 1:rk_maxButcherTableauDimension)

  !$omp threadprivate (rk_a, rk_b, rk_c, rk_initialFunction, rk_kVectors)

end Module RungeKutta_data
