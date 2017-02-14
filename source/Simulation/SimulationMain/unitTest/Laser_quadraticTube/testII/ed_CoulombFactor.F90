!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/ed_CoulombFactor
!!
!! NAME
!!
!!  ed_CoulombFactor
!!
!! SYNOPSIS
!!
!!  ed_CoulombFactor (real (in) :: Z,
!!                    real (in) :: e,
!!                    real (in) :: k,
!!                    real (in) :: T,
!!                    real (in) :: Ne)
!!
!! DESCRIPTION
!!
!!  Special overriding version for computation of the Coulomb factor 'ln (Lambda)'. It is
!!  simply set equal to 1. in order to compare to the analytical solution, which was derived
!!  based on the assumption that the Coulomb term Lambda is independent of the electron
!!  number density and the electron temperature.
!!
!!  This routine is needed to override the version in the Energy Deposition unit. Otherwise
!!  the ray power deposition is not correct.
!!
!! ARGUMENTS
!!
!!  Z  : Average ionization number
!!  e  : electron charge            (in esu   = g^(1/2) cm^(3/2) / s)
!!  k  : Boltzmann constant         (in erg/K = g cm^2 / s^2 K)
!!  T  : electron Temperature       (in K)
!!  Ne : electron density           (in cm^-3)
!!
!! NOTES
!!
!!***

real function ed_CoulombFactor (Z,e,k,T,Ne)

  implicit none

#include "constants.h"

  real, intent (in) :: Z
  real, intent (in) :: e
  real, intent (in) :: k
  real, intent (in) :: T
  real, intent (in) :: Ne

  real :: dummy1
  real :: dummy2
  real :: dummy3
  real :: dummy4
  real :: dummy5
!
!
!     ...Perform dummy operation to avoid compiler warning messages of unused variables.
!
!
  dummy1 = Z
  dummy2 = e
  dummy3 = k
  dummy4 = T
  dummy5 = Ne
!
!
!     ...Set Coulomb factor equal to 1.
!
!
  ed_CoulombFactor = 1.0
!
!
!     ...Ready!
!
!
  return
end function ed_CoulombFactor
