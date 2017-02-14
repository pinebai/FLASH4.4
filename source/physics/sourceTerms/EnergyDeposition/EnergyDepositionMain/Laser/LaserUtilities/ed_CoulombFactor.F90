!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_CoulombFactor
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
!!  Computes the Coulomb factor 'ln (Lambda)', where Lambda is given by:
!!
!!                 Lambda = (3/2Ze^3) * sqrt (2(kT)^3/(pi*Ne))
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

  real :: kT
  real :: Lambda
  real :: lnLambda
!
!
!     ...Calculate the Coulomb factor and floor it to 1.0 just in case.
!
!
  kT = k * T
  Lambda   = (1.5 / (Z * e * e * e)) * sqrt (kT * kT * kT / (PI * Ne))
  lnLambda = log (Lambda)

  ed_CoulombFactor = max (lnLambda , 1.0)
!
!
!     ...Ready!
!
!
  return
end function ed_CoulombFactor
