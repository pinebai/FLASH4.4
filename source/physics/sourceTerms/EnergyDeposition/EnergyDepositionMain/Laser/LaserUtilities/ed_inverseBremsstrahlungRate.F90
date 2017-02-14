!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_inverseBremsstrahlungRate
!!
!! NAME
!!
!!  ed_inverseBremsstrahlungRate
!!
!! SYNOPSIS
!!
!!  ed_inverseBremsstrahlungRate (real (in) :: Z,
!!                                real (in) :: e,
!!                                real (in) :: Me,
!!                                real (in) :: k,
!!                                real (in) :: T,
!!                                real (in) :: Ne,
!!                                real (in) :: Nc,
!!                                real (in) :: lnLambda)
!!
!! DESCRIPTION
!!
!!  Computes the inverse-Bremsstrahlung rate:
!!
!!            nu = (4/3) * sqrt (2pi/Me) * (Ze^4/Nc) * (Ne^2/(kT)^(3/2)) * lnLambda
!!
!!  The units of this rate are: # of electrons / s
!!
!! ARGUMENTS
!!
!!  Z        : Average ionization number
!!  e        : electron charge            (in esu   = g^(1/2) cm^(3/2) / s)
!!  Me       : electron mass              (in g)
!!  k        : Boltzmann constant         (in erg/K = g cm^2 / s^2 K)
!!  T        : electron Temperature       (in K)
!!  Ne       : electron density           (in cm^-3)
!!  Nc       : electron critical density  (in cm^-3)
!!  lnLambda : the Coulomb factor         (dimensionless)
!!
!! NOTES
!!
!!***

real function ed_inverseBremsstrahlungRate (Z,e,Me,k,T,Ne,Nc,lnLambda)

  implicit none

#include "constants.h"

  real, intent (in) :: Z
  real, intent (in) :: e
  real, intent (in) :: Me
  real, intent (in) :: k
  real, intent (in) :: T
  real, intent (in) :: Ne
  real, intent (in) :: Nc
  real, intent (in) :: lnLambda

  real :: e2Ne
  real :: kT
  real :: twoe2Ne
!
!
!     ...Calculate the rate.
!
!
  kT   = k * T
  e2Ne = e * e * Ne
  twoe2Ne = e2Ne + e2Ne

  ed_inverseBremsstrahlungRate =  (Z * twoe2Ne * twoe2Ne / (Nc + Nc + Nc)) &
                                 * sqrt ((PI + PI) / (Me * kT * kT * kT))    &
                                 * lnLambda
!
!
!     ...Ready!
!
!
  return
end function ed_inverseBremsstrahlungRate
