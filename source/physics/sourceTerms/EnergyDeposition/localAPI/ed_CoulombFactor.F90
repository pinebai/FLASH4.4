!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_CoulombFactor
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

  real, intent (in) :: Z
  real, intent (in) :: e
  real, intent (in) :: k
  real, intent (in) :: T
  real, intent (in) :: Ne

  ed_CoulombFactor = 0.0

  return
end function ed_CoulombFactor
