!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_azbar
!!
!! NAME
!!
!! pchem_azbar
!!
!! SYNOPSIS
!! 
!! pchem_azbar()
!!
!! DESCRIPTION
!!
!! routine azbar computes composition variables from the mass fractions
!!
!! Given the mass fractions in xmass(i), return the molar abundances ymass(i),
!! total nucleon charge, zbar, mean numcleon charge squared z2bar, and the
!! electron mole number ye
!!
!! NOTES
!! 
!! the output varibles are stored in data structure PrimordialChemistry_dataEOS (write)
!!
!!***

subroutine pchem_azbar()
   
  use PrimordialChemistry_dataEOS
  use PrimordialChemistry_data

  implicit none

#include "constants.h"
#include "Flash.h"


  !! local declarations
  integer       i
  real          zbarxx,z2barxx

  zbarxx  = 0.0e0
  z2barxx = 0.0e0
  ytot1   = 0.0e0

  do i=1,NSPECIES
     ymass(i) = xmass(i) * aioninv(i)
     zbarxx   = zbarxx + zion(i) * ymass(i)
     z2barxx  = z2barxx + zionsq(i) * ymass(i)
     ytot1    = ytot1 + ymass(i)
  enddo

  abar  = 1.0e0/ytot1
  zbar  = zbarxx * abar
  z2bar = z2barxx * abar
  ye    = zbar * ytot1

  return
end subroutine pchem_azbar
