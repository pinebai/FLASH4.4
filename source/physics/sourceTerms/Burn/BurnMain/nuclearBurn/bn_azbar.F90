!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bn_azbar
!!
!! NAME
!! 
!! bn_azbar
!!
!! SYNOPSIS
!!
!! bn_azbar()
!!
!! DESCRIPTION
!!
!!  routine azbar computes composition variables from the mass fractions 
!!  
!!  given the mass fractions in xmass(i), return the molar abundances ymass(i), 
!!  total number of moles per gram ytot1, the mean number of nucleons abar, 
!!  mean nucleon charge zbar, mean nucleon charge squared z2bar, and the
!!  electron mole number ye.
!!
!! NOTES
!!   the output variables are stored in data structure Burn_dataEOS
!!  
!!***

subroutine bn_azbar()

  use Burn_dataEOS, ONLY:  abar,zbar,z2bar,ytot1,ye
  use Burn_data, ONLY: xmass,ymass,aion,zion,bion,aioninv,zionsq

  implicit none

#include "constants.h"
#include "Flash.h"

  !!  local declarations
  integer          i
  real             zbarxx,z2barxx

  zbarxx  = 0.0e0
  z2barxx = 0.0e0
  ytot1   = 0.0e0

  do i=1,NSPECIES
     ymass(i) = xmass(i) * aioninv(i)
     zbarxx   = zbarxx + zion(i) * ymass(i)
     z2barxx  = z2barxx + zionsq(i) * ymass(i)
     ytot1    = ytot1 + ymass(i)
  enddo

  abar   = 1.0e0/ytot1
  zbar   = zbarxx * abar
  z2bar  = z2barxx * abar
  ye     = zbar * ytot1

  return
end subroutine bn_azbar
