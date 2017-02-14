!!****if* source/physics/materialProperties/Opacity/localAPI/op_computeIonNumberDensities
!!
!! NAME
!!
!!  op_computeIonNumberDensities
!!
!! SYNOPSIS
!!
!!  call op_computeIonNumberDensities ()
!!
!! DESCRIPTION
!!
!!  Computes the ion number densities for all the species in the current cell
!!  and the total sum over all species. The individual ion number densities for
!!  each species is given by the following expression:
!!
!!          # ion / cm^3  =  Mf * rho * Na / A
!!
!!  where Mf and A are the mass fraction and atomic weight of the species and
!!  rho and Na are the total cell mass density and Avogadro's number.
!!
!! ARGUMENTS
!!
!!***
subroutine op_computeIonNumberDensities ()

  implicit none

  return
end subroutine op_computeIonNumberDensities
