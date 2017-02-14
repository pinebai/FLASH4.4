!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Numerics/op_initNumerics
!!
!! NAME
!!
!!  op_initNumerics
!!
!! SYNOPSIS
!!
!!  call op_initNumerics ()
!!
!! DESCRIPTION
!!
!!  Inititalizes the section of the numerics for the opacity unit. It is here where
!!  several constants are set to enable smooth integrations without dangers of
!!  computational over- or underflow.
!!
!! ARGUMENTS
!!
!!***
subroutine op_initNumerics ()

  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  use Opacity_data,                ONLY : op_writeOpacityInfo

  use op_numericsData,             ONLY : one,two,                       &
                                          op_initializedNumerics,        &
                                          op_Pi,                         &
                                          op_Avogadro,                   &
                                          op_Boltzmann,                  &
                                          op_KleinNishinaPrefactor,      &
                                          op_classicalElectronRadius,    &
                                          op_electronCharge,             &
                                          op_electronRestMass,           &
                                          op_electronRestMassEnergy,     &
                                          op_speedLight,                 &
                                          op_eV2erg,                     &
                                          op_erg2eV,                     &
                                          op_keV2erg,                    &
                                          op_erg2keV,                    &
                                          op_eBaseLowestExponent,        &
                                          op_eBaseLargestExponent,       &
                                          op_eBaseIgnorePlusOneExponent, &
                                          op_smallestPositiveNumber

  use op_interface,                ONLY : op_writeNumerics

  implicit none

  integer :: base
  integer :: decimalPrecision

  real    :: ln10,lnbase
!
!
!   ...Get the needed physical constants.
!
!
  call PhysicalConstants_get ("pi",               op_Pi)
  call PhysicalConstants_get ("Avogadro",         op_Avogadro)
  call PhysicalConstants_get ("Boltzmann",        op_Boltzmann)
  call PhysicalConstants_get ("speed of light",   op_speedLight)
  call PhysicalConstants_get ("electron charge",  op_electronCharge)
  call PhysicalConstants_get ("electron mass",    op_electronRestMass)
!
!
!   ...Set the following constants:
!
!       op_erg2keV                    = inverse of op_keV2erg
!       op_classicalElectronRadius    = classical relativistic Thompson scattering length
!       op_KleinNishinaPrefactor      = 2 * pi * (electron radius)^2 * (Avogadro number)
!       op_electronRestMassEnergy     = electron rest mass energy (rest mass * c^2)
!       op_eBaseLowestExponent        = the lowest -x for which e^(-x) is still representable on the machine
!       op_eBaseLargestExponent       = the largest x for which e^(+x) is still representable on the machine
!       op_eBaseIgnorePlusOneExponent = the lowest  x for which e^(+x) + 1 = e^(+x) within computer accuracy
!       op_smallestPositiveNumber     = the smallest representable positive number on the machine
!
!
  base             = radix (one)
  ln10             = log (real (10))
  lnbase           = log (real (base))
  decimalPrecision = precision (one)

  op_erg2eV                     = one / op_eV2erg
  op_erg2keV                    = one / op_keV2erg
  op_electronRestMassEnergy     = op_electronRestMass * op_speedLight ** 2
  op_classicalElectronRadius    = op_electronCharge ** 2 / op_electronRestMassEnergy
  op_KleinNishinaPrefactor      = two * op_Pi * (op_classicalElectronRadius ** 2) * op_Avogadro
  op_eBaseLowestExponent        = real (int (lnbase * real (minexponent (one)) ))
  op_eBaseLargestExponent       = real (int (lnbase * real (maxexponent (one)) ))
  op_eBaseIgnorePlusOneExponent = real (int (  ln10 * real (decimalPrecision)  )) + one   ! + 1 added for safety
  op_smallestPositiveNumber     = tiny (one)
!
!
!   ...Write out the numerical constants (if requested).
!
!
  if (op_writeOpacityInfo) then
      call op_writeNumerics ()
  end if
!
!
!   ...Set initialization status.
!
!
  op_initializedNumerics = .true.
!
!
!   ...Ready! 
!
!
  return
end subroutine op_initNumerics
