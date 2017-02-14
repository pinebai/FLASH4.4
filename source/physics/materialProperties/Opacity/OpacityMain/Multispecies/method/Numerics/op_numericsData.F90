!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Numerics/op_numericsData
!!
!! NAME
!!
!!  op_numericsData
!!
!! SYNOPSIS
!!
!!  use op_numericsData
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for the numerics opacity section. This section
!!  defines all machine dependent numerical values to be used within the opacity unit.
!!  If also defines numerical constants like one, zero, etc., which are used all over
!!  the place.
!!
!!  Meaning of data (some will be set in the op_initNumerics routine):
!!
!!    op_Pi                         = value of pi
!!    op_Avogadro                   = Avogadro's constant in # of particles / mole
!!    op_Boltzmann                  = Boltzmann's constant in erg / K
!!    op_classicalElectronRadius    = classical relativistic Thompson scattering length in cm
!!    op_electronCharge             = electron charge in esu (cgs unit of charge = 1statC = erg^1/2 cm^1/2)
!!    op_electronRestMass           = electron rest mass in g
!!    op_electronRestMassEnergy     = electron rest mass energy (rest mass * c^2) in erg
!!    op_speedLight                 = the vacuum speed of light in cm/s
!!    op_KleinNishinaPrefactor      = 2*pi * (electron radius)^2 * (Avogadro number) in cm^2/(mole of electrons)
!!    op_eV2Kelvin                  = temperature conversion factor from eV to K
!!    op_eV2erg                     = energy conversion factor from eV to erg (cgs energy unit)
!!    op_erg2eV                     = energy conversion factor from erg to eV
!!    op_keV2erg                    = energy conversion factor from keV to erg (cgs energy unit)
!!    op_erg2keV                    = energy conversion factor from erg to keV
!!    op_eV2keV                     = energy conversion factor from eV to keV
!!
!!    op_eBaseLowestExponent        = the lowest -x for which e^(-x) is still representable on the machine
!!    op_eBaseLargestExponent       = the largest x for which e^(+x) is still representable on the machine
!!    op_eBaseIgnorePlusOneExponent = the lowest  x for which e^(+x) + 1 = e^(+x) within computer accuracy
!!
!!    op_smallestPositiveNumber     = the smallest positive number representable on the machine
!!  
!!***
module op_numericsData
  
  implicit none

  logical, save :: op_initializedNumerics = .false.

  real,    save :: op_Avogadro
  real,    save :: op_Boltzmann
  real,    save :: op_classicalElectronRadius
  real,    save :: op_electronCharge
  real,    save :: op_electronRestMass
  real,    save :: op_electronRestMassEnergy
  real,    save :: op_eBaseLowestExponent
  real,    save :: op_eBaseLargestExponent
  real,    save :: op_eBaseIgnorePlusOneExponent
  real,    save :: op_erg2eV
  real,    save :: op_erg2keV
  real,    save :: op_KleinNishinaPrefactor
  real,    save :: op_Pi
  real,    save :: op_speedLight
  real,    save :: op_smallestPositiveNumber

  real, parameter :: zero  = 0.0
  real, parameter :: half  = 0.5
  real, parameter :: one   = 1.0
  real, parameter :: two   = 2.0
  real, parameter :: three = 3.0
  real, parameter :: four  = 4.0
  real, parameter :: ten   = 10.0

  real, parameter :: op_eV2Kelvin = 11604.5221
  real, parameter :: op_eV2keV    = 1.E-3
  real, parameter :: op_eV2erg    = 1.6021766208E-12
  real, parameter :: op_keV2erg   = 1.6021766208E-9

end module op_numericsData
