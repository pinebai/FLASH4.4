!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Numerics/op_writeNumerics
!!
!! NAME
!!
!!  op_writeNumerics
!!
!! SYNOPSIS
!!
!!  call op_writeNumerics ()
!!
!! DESCRIPTION
!!
!!  Prints out all the numerical data for the opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_writeNumerics ()

  use op_numericsData, ONLY : op_Pi,                         &
                              op_Avogadro,                   &
                              op_Boltzmann,                  &
                              op_KleinNishinaPrefactor,      &
                              op_classicalElectronRadius,    &
                              op_electronCharge,             &
                              op_electronRestMass,           &
                              op_electronRestMassEnergy,     &
                              op_speedLight,                 &
                              op_eV2Kelvin,                  &
                              op_eV2erg,                     &
                              op_erg2eV,                     &
                              op_keV2erg,                    &
                              op_erg2keV,                    &
                              op_eBaseLowestExponent,        &
                              op_eBaseLargestExponent,       &
                              op_eBaseIgnorePlusOneExponent, &
                              op_smallestPositiveNumber

  implicit none

#include "Opacity.h"

  integer :: fileUnit, ut_getFreeFileUnit
!
!
!   ...Open the printout file. 
!
!
  fileUnit = ut_getFreeFileUnit ()

  open (unit = fileUnit, &
        file = "opacity_printout_numerics.txt", &
        form = 'formatted')
!
!
!   ...Print out the title. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) "   OPACITY NUMERICS PRINTOUT"
  write (fileUnit,*)
!
!
!   ...Print all numerical constants. 
!
!
  write (fileUnit,'(1X,A42,ES20.10)') '                            Value of pi = ',op_Pi
  write (fileUnit,'(1X,A42,ES20.10)') ' Avogadro constant (# particles / mole) = ',op_Avogadro
  write (fileUnit,'(1X,A42,ES20.10)') '           Boltzmann constant (erg / K) = ',op_Boltzmann
  write (fileUnit,'(1X,A42,ES20.10)') ' Klein-Nishina prefactor (cm^2/mole e-) = ',op_KleinNishinaPrefactor
  write (fileUnit,'(1X,A42,ES20.10)') '                  Electron charge (esu) = ',op_electronCharge
  write (fileUnit,'(1X,A42,ES20.10)') '                      Electron mass (g) = ',op_electronRestMass
  write (fileUnit,'(1X,A42,ES20.10)') '         Classical electron radius (cm) = ',op_classicalElectronRadius
  write (fileUnit,'(1X,A42,ES20.10)') '        Electron rest mass energy (erg) = ',op_electronRestMassEnergy
  write (fileUnit,'(1X,A42,ES20.10)') '                  Speed of light (cm/s) = ',op_speedLight
  write (fileUnit,'(1X,A42,F9.3)'   ) '  Temperature conversion factor eV -> K = ',op_eV2Kelvin
  write (fileUnit,'(1X,A42,ES20.10)') '     Energy conversion factor eV -> erg = ',op_eV2erg
  write (fileUnit,'(1X,A42,ES20.10)') '     Energy conversion factor erg -> eV = ',op_erg2eV
  write (fileUnit,'(1X,A42,ES20.10)') '    Energy conversion factor keV -> erg = ',op_keV2erg
  write (fileUnit,'(1X,A42,ES20.10)') '    Energy conversion factor erg -> keV = ',op_erg2keV
  write (fileUnit,'(1X,A42,ES20.10)') '    Smallest positive number on machine = ',op_smallestPositiveNumber
  write (fileUnit,'(1X,A42,F6.0)'   ) '                   Lowest -x for e^(-x) = ',op_eBaseLowestExponent
  write (fileUnit,'(1X,A42,F6.0)'   ) '                  Largest  x for e^(+x) = ',op_eBaseLargestExponent
  write (fileUnit,'(1X,A42,F6.0)'   ) ' Lowest  x for e^(+x) + 1 equals e^(+x) = ',op_eBaseIgnorePlusOneExponent
!
!
!   ...Close the printout file. 
!
!
  close (fileUnit)
!
!
!    ...Ready!
!
!
  return
end subroutine op_writeNumerics
