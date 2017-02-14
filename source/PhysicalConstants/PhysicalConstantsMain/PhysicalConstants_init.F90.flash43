!!****if* source/PhysicalConstants/PhysicalConstantsMain/PhysicalConstants_init
!!
!! NAME
!!  PhysicalConstants_init
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_init()
!!
!! DESCRIPTION
!!
!! This subroutine initializes the Physical Constants databases for
!! units and constants.  Must be called in the Simulation_init
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  pc_unitsBase   [default "CGS"] set the default system of units, either "CGS" or "MKS"
!!                 CGS:  centimeters, grams, seconds, charge=esu
!!                 MKS:  meters,  kilograms, seconds, charge=Coloumb
!!                     both systems have temperature in Kelvin
!!
!! NOTES
!!
!!   If you see a constant named "Newton", what is meant is probably not the SI unit
!!   of force which has dimensions LENGTH^1 * TIME^(-2) * MASS^1, but 
!!   Newton's gravitational constant G.
!!
!!   Using the "physical constant" "pi" is deprecated, use PI from "constants.h" instead.
!!***!            

!!  Newer data are available! See
!!Peter J. Mohr and Barry N. Taylor (January 2005). "CODATA recommended values of the
!!fundamental physical constants: 2002". Reviews of Modern Physics 77: 1-107. An in-depth 
!!discussion of how the CODATA constants were selected and determined.

subroutine PhysicalConstants_init()

  use PhysicalConstants_data, ONLY: pc_sizeConstant, pc_sizeUnit,         &
       &              pc_initialized, pc_globalMe, pc_globalNumProcs
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use pc_interface, ONLY: pc_checkCGSMKS, pc_addUnit, pc_addConstant

  implicit none 

#include "constants.h" 

  character(len=MAX_STRING_LENGTH)        :: cgsORmks, errorstring
  integer                                 :: isError
  
  ! Large constants that require t.xDyz notation should be declared
  ! here and then used as arguments to pc_addUnit subroutine.  This
  ! is because xlf compilers promote t.xDyz to quad-precision when
  ! given the compile option -qrealsize=8, but we want to represent
  ! the constant in double-precision.  The notation t.xDyz is
  ! selected over t.xEyz so that large constants can also be declared
  ! in an absoft compiled application.
  real, parameter :: mass_MFLY = 9.8847D45, mass_clMass = 1.9889225D48
  real :: c, sb

!-------------------------------------------------------------------------
  ! Everybody should know these
  call Driver_getMype(GLOBAL_COMM,pc_globalMe)


  !  check if initialization has already occurred
  if (pc_initialized) return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Set the length of the databases to zero
  pc_sizeConstant = 0
  pc_sizeUnit = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Determine whether the user wants default values returned in CGS or MKS
  call RuntimeParameters_get("pc_unitsBase",cgsORmks)
  call pc_checkCGSMKS(cgsORmks,isError)
  if (isError /= 0) then
     write(errorString,900)cgsORmks
     call Driver_abortFlash(errorString)
     return
  endif

  ! Physical constants are taken from the 1998 Review of Particle
  ! Properties, Eur. Phys. J. C 3, 1 (1998), which in turn takes
  ! most of its values from Cohen, E. R. and Taylor, B. N.,
  ! Rev. Mod. Phys. 59, 1121 (1987).

  ! Set up the units of measurement.
  !! Note!  If you add new units, you'll need to change PhysicalConstants_interface

  call pc_addUnit ("length",         "cm",   1.)
  call pc_addUnit ("time",           "s",    1.)
  call pc_addUnit ("mass",           "g",    1.)
  call pc_addUnit ("temperature",    "K",    1.)
  call pc_addUnit ("charge",         "esu",  1.)
  call pc_addUnit ("substance amount","mol",  1.)

  call pc_addUnit ("length",         "m",    1.E2)    !meter
  call pc_addUnit ("length",         "km",   1.E5)    !kilometer
  call pc_addUnit ("length",         "pc",   3.0856775807E18)
  call pc_addUnit ("length",         "kpc",  3.0856775807E21)
  call pc_addUnit ("length",         "Mpc",  3.0856775807E24)
  call pc_addUnit ("length",         "Gpc",  3.0856775807E27)
  call pc_addUnit ("length",         "Rsun", 6.96E10)
  call pc_addUnit ("length",         "AU",   1.49597870662E13)
  call pc_addUnit ("time",           "yr",   3.15569252E7) !year
  call pc_addUnit ("time",           "Myr",  3.15569252E13)
  call pc_addUnit ("time",           "Gyr",  3.15569252E16)
  call pc_addUnit ("mass",           "kg",   1.E3) !kilogram
  call pc_addUnit ("mass",           "Msun", 1.9889225E33)
!  call pc_addUnit ("mass",           "amu",  1.660540210E-24) !OLD
  call pc_addUnit ("mass",           "amu",  1.660538782E-24)
  call pc_addUnit ("temperature",    "eV",   11604.505)
  call pc_addUnit ("charge",         "C",    2.99792458E9)  ! Coulomb

  ! cosmology-friendly units (after Vincenzo Antonuccio-Delogu):
  !
  !                          length_unit = 1 Mpc
  !                          time_unit   = 2/(3H_0)
  !                          mass_unit   = 5.23e12 Msun
  !                          H_0         = 100 km/sec/Mpc

  call pc_addUnit ("length",         "LFLY",  3.0856775807E24)
  call pc_addUnit ("time",           "TFLY",  2.05759E17)
  call pc_addUnit ("mass",           "MFLY",  mass_MFLY)

  ! Units appropriate for simulations of galaxy clusters

  call pc_addUnit ("length",         "clLength", 3.0856775807E24) ! Mpc
  call pc_addUnit ("time",           "clTime",   3.15569252E16) ! Gyr
  call pc_addUnit ("mass",           "clMass",   mass_clMass) ! 10^15 Msun
  call pc_addUnit ("temperature",    "clTemp",   11604440.207109345) ! keV

  ! Set up the physical constants.
  !! Note!  If you add new constants, you'll need to change PhysicalConstants_interface
  ! Exponent ordering:  len, time, mass, temp, charge (, substance amount)
  ! MUST adhere to electrostatic CGS units

  ! Note that Newton here is NOT the SI unit of force 1 N = 1 kg * m * s^(-2).
  ! It is Newton's gravitational constant (G).
  call pc_addConstant ("Newton",              6.67428E-8,               & 
       &                                                  3., -2., -1.,           & 
       &                                                  0., 0.)
  c = 2.99792458E10
  call pc_addConstant ("speed of light",      c,              &  
       &                                                  1., -1., 0.,            & 
       &                                                  0., 0.)
  call pc_addConstant ("Planck",              6.62606896E-27,            & 
       &                                                  2., -1., 1.,            & 
       &                                                  0., 0.)
  call pc_addConstant ("electron charge",     4.80320427E-10,            & 
       &                                                  0., 0., 0.,             & 
       &                                                  0., 1.)
  call pc_addConstant ("electron mass",       9.10938215E-28,            & 
       &                                                  0., 0., 1.,             & 
       &                                                  0., 0.)
  call pc_addConstant ("proton mass",        1.672621637E-24,            & 
       &                                                  0., 0., 1.,             & 
       &                                                  0., 0.)
  call pc_addConstant ("fine-structure",      7.2973525376E-3,            & 
       &                                                  0., 0., 0.,             & 
       &                                                  0., 0.)
  call pc_addConstant ("Avogadro",            6.02214179E23,             &  
       &                                                  0., 0., 0.,             & 
       &                                                  0., 0., -1.0)
  call pc_addConstant ("Boltzmann",           1.3806504E-16,             & 
       &                                                  2., -2., 1.,            & 
       &                                                  -1., 0.)
  call pc_addConstant ("ideal gas constant",  8.3144725E7,             & 
       &                                                  2., -2., 1.,            & 
       &                                                  -1., 0., -1.0)
  call pc_addConstant ("Wien",                2.8977685E-1,              & 
       &                                                  1., 0., 0.,             & 
       &                                                  1., 0.)

  sb = 5.670400E-5
  call pc_addConstant ("Stefan-Boltzmann",    sb,               & 
       &                                                  0., -3., 1.,            & 
       &                                                  -4., 0.)
  call pc_addConstant ("Radiation Constant",  4.0*sb/c, &
       -1.0, -2.0, 1.0, -4.0, 0.0)

  ! Note that using the "physical constant" "pi" is deprecated - use PI from constants.h instead.
  call pc_addConstant ("pi",         3.141592653589793238E0,              & 
       &                                         0., 0., 0., 0., 0.)
  call pc_addConstant ("e",          2.718281828459045235E0,              & 
       &                                         0., 0., 0., 0., 0.)
  call pc_addConstant ("Euler",      5.77215664901532861E-1,              & 
       &                                         0., 0., 0., 0., 0.)
  !------------------------------------------------------------------------
  !  initialization done
  pc_initialized = .true.

  return
900 format("PhysicalConstants_init: pc_unitsBase ", A6, &
       &         " invalid, must be CGS/MKS")
end subroutine PhysicalConstants_init


       
