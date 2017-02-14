!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/op_setAtomNames
!!
!! NAME
!!
!!  op_setAtomNames
!!
!! SYNOPSIS
!!
!!  call op_setAtomNames ()
!!
!! DESCRIPTION
!!
!!  This routine sets the atomic element names.
!!
!! ARGUMENTS
!!
!!***
subroutine op_setAtomNames ()

  use Opacity_data,      ONLY : op_atomName
  use Driver_interface,  ONLY : Driver_abortFlash

  implicit none

  logical :: goodShape
!
!
!   ...Check dimensions of the atomic element name array.
!
!
  goodShape = size (op_atomName) == 100

  if (.not.goodShape) then
       call Driver_abortFlash ('[op_setAtomNames] ERROR: array op_atomName has bad shape')
  end if
!
!
!   ...Set the atomic names.
!
!
  op_atomName (  1) = 'HYDROGEN'
  op_atomName (  2) = 'HELIUM'
  op_atomName (  3) = 'LITHIUM'
  op_atomName (  4) = 'BERRYLLIUM'
  op_atomName (  5) = 'BORON'
  op_atomName (  6) = 'CARBON'
  op_atomName (  7) = 'NITROGEN'
  op_atomName (  8) = 'OXYGEN'
  op_atomName (  9) = 'FLOURINE'
  op_atomName ( 10) = 'NEON'

  op_atomName ( 11) = 'SODIUM'
  op_atomName ( 12) = 'MAGNESIUM'
  op_atomName ( 13) = 'ALUMINIUM'
  op_atomName ( 14) = 'SILICON'
  op_atomName ( 15) = 'PHOSPHORUS'
  op_atomName ( 16) = 'SULFUR'
  op_atomName ( 17) = 'CHLORINE'
  op_atomName ( 18) = 'ARGON'
  op_atomName ( 19) = 'POTASSIUM'
  op_atomName ( 20) = 'CALCIUM'

  op_atomName ( 21) = 'SCANDIUM'
  op_atomName ( 22) = 'TITANIUM'
  op_atomName ( 23) = 'VANADIUM'
  op_atomName ( 24) = 'CHROMIUM'
  op_atomName ( 25) = 'MANGANESE'
  op_atomName ( 26) = 'IRON'
  op_atomName ( 27) = 'COBALT'
  op_atomName ( 28) = 'NICKEL'
  op_atomName ( 29) = 'COPPER'
  op_atomName ( 30) = 'ZINC'

  op_atomName ( 31) = 'GALLIUM'
  op_atomName ( 32) = 'GERMANIUM'
  op_atomName ( 33) = 'ARSENIC'
  op_atomName ( 34) = 'SELENIUM'
  op_atomName ( 35) = 'BROMINE'
  op_atomName ( 36) = 'KRYPTON'
  op_atomName ( 37) = 'RUBIDIUM'
  op_atomName ( 38) = 'STRONTIUM'
  op_atomName ( 39) = 'YTTRIUM'
  op_atomName ( 40) = 'ZIRCONIUM'

  op_atomName ( 41) = 'NIOBIUM'
  op_atomName ( 42) = 'MOLYBDENUM'
  op_atomName ( 43) = 'TECHNETIUM'
  op_atomName ( 44) = 'RUTHENIUM'
  op_atomName ( 45) = 'RHODIUM'
  op_atomName ( 46) = 'PALLADIUM'
  op_atomName ( 47) = 'SILVER'
  op_atomName ( 48) = 'CADMIUM'
  op_atomName ( 49) = 'INDIUM'
  op_atomName ( 50) = 'TIN'

  op_atomName ( 51) = 'ANTIMONY'
  op_atomName ( 52) = 'TELLURIUM'
  op_atomName ( 53) = 'IODINE'
  op_atomName ( 54) = 'XENON'
  op_atomName ( 55) = 'CESIUM'
  op_atomName ( 56) = 'BARIUM'
  op_atomName ( 57) = 'LANTHANUM'
  op_atomName ( 58) = 'CERIUM'
  op_atomName ( 59) = 'PRASEODYMIUM'
  op_atomName ( 60) = 'NEODYMIUM'

  op_atomName ( 61) = 'PROMETHIUM'
  op_atomName ( 62) = 'SAMARIUM'
  op_atomName ( 63) = 'EUROPIUM'
  op_atomName ( 64) = 'GADOLINIUM'
  op_atomName ( 65) = 'TERBIUM'
  op_atomName ( 66) = 'DYSPROSIUM'
  op_atomName ( 67) = 'HOLMIUM'
  op_atomName ( 68) = 'ERBIUM'
  op_atomName ( 69) = 'THULIUM'
  op_atomName ( 70) = 'YTTERBIUM'

  op_atomName ( 71) = 'LUTETIUM'
  op_atomName ( 72) = 'HAFNIUM'
  op_atomName ( 73) = 'TANTALUM'
  op_atomName ( 74) = 'TUNGSTEN'
  op_atomName ( 75) = 'RHENIUM'
  op_atomName ( 76) = 'OSMIUM'
  op_atomName ( 77) = 'IRIDIUM'
  op_atomName ( 78) = 'PLATINUM'
  op_atomName ( 79) = 'GOLD'
  op_atomName ( 80) = 'MERCURY'

  op_atomName ( 81) = 'THALLIUM'
  op_atomName ( 82) = 'LEAD'
  op_atomName ( 83) = 'BISMUTH'
  op_atomName ( 84) = 'POLONIUM'
  op_atomName ( 85) = 'ASTATINE'
  op_atomName ( 86) = 'RADON'
  op_atomName ( 87) = 'FRANCIUM'
  op_atomName ( 88) = 'RADIUM'
  op_atomName ( 89) = 'ACTINIUM'
  op_atomName ( 90) = 'THORIUM'

  op_atomName ( 91) = 'PROTACTINIUM'
  op_atomName ( 92) = 'URANIUM'
  op_atomName ( 93) = 'NEPTUNIUM'
  op_atomName ( 94) = 'PLUTONIUM'
  op_atomName ( 95) = 'AMERICIUM'
  op_atomName ( 96) = 'CURIUM'
  op_atomName ( 97) = 'BERKELIUM'
  op_atomName ( 98) = 'CALIFORNIUM'
  op_atomName ( 99) = 'EINSTEINIUM'
  op_atomName (100) = 'FERMIUM'
!
!
!   ...Ready! 
!
!
  return
end subroutine op_setAtomNames
