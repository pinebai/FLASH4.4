!!****if* source/physics/materialProperties/Opacity/localAPI/op_calculateLowTempOpacities
!!
!! NAME
!!
!!  op_calculateLowTempOpacities
!!
!! SYNOPSIS
!!
!!  call op_calculateLowTempOpacities ()
!!
!! DESCRIPTION
!!
!!  This routine constructs the low temperature Planck and Rosseland energy group opacities for all
!!  temperature grid/species/energy group triple combinations and places them in the corresponding
!!  tables. Currently it considers the photoelectronic (PE) and Klein-Nishina (KN) cross sections only.
!!  The PE data is extracted from the Biggs & Lighthill tables. 
!!
!! ARGUMENTS
!!
!!***
subroutine op_calculateLowTempOpacities ()

  implicit none

  return
end subroutine op_calculateLowTempOpacities
