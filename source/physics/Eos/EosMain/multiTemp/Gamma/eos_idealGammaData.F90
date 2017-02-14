!!****if* source/physics/Eos/EosMain/multiTemp/Gamma/eos_idealGammaData
!!
!! NAME
!!
!!  eos_idealGammaData
!!
!! 
!! SYNOPSIS
!!
!! use eos_idealGammaData
!!
!! DESCRIPTION
!!
!!  This is the data module for the Gamma law Eos implementation.
!!  It stores all the runtime parameters, and all the unit scope
!!  data. Some of the unit scope data is fecthed by the wrapper layer
!!  from elsewhere in the code and some is local unit data common
!!  multiple functions in the unit 
!! 
!! PARAMETERS
!!  
!!   These are the runtime parameters used by the Gamma law implementation
!!   of the Eos unit.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.  
!!
!!   gamma[Real]   --- The ideal gas gamma from runtime parameters
!!   smalle[Real]  --- the smallest value for energy (runtime Paramters)
!!
!!***


module eos_idealGammaData

#include "Flash.h"

  real, save :: eos_gasConstant
  real, save :: eos_gamma, eos_gammaIon, eos_gammaEle, eos_gammaRad
  real, save :: eos_singleSpeciesA
  real, save :: eos_singleSpeciesZ
  real, save :: eos_gammam1, eos_gammam1Ion, eos_gammam1Ele, eos_gammam1Rad

  real, save :: eos_eMass
  real, save :: eos_eMassInUAmu !electron mass in unified atomic mass units (aka daltons)


end module eos_idealGammaData
