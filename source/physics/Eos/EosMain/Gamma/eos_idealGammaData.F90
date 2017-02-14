!!****if* source/physics/Eos/EosMain/Gamma/eos_idealGammaData
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

!  real, save :: eos_gamma
!  real, save :: eos_singleSpeciesA
!  real, save :: eos_singleSpeciesZ
  real, save :: eos_gammam1


end module eos_idealGammaData
