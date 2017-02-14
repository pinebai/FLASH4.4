!!****if* source/physics/Eos/EosMain/Multigamma/eos_mgammaData
!!
!! NAME
!!
!!  eos_mgammaData
!!
!! 
!! SYNOPSIS
!!
!!  use eos_mgammaData
!!
!! DESCRIPTION
!!
!!  This is the data module for the Gamma law Eos implementation with 
!!  multiple fluids/species. 
!!  It stores all the runtime parameters, and all the unit scope
!!  data. Some of the unit scope data is fetched by the wrapper layer
!!  from elsewhere in the code and some is local unit data common to
!!  multiple functions in the unit.
!! 
!! PARAMETERS
!!  
!!   These are the runtime parameters used in Multigamma Eos.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.  
!!
!!   smalle[Real]  --- the smallest value for energy 
!!
!!***

module eos_mgammaData

#include "Flash.h"

  real, save, dimension(NSPECIES) ::  eos_gc, eos_gammam1j, eos_ggprodj, eos_ggprodinvj, eos_gam1invj

end module eos_mgammaData
