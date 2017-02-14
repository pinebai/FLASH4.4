!!****if* source/physics/sourceTerms/Polytrope/PolytropeMain/Polytrope_data
!!
!! NAME
!!  Polytrope_data
!!
!! SYNOPSIS
!!  Polytrope_data()
!!
!! DESCRIPTION
!!  Stores the local data for Source Term: Polytrope
!!
!! PARAMETERS
!!
!!   These are the runtime parameters used in the Polytrope unit.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.
!!
!!    poly_usePolytrope [BOOLEAN]
!!        Switch to turn the polytropic eos on or off at runtime.
!!
!! WRITTEN BY
!!  Christoph Federrath 2007-2012
!!
!!***

Module Polytrope_data

  logical, save :: poly_usePolytrope
  real, save    :: polytropeKonst
  real, save    :: polytropeDens1, polytropeDens2, polytropeDens3, &
                 & polytropeDens4, polytropeDens5
  real, save    :: polytropeGamma1, polytropeGamma2, polytropeGamma3, &
                 & polytropeGamma4, polytropeGamma5

end Module Polytrope_data
