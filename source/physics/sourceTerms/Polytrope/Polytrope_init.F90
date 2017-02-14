!!****f* source/physics/sourceTerms/Polytrope/Polytrope_init
!!
!! NAME
!!  Polytrope_init
!!
!! SYNOPSIS
!!  Polytrope_init(integer(in) :: myPE)
!!
!! DESCRIPTION
!!  Implement the polytropic eos as source term
!!
!! ARGUMENTS
!!   myPE    - current processor number
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
!!    poly_polytropicGamma1 [REAL]
!!        polytropic exponent 1
!!    poly_usePolytrope [BOOLEAN]
!!        runtime switch to turn Polytrope on/off
!!
!! WRITTEN BY
!!   Christoph Federrath 2007
!!
!!***

subroutine Polytrope_init()
  implicit none
end subroutine Polytrope_init
