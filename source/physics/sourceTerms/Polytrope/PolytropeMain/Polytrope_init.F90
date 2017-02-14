!!****if* source/physics/sourceTerms/Polytrope/PolytropeMain/Polytrope_init
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
!!    poly_usePolytrope [BOOLEAN]
!!        runtime switch to turn Polytrope on/off
!!
!! WRITTEN BY
!!   Christoph Federrath 2007
!!
!!***

subroutine Polytrope_init()
  use Polytrope_data
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype
  implicit none
#include "constants.h"
#include "Flash.h"

  integer :: mype
  
  call Driver_getMype(GLOBAL_COMM, mype)
  
  call RuntimeParameters_get('usePolytrope',   poly_usePolytrope)
  call RuntimeParameters_get('polytropeKonst', polytropeKonst)

  call RuntimeParameters_get('polytropeGamma1', polytropeGamma1)
  call RuntimeParameters_get('polytropeGamma2', polytropeGamma2)
  call RuntimeParameters_get('polytropeGamma3', polytropeGamma3)
  call RuntimeParameters_get('polytropeGamma4', polytropeGamma4)
  call RuntimeParameters_get('polytropeGamma5', polytropeGamma5)

  call RuntimeParameters_get('polytropeDens1', polytropeDens1)
  call RuntimeParameters_get('polytropeDens2', polytropeDens2) 
  call RuntimeParameters_get('polytropeDens3', polytropeDens3)
  call RuntimeParameters_get('polytropeDens4', polytropeDens4)
  call RuntimeParameters_get('polytropeDens5', polytropeDens5)

  if (poly_usePolytrope .and. (mype .eq. MASTER_PE)) then
    print*, "Initializing Polytropic Equation of State"
  endif
end subroutine Polytrope_init
