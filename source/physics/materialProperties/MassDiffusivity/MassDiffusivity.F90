!!****f* source/physics/materialProperties/MassDiffusivity/MassDiffusivity
!!
!!  NAME    
!!   MassDiffusivity
!!
!!  SYNOPSIS
!!   MassDiffusivity(real(in) :: xtemp, 
!!                    real(in) :: xden, 
!!                    real(in) :: massfrac, 
!!                    real(out) :: diffusivity)
!!
!!
!!  DESCRIPTION
!!   A generic mass diffusivity routine.  Just returns a constant 
!!   diffusivity from the run-time parameter, 'diff_spec_D'
!!  
!!  ARGUMENTS
!!   xtemp : temperature (in K)
!!   xden : density (in g/cm**3)
!!   massfrac : mass fractions of the composition
!!   diffusivity : mass diffusivity
!!
!!***

subroutine MassDiffusivity(xtemp,xden,massfrac,diffusivity)
!! True Stub  
  implicit none
#include "Flash.h"

  real,INTENT(in)    :: xtemp
  real,INTENT(in)    :: xden
  real,INTENT(in)    :: massfrac(NSPECIES)
  real,INTENT(out)   :: diffusivity

  !dummy value assigned in stub
  diffusivity = 0    
  
  return 
end subroutine MassDiffusivity


