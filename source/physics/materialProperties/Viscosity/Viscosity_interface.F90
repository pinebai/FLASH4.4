!!****h* source/physics/materialProperties/Viscosity/Viscosity_interface
!!
!! NAME    
!!   Viscosity_interface
!!
!! SYNOPSIS
!!
!!  use Viscosity_interface, ONLY : Viscosity
!!
!! DESCRIPTION
!!
!! This is the header file for the Viscosity Unit that defines its
!! public interfaces.
!!***

#include "Flash.h"

module Viscosity_interface
  interface

     subroutine Viscosity_init()
        
     end subroutine Viscosity_init
  end interface

  interface
     subroutine Viscosity_finalize()
     end subroutine Viscosity_finalize
  end interface
  
  interface
     subroutine Viscosity(xtemp,xden,massfrac,viscDynamic,viscKinematic)
       real,INTENT(in)    :: xtemp
       real,INTENT(in)    :: xden
       real,INTENT(in)    :: massfrac(NSPECIES)
       real,INTENT(out)   :: viscDynamic, viscKinematic
     end subroutine Viscosity
  end interface

end module Viscosity_interface
