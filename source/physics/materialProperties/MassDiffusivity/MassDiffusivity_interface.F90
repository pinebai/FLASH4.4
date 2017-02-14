!!****h* source/physics/materialProperties/MassDiffusivity/MassDiffusivity_interface
!!
!! NAME
!!  MassDiffusivity_interface
!!
!! SYNOPSIS
!!  use MassDiffusivity_interface
!!
!! DESCRIPTION
!!  This is the header file for the MassDiffusivity Unit that defines its
!!  public interfaces.
!!***

module MassDiffusivity_interface
  implicit none

  interface
     subroutine MassDiffusivity_init()
       
     end subroutine MassDiffusivity_init
  end interface

  interface
     subroutine MassDiffusivity_finalize()
     end subroutine MassDiffusivity_finalize
  end interface
  
  interface

     subroutine MassDiffusivity(xtemp,xden,massfrac,diffusivity)
       implicit none
#include "Flash.h"  
       real,INTENT(in)    :: xtemp
       real,INTENT(in)    :: xden
       real,INTENT(in)    :: massfrac(NSPECIES)
       real,INTENT(out)   :: diffusivity
       
     end subroutine MassDiffusivity
  end interface

end module MassDiffusivity_interface
