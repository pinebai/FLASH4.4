!!****h* source/physics/materialProperties/MagneticResistivity/MagneticResistivity_interface
!!
!! NAME
!!  MagneticResistivity_interface
!!
!! SYNOPSIS
!!  MagneticResistivity_interface()
!!
!! DESCRIPTION
!!  This is the header file for the MagneticResistivity Unit that defines its
!!  public interfaces.
!!
!!***

module MagneticResistivity_interface

  implicit none

#include "Flash.h"

  interface
     subroutine MagneticResistivity_init()
        
     end subroutine MagneticResistivity_init
  end interface

  interface
     subroutine MagneticResistivity_finalize()
     end subroutine MagneticResistivity_finalize
  end interface
  
  interface MagneticResistivity
  
     subroutine MagneticResistivity(temp,dens,xn,magResist)
       real, intent(IN)  :: temp, dens
       real, intent(IN), dimension(NSPECIES)  :: xn
       real, intent(OUT) :: magResist
     end subroutine MagneticResistivity

     subroutine MagneticResistivity_fullState(solnVec,resPar, resPerp)
       real, intent(in)  :: solnVec(NUNK_VARS)
       real, intent(out) :: resPar
       real, OPTIONAL, intent(out) :: resPerp
     end subroutine MagneticResistivity_fullState

 
  end interface

     
end module MagneticResistivity_interface
