!!****h* source/physics/materialProperties/Conductivity/Conductivity_interface
!! NAME
!!  Conductivity_interface
!!
!! SYNOPSIS
!!  Conductivity_interface()
!!
!! DESCRIPTION
!!  This is the header file for the Conductivity Unit that defines its
!!  public interfaces.
!!***

#include "Flash.h"

module Conductivity_interface
  interface

     subroutine Conductivity_init()
     
     end subroutine Conductivity_init
  end interface

  interface
     subroutine Conductivity_finalize()
     end subroutine Conductivity_finalize
  end interface
  
  interface Conductivity

     subroutine Conductivity(xtemp,xden,massfrac,isochoricCond,diff_coeff,component)

       real, intent(IN)   :: xtemp
       real, intent(IN)   :: xden
       real, intent(OUT)  :: diff_coeff
       real, intent(OUT)  :: isochoricCond
       real, intent(IN) :: massfrac(NSPECIES)
       integer, intent(IN) :: component

     end subroutine Conductivity

     subroutine Conductivity_fullState(solnVec,isochoricCond,diffCoeff,component)

       real,              intent(IN) :: solnVec(NUNK_VARS)
       real,    OPTIONAL, intent(OUT)  :: diffCoeff
       real,    OPTIONAL, intent(OUT)  :: isochoricCond
       integer, OPTIONAL, intent(IN) :: component

     end subroutine Conductivity_fullState
  end interface

end module Conductivity_interface
