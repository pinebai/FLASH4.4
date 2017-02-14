!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/MagneticResistivity_fullState
!!
!! NAME
!!  MagneticResistivity_fullState
!!
!! SYNOPSIS
!!  call MagneticResistivity_fullState(real(in)    :: solnVec(NUNK_VARS),
!!                                     real(out)   :: resPar,
!!                            OPTIONAL,real(out) :: resPerp)
!!
!! DESCRIPTION
!!
!! Just returns on resPar the value from MagneticResistivity call. resPerp is ZERO 
!!
!! ARGUMENTS
!!
!!***

#include "Flash.h"

subroutine MagneticResistivity_fullState(solnVec,resPar, resPerp)
  
  use MagneticResistivity_interface, ONLY: MagneticResistivity
  implicit none

  real, intent(in)  :: solnVec(NUNK_VARS)
  real, intent(out) :: resPar
  real, OPTIONAL, intent(out) :: resPerp
  
  ! Stub
  resPar  = 0.0
  if (present(resPerp)) resPerp = 0.0
  
#if defined(DENS_VAR) && defined(TEMP_VAR)
  call MagneticResistivity(solnVec(TEMP_VAR),solnVec(DENS_VAR), &
       solnVec(SPECIES_BEGIN:SPECIES_END), resPar)
#endif
    
end subroutine MagneticResistivity_fullState

