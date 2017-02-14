!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabZeroOutsideCounts
!!
!! NAME
!!
!!  eos_tabZeroOutsideCounts
!!
!! SYNOPSIS
!!
!!  call eos_tabZeroOutsideCounts()
!!
!! DESCRIPTION
!! 
!! 
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***

subroutine eos_tabZeroOutsideCounts()

  use eos_tabData,                ONLY : eos_tabAllDiag

  implicit none

#include "Flash.h"

  integer :: species

  do species = 1,NSPECIES
     eos_tabAllDiag(species)%highTempCount = 0
     eos_tabAllDiag(species)%highDensCount = 0
     eos_tabAllDiag(species)%highestTemp = -999.0
     eos_tabAllDiag(species)%highestDens = -999.0
     eos_tabAllDiag(species)%highTempVarsLookedUp(:) = .FALSE.
     eos_tabAllDiag(species)%highDensVarsLookedUp(:) = .FALSE.
  end do
  

end subroutine eos_tabZeroOutsideCounts
