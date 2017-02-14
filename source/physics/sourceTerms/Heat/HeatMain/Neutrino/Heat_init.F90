!!****if* source/physics/sourceTerms/Heat/HeatMain/Neutrino/Heat_init
!!
!! NAME
!!  
!!  Heat_init
!!
!!
!! SYNOPSIS
!! 
!!  call Heat_init()
!!
!!  
!! DESCRIPTION
!!
!!  Perform various initializations (apart from the problem-dependent ones)
!!  for the heat module.
!!
!!
!! ARGUMENTS
!!
!!   
!!
!!***
subroutine Heat_init()

  use Heat_data

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype,Driver_getNumProcs

  implicit none

#include "constants.h"

  ! Everybody should know these

  call Driver_getMype(MESH_COMM,ht_meshMe)
  call Driver_getNumProcs(MESH_COMM, ht_numProcs)

  call RuntimeParameters_get("useHeat", useHeat)
  call RuntimeParameters_get("Lneut", ht_Lneut)
  call RuntimeParameters_get("Tneut", ht_Tneut)
  call RuntimeParameters_get("postBounce", ht_postBounce)
  call RuntimeParameters_get("useEntr", ht_useEntr)
  call RuntimeParameters_get("heatTimeFac", ht_heatTimeFac)

  return
end subroutine Heat_init
  
