!!****f* source/Simulation/Simulation_freeUserArrays
!!
!! NAME
!!  Simulation_freeUserArrays
!!
!! SYNOPSIS
!!  Simulation_freeUserArrays()
!!
!! DESCRIPTION
!!  This is where the user should place code for a setup that needs to
!!  free memory created in Simulation_init that may have been used in
!!  other aspects of initialization but is no longer needed once
!!  the initialization stage is finished. 
!!
!! ARGUMENTS
!!
!!***
subroutine Simulation_freeUserArrays()

  implicit none

  return

end subroutine Simulation_freeUserArrays
