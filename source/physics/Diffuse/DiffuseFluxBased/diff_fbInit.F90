!!****if* source/physics/Diffuse/DiffuseFluxBased/diff_fbInit
!!
!! NAME
!!
!!  diff_fbInit
!!
!!
!! SYNOPSIS
!!
!!  call diff_fbInit()
!!
!! Description
!!
!!  Initializes local data for Unit Diffuse defined in Module Diffuse_data.
!!  All the variables here are initialized by calling the
!!  RuntimeParameters_get subroutine. These data variables are for
!!  Unit Scope ->  Diffuse, and used by all three major diffuse subroutines,
!!  Diffuse_species(for mass diffusivity), useDiffusetherm (for 
!!  thermal Conductivity), and useDiffuseVisc (for viscosity).
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!    useDiffuse
!!        whether any method of the Diffuse unit should contribute to fluxes
!!    useDiffuseTherm
!!        whether Diffuse_therm should contribute to fluxes
!!    useDiffuseVisc
!!        whether Diffuse_visc should contribute to fluxes
!!    useDiffuseSpecies
!!        whether Diffuse_therm should contribute to fluxes [TO BE IMPLEMENTED]
!!    geometry [STRING]
!!        Grid geometry
!!    dt_diff_factor
!!        factor that scales the timestep returned by Diffuse_computeDt
!!    diffusion_cutoff_density
!!        density below which we no longer diffuse
!!    thermal_diff_method
!!***


subroutine diff_fbInit()

  use Diffuse_data, ONLY: useDiffuse, thermal_diff_method, diff_geometricMeanDiff, &
       diff_scaleFactThermFlux
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Grid_interface, ONLY:   Grid_setFluxHandling

  implicit none

#include "constants.h"
#include "Flash.h"


  integer :: istat

  call RuntimeParameters_get('geometric_mean_diff', diff_geometricMeanDiff)
  call RuntimeParameters_get('thermal_diff_method', thermal_diff_method)
  if (.not. useDiffuse) then
     write(6,*)'WARNING:  You have included the flux-based Diffuse implementation '
     write(6,*)'   but have set the runtime parameter useDiffuse to FALSE'
     write(6,*)'   No Diffusion will occur but diff_fbInit will continue.'
  end if

  call RuntimeParameters_get('diff_scaleFactThermFlux',diff_scaleFactThermFlux)


  return
end subroutine diff_fbInit
