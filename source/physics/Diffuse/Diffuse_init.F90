!!****f* source/physics/Diffuse/Diffuse_init
!!
!! NAME
!!
!!  Diffuse_init
!!
!! SYNOPSIS
!!
!!  call Diffuse_init()
!!
!! Description
!!
!!  Initializes local data for Unit Diffuse defined in Module Diffuse_data.
!!  All the variables here are initialized by calling the
!!  RuntimeParameters_get subroutine. These data variables are for
!!  Unit Scope ->  Diffuse, and used by all three major diffuse subroutines,
!!  Diffuse_species(for mass diffusivity), Diffuse_therm (for 
!!  thermal Conductivity), and Diffuse_visc (for viscosity).
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

!!  This is an empty subroutine (stub), an implementation is
!!  in DiffuseMain/Diffuse_init.

subroutine Diffuse_init()

    implicit none

  return
end subroutine Diffuse_init
