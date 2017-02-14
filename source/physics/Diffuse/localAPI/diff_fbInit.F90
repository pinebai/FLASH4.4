!!****if* source/physics/Diffuse/localAPI/diff_fbInit
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
!!    diffusion_cutoff_density
!!        density below which we no longer diffuse
!!    thermal_diff_method
!!***
subroutine diff_fbInit()
  implicit none
end subroutine diff_fbInit
