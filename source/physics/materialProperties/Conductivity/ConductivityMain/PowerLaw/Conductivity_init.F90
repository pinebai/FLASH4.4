!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/PowerLaw/Conductivity_init
!!
!! NAME
!!
!!  Conductivity_init
!!
!! SYNOPSIS
!!
!!  Conductivity_init()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!
!!
!!***

subroutine Conductivity_init

  use Conductivity_data, ONLY: cond_useConductivity,&
      cond_TemperatureExponent, cond_K0, cond_alpha, &
      cond_DensityExponent
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  


  call RuntimeParameters_get("useConductivity", cond_useConductivity)

  call RuntimeParameters_get("cond_TemperatureExponent", cond_TemperatureExponent)
  call RuntimeParameters_get("cond_K0", cond_K0)
  call RuntimeParameters_get("cond_DensityExponent", cond_DensityExponent)

  cond_alpha = cond_K0


end subroutine Conductivity_init

