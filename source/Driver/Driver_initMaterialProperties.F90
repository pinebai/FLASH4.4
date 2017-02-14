!!****f* source/Driver/Driver_initMaterialProperties
!!
!! NAME
!!  Driver_initMaterialProperties
!!
!! SYNOPSIS
!!
!!   call Driver_initMaterialProperties()
!!
!! DESCRIPTION
!!
!!   Initializes material properties units by
!!   calling their respective initialization routines 
!!    
!! ARGUMENTS
!!
!!   
!!
!!***

subroutine Driver_initMaterialProperties()
  use MagneticResistivity_interface, ONLY : MagneticResistivity_init
  use Conductivity_interface, ONLY : Conductivity_init
  use MassDiffusivity_interface, ONLY : MassDiffusivity_init
  use Viscosity_interface, ONLY : Viscosity_init

  implicit none


  call Conductivity_init()
  call MassDiffusivity_init()
  call Viscosity_init()
  call MagneticResistivity_init()


end subroutine Driver_initMaterialProperties
