!!****if* source/physics/materialProperties/Viscosity/ViscosityMain/Spitzer/Viscosity_init
!!
!! NAME
!!
!!  Viscosity_init
!!
!! SYNOPSIS
!!
!!  Viscosity_init() 
!!                  
!!                
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!  
!!  
!!  
!!
!!
!!
!!***

subroutine Viscosity_init()

  use Viscosity_data, ONLY: &
        visc_useViscosity, &
       viscTempLow, viscTempHigh, viscSuppressFactor
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
  

  

  call RuntimeParameters_get("useViscosity", visc_useViscosity)
  call RuntimeParameters_get("viscTempLow", viscTempLow)
  call RuntimeParameters_get("viscTempHigh", viscTempHigh)
  call RuntimeParameters_get("viscSuppressFactor", viscSuppressFactor)

end subroutine Viscosity_init
