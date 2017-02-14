!!****if* source/physics/materialProperties/Viscosity/ViscosityMain/Constant/Viscosity_init
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
!!   Initializes viscosity data to a constant value in Viscosity_data.
!!
!! ARGUMENTS
!!
!!  
!!  
!!  
!!
!! PARAMETERS
!!
!!   visc_whichCoefficientIsConst which kind of coefficient to keep constant in Constant Viscosity implementation;
!!                                set to 1 for constant dynamic viscosity (the value of diff_visc_mu is used);
!!                                set to 2 for constant kinematic viscosity (the value of diff_visc_nu is used).
!!   diff_visc_mu  constant dynamic viscosity (used in Constant Viscosity if visc_whichCoefficientIsConst is 1)
!!   diff_visc_nu  constant kinematic viscosity (used in Constant Viscosity if visc_whichCoefficientIsConst is 2)
!!***

subroutine Viscosity_init()

  use Viscosity_data, ONLY : visc_useViscosity,visc_whichCoefficientIsConst,&
       visc_diffNu,visc_diffMu
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none


  call RuntimeParameters_get("useViscosity", visc_useViscosity)
  call RuntimeParameters_get("diff_visc_nu", visc_diffNu)
  call RuntimeParameters_get("diff_visc_Mu", visc_diffMu)
  call RuntimeParameters_get("visc_whichCoefficientIsConst", visc_whichCoefficientIsConst)

end subroutine Viscosity_init
