!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Constant-diff/Conductivity_init
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
!! Initializes the data in Conductivity data to a constant value "diff_constant"
!! from RuntimeParameters
!!
!! ARGUMENTS
!!
!!
!!
!!
!!***

subroutine Conductivity_init()

  use Conductivity_data, ONLY : cond_diffConstant, cond_useConductivity
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
  


  ! Everybody should know this

  
  call RuntimeParameters_get("diff_constant", cond_diffConstant)
  call RuntimeParameters_get("useConductivity", cond_useConductivity)

end subroutine Conductivity_init
