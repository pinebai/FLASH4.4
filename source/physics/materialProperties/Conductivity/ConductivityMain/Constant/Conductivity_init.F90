!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Constant/Conductivity_init
!!
!! NAME
!!
!!  Conductivity_init
!!
!! SYNOPSIS
!!
!!  call Conductivity_init()
!!                          
!!                         
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

subroutine Conductivity_init

  use Conductivity_data, ONLY : cond_constantIsochoric, cond_useConductivity
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  

  ! Everybody should know this

  call RuntimeParameters_get("cond_constantIsochoric", cond_constantIsochoric)
  call RuntimeParameters_get("useConductivity", cond_useConductivity)

end subroutine Conductivity_init
