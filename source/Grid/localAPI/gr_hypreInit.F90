!!****if* source/Grid/localAPI/gr_hypreInit
!!
!! NAME
!!
!!  gr_hypreInit
!!
!!
!! SYNOPSIS
!!
!!  call gr_hypreInit()
!!
!! Description
!!
!!  Initializes local data for Unit Diffuse defined in Module diff_saData.
!!  All the variables here are initialized by calling the
!!  RuntimeParameters_get subroutine. These data variables are for
!!  Unit Scope ->  Diffuse.
!!
!! ARGUMENTS
!!
!!  none  
!!
!! PARAMETERS
!!
!!    diff_scaleFactThermSaTempDiff
!!        factor by which the solution is scaled.
!!    diff_scaleFactThermSaTime
!!        factor by which diffusion time step is scaled.
!!***

subroutine gr_hypreInit() 
  
  implicit none  
  
  
  return
  
end subroutine gr_hypreInit
