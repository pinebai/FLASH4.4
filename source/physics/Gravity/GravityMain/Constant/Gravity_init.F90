!!****if* source/physics/Gravity/GravityMain/Constant/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!  
!! SYNOPSIS
!!
!!  Gravity_init()
!!
!! DESCRIPTION
!!
!!  This is initialization routine for constant gravity. 
!!  This implementation of gravity unit uses only two
!!  runtime Parameters, as described in Gravity_data
!!
!! ARGUMENTS
!!
!!  
!!
!!
!!***

#include "constants.h"
subroutine Gravity_init ()

  use Gravity_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs

  implicit none

  


!==============================================================================

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getNumProcs(MESH_COMM,grv_meshNumProcs)

  call RuntimeParameters_get('gconst', grv_const)
  call RuntimeParameters_get('gdirec', grv_direc)

  call RuntimeParameters_get("useGravity", useGravity)
  grv_vector(:) = 0.0

  if ( (grv_direc == "z") .or. (grv_direc == "Z") ) then
     grv_vector(3) = grv_const
        
  elseif ( (grv_direc == "y") .or. (grv_direc == "Y") ) then
     grv_vector(2) = grv_const
     
  else 

! x is default dir if gdirec is unintelligible
     grv_vector(1) = grv_const

  endif

!==============================================================================

  return

end subroutine Gravity_init
