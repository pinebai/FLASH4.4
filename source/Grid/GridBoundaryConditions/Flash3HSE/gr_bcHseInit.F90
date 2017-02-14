!!****if* source/Grid/GridBoundaryConditions/Flash3HSE/gr_bcHseInit
!!
!! NAME
!!
!!  gr_bcHseInit
!!
!! SYNOPSIS
!!
!!  call gr_bcHseInit()
!!
!! DESCRIPTION
!!
!!  Initialization for a hydrostatic equilibrium (HSE) implementation
!!
!! ARGUMENTS
!!
!!
!!***

subroutine gr_bcHseInit()

  use gr_bcHseData, ONLY : gr_bcHseGravDirec, gr_bcHseDirection, gr_bcHseGravConst

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

#include "constants.h"
  implicit none 

  call RuntimeParameters_get('gconst', gr_bcHseGravConst)
  call RuntimeParameters_get("gdirec",gr_bcHseGravDirec)
  select case (gr_bcHseGravDirec)
  case('x','X')
     gr_bcHseDirection = IAXIS
  case('y','Y')
     gr_bcHseDirection = JAXIS
  case('z','Z')
     gr_bcHseDirection = KAXIS
  case default
     call Driver_abortFlash('Runtime parameter "gdirec" only allows "x", "y", and "z".')
  end select

  return
end subroutine gr_bcHseInit
