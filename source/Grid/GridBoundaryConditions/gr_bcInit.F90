!!****if* source/Grid/GridBoundaryConditions/gr_bcInit
!!
!! NAME
!!
!!  gr_bcInit
!!
!! SYNOPSIS
!!
!!  gr_bcInit()
!!
!! DESCRIPTION
!!
!!  Initialize values for all data in the module gr_ptData,
!!  and allocate the scratch buffers
!!
!! ARGUMENTS
!!
!!
!!***

subroutine gr_bcInit()
  use gr_bcData, ONLY : gr_bcUseGravity, gr_bcEintSwitch

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

#include "Flash.h"

  implicit none 

  gr_bcUseGravity = .false.
#ifdef GRAVITY
  call RuntimeParameters_get("useGravity",gr_bcUseGravity)
#endif

  gr_bcEintSwitch = 0.0
#ifdef FLASH_EOS
  call RuntimeParameters_get("eintSwitch",gr_bcEintSwitch)
#endif

  call gr_bcHseInit()

  return
end subroutine gr_bcInit
