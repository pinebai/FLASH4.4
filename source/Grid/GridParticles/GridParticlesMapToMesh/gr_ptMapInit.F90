!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/gr_ptMapInit
!!
!! NAME
!!
!!  gr_ptMapInit
!!
!! SYNOPSIS
!!
!!  gr_ptMapInit()
!!
!! DESCRIPTION
!!
!!  Initialize values for all data in the module gr_ptMapData
!!
!! ARGUMENTS
!!
!!
!!***

subroutine gr_ptMapInit()

  use gr_ptMapData, ONLY : gr_ptSmearLen
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none 


  call RuntimeParameters_get("smearLen", gr_ptSmearLen)


  if (gr_ptSmearLen < 0) then
     call Driver_abortFlash("Variable smearLen must be at least 0")
  end if

end subroutine gr_ptMapInit
