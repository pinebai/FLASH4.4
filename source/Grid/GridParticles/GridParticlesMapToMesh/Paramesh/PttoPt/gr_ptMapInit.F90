!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/PttoPt/gr_ptMapInit
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

  use gr_ptMapData, ONLY : gr_ptSmearLen, gr_ptRecvSpecifier, &
       gr_ptRecvSpecifierTmp, gr_ptRecvTotalTmp, gr_ptRecvTotal
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_meshNumProcs

  implicit none 
  integer :: error

  call RuntimeParameters_get("smearLen", gr_ptSmearLen)

  if (gr_ptSmearLen < 0) then
     call Driver_abortFlash("[gr_ptMapInit]: Variable smearLen must be at least 0")
  end if


  allocate(gr_ptRecvSpecifier(0:gr_meshNumProcs-1), & 
       gr_ptRecvSpecifierTmp(0:gr_meshNumProcs-1), &
       gr_ptRecvTotalTmp(0:gr_meshNumProcs-1), &
       gr_ptRecvTotal(0:gr_meshNumProcs-1), STAT=error)
  if (error /= 0) then
     call Driver_abortFlash("[gr_ptMapInit]: Memory cannot be allocated!")
  end if

end subroutine gr_ptMapInit
