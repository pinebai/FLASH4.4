!!****ih* source/Grid/GridMain/paramesh/paramesh4/gr_pmFlashHookData
!!
!! NAME
!!
!!  gr_pmFlashHookData
!!
!! SYNOPSIS
!!
!!  use gr_pmFlashHookData
!!
!! DESCRIPTION
!!
!!  Module for some data items used by FLASH to customize PARAMESH.
!!
!!***

! Modification history:
!     Created   October 2016  KW

module gr_pmFlashHookData

  implicit none
      
  logical,parameter :: gr_pmAlwaysFillFcGcAtDomainBC = .FALSE.

end module gr_pmFlashHookData
