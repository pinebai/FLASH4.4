!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Collective/gr_pfftReconfigData
!!
!! NAME
!!  gr_pfftReconfigData
!!
!! SYNOPSIS
!!  use gr_pfftReconfigData
!!  This includes subunit scope data for the pfft data movement!!  
!!
!!***

module gr_pfftReconfigData
#include "constants.h"
#include "Flash.h"
  implicit none

  integer, save, allocatable, dimension(:,:), target :: &
       pfft_sendJMap,pfft_sendKMap, pfft_recvJMap, pfft_recvKMap
  integer, save, allocatable, dimension(:) :: pfft_fragmentPtr
  real, save, allocatable, dimension(:,:) :: pfft_sendBuf, pfft_recvBuf
  integer, save :: pfft_maxProcs, pfft_maxProcData

  !We have devised pfft_procLookup which is an array of pointers to a 
  !a small datatype consisting of 2 integers. It allows us to extract 
  !information about a PFFT processor and its neighbors directly 
  !(i.e. no extensive search required, as is the case with pure arrays).

  !Example usage:
  !pfft_procLookup(JAXIS) % procInfo(3) % globalStartGridPoint
  !corresponds to information for the 4th (it is zero based) processor 
  !along the JAXIS.
  type TProcInfo
     integer :: globalStartGridPoint
     integer :: globalEndGridPoint
  end type TProcInfo

  type TPtrProcInfo
     type(TProcInfo), dimension(:), pointer :: procInfo
  end type TPtrProcInfo

  type(TPtrProcInfo), save, dimension(1:NDIM) :: pfft_procLookup

end module gr_pfftReconfigData
