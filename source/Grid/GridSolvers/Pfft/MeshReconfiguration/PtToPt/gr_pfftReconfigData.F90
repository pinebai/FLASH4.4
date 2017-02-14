!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftReconfigData
!!
!! NAME
!!  gr_pfftReconfigData
!!
!! SYNOPSIS
!!  use gr_pfftReconfigData
!!
!!  This includes subunit scope data for the pfft data movement
!!  
!!***

module gr_pfftReconfigData
#include "constants.h"
#include "Flash.h"
  use gr_pfftListObject, ONLY : list
  implicit none

  !We have devised pfft_procLookup which is an array of pointers to a 
  !a small datatype consisting of 2 integers. It allows us to extract 
  !information about a PFFT processor and its neighbors directly 
  !(i.e. no extensive search required, as is the case with pure arrays).

  !Example usage:
  !pfft_procLookup(JAXIS) % procInfo(3) % globalStartGridPoint
  !corresponds to information for the 4th (it is zero based) 
  !processor along the JAXIS.
  type TProcInfo
     integer :: globalStartGridPoint
     integer :: globalEndGridPoint
  end type TProcInfo

  type TPtrProcInfo
     type(TProcInfo), dimension(:), pointer :: procInfo
  end type TPtrProcInfo

  type(TPtrProcInfo), save, dimension(1:NDIM) :: pfft_procLookup


  !----------------------- PtToPt specific ----------------------- 
  !Integers representing FLASH grid unk variables.
  integer, save :: pfft_srcGridVar, pfft_solnGridVar

  !Pointer to a particular 1D PFFT array:
  real, save, dimension(:), pointer :: pfft_pfftBuf

  !Pointers to maintain our lists:
  type(list), save, pointer :: pfft_listFG, pfft_listPG

  !Array holding number of messages to send to each proc:
  integer, save, allocatable, dimension(:) :: pfft_numMsgSendToEachProc

  integer, save :: pfft_pencilSize !Size of pencil assigned to pfft_myPE.
  integer, save :: pfft_numFlashNodes !Num. nodes in our FLASH list.
  integer, save :: pfft_numPfftNodes !Num. nodes in our PFFT list.
  integer, save :: pfft_metaType !MPI derived data type for node metadata.
  integer, save :: pfft_logUnit !Log file unit number (for debugging purposes).

end module gr_pfftReconfigData
