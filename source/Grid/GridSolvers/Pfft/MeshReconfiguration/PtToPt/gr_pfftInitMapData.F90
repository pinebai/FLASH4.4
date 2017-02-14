!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftInitMapData
!!
!! NAME 
!!
!! gr_pfftInitMapData
!!
!! SYNOPSIS
!!
!! gr_pfftInitMapData()
!!
!! DESCRIPTION 
!!
!! Intialises variables, arrays and data structure which will be used
!! when we map from FLASH grid -> Pencil grid and Pencil grid -> FLASH grid. 
!!
!! ARGUMENTS
!!
!! SIDE EFFECTS 
!!
!! Initialisation of module level variables, arrays, lists:
!!
!! NOTES
!! 
!! We construct an MPI datatype in this routine which describes the metadata 
!! that we exchange between nodes.
!!
!!***

subroutine gr_pfftInitMapData()
#include "constants.h"
  use Driver_interface, ONLY : Driver_checkMPIErrorCode
  use Logfile_interface, ONLY : Logfile_open
  use gr_pfftListObject, ONLY : initialise_list
  use gr_pfftReconfigData, ONLY : pfft_metaType, pfft_numFlashNodes, &
       pfft_numPfftNodes, pfft_srcGridVar, pfft_solnGridVar, pfft_pfftBuf, &
       pfft_numMsgSendToEachProc, pfft_listFG, pfft_listPG, pfft_logUnit
  use gr_pfftData, ONLY : pfft_numProcs

  implicit none
  include "Flash_mpi.h"
  integer, dimension(1) :: oldTypes, blockCounts, offsets
  integer :: ierr, typeSize, typeExtent

  pfft_metaType = MPI_DATATYPE_NULL
  pfft_numFlashNodes = 0; pfft_numPfftNodes = 0
  pfft_srcGridVar = NONEXISTENT; pfft_solnGridVar = NONEXISTENT
  nullify(pfft_pfftBuf)

  allocate( pfft_numMsgSendToEachProc(0:pfft_numProcs-1) )
  pfft_numMsgSendToEachProc(0:pfft_numProcs-1) = 0

  call initialise_list(pfft_listFG)  !A list for FLASH grid fragments.
  call initialise_list(pfft_listPG)  !A list for PFFT grid fragments.


  !Create an MPI datatype describing the metadata portion of a node.
  !There are 4 integers and 4 integer arrays of size MDIM stored contiguously:
  !i.e. node_metadata looks like:
  !
  !integer :: flashProcID, flashBlockID, pfftProcID, tagID
  !integer, dimension(1:MDIM) :: flashStartPos, flashEndPos, &
  !      pfftStartPos, pfftEndPos
  offsets(1) = 0
  oldtypes(1) = FLASH_INTEGER
  blockcounts(1) = 4 + (MDIM * 4)

  !Note: Some systems do not have an implementation for MPI_Type_create_struct() 
  !which is part of MPI-2.  As such, we use MPI_Type_struct() which is safe 
  !in this particular case because we never refer to the addresses of the struct 
  !elements.  We could not (safely) get away with this on a 64 bit system if our
  !struct consisted of different datatypes (here we must refer to addresses).
  call MPI_Type_struct(1, blockCounts, offsets, oldTypes, &
       pfft_metaType, ierr) 
  call Driver_checkMPIErrorCode(ierr)

  call MPI_Type_commit(pfft_metaType, ierr)
  call Driver_checkMPIErrorCode(ierr)

#ifdef DEBUG_PFFT
  call Logfile_open(pfft_logUnit,.true.)
#endif
end subroutine gr_pfftInitMapData
