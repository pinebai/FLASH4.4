!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Collective/gr_pfftGenMapHelper
!!
!! NAME
!!
!!  gr_pfftGenMapHelper
!!
!! SYNOPSIS
!!  
!!  gr_pfftGenMapHelper(integer, intent(IN) :: axis, &
!!                      integer, dimension(:), intent(IN) :: fragmentPtr, &
!!                      integer, intent(OUT) :: maxSingleProcData)
!!
!! DESCRIPTION
!!
!!  This subroutine is a helper subroutine for gr_pfftGenMap.
!!  It calculate a metadata quantity named PFFT_NUMBLKS in the map structure. 
!!  This describes the number of blocks moving to each processor.
!!  We also return maxSingleProcData, which is an integer variable used 
!!  to size our send/recv buffer containing the actual block data.
!!
!! ARGUMENTS
!!
!!  axis:  Specifies whether we wish to update the JAXIS or KAXIS
!!         send map data structure.
!!  fragmentPtr:  Array describing the next free space in the 
!!                send map data structure.
!!  maxSingleProcData:  The maximum data that is sent to a single
!!                      processor along the current axis.
!!
!!***

subroutine gr_pfftGenMapHelper(axis, fragmentPtr, maxSingleProcData)

  use gr_pfftData, ONLY : pfft_comm
  use gr_pfftReconfigData, ONLY : pfft_maxProcs, pfft_sendJMap, pfft_sendKMap
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Pfft.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN) :: axis
  integer, dimension(:), intent(IN) :: fragmentPtr
  integer, intent(OUT) :: maxSingleProcData

  integer, dimension(:,:), pointer :: sendMap
  integer :: maxLocalProcData, bsize, nblks, i, j, ierr

  if (axis == JAXIS) then
     sendMap => pfft_sendJMap
  else if (axis == KAXIS) then
     sendMap => pfft_sendKMap
  else
     call Driver_abortFlash("[gr_pfftGenMapHelper]: Invalid AXIS passed!")
  end if

  if (ubound(fragmentPtr,1) /= pfft_maxProcs) then
     call Driver_abortFlash("[gr_pfftGenMapHelper]: Unexpected fragmentPtr size!")
  end if

  maxLocalProcData = 0
  do i = 1,pfft_maxProcs
     nblks = ((fragmentPtr(i)-PFFT_MAPEND) / PFFT_MAPEND) !Exclude metadata.
     bSize = 0

     do j = 1,nblks
        bSize = bSize + sendMap((j*PFFT_MAPEND) + PFFT_BSIZE,i)
     end do

     sendMap(PFFT_NUMBLKS,i) = nblks
     maxLocalProcData = max(bsize,maxLocalProcData)
  end do


  call MPI_ALLREDUCE(maxLocalProcData,maxSingleProcData,1,FLASH_INTEGER,MPI_MAX,&
       pfft_comm(IAXIS),ierr)

  nullify(sendMap)

end subroutine gr_pfftGenMapHelper
