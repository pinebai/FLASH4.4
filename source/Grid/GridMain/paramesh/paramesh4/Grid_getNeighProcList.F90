!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_getNeighProcList
!!
!! NAME
!!  Grid_getNeighProcList
!!
!! SYNOPSIS
!!
!!  call Grid_getNeighProcList(logical(IN)             :: includeMyProc, 
!!                             integer(INOUT), pointer :: neighProcList(:),
!!                             integer(OUT)            :: numNeigh)
!!
!! DESCRIPTION 
!!
!!  Creates a pointer array containing the neighboring processor IDs of all
!!  LEAF blocks on this processor.
!!
!! ARGUMENTS 
!!  
!!
!!  includeMyProc - Whether the array should include my processor ID.
!!
!!  neighProcList - The processor IDs of all neighboring LEAF blocks. 
!!
!!  numNeigh      - The number of entries in neighProcList.
!!
!!
!! NOTES
!!
!!  Currently only implemented for Paramesh4 Grid implementations.
!!
!!  It is the users resposibility to
!!    1. deallocate neighProcList
!!    2. obtain a new neighProcList when the mesh changes
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_getNeighProcList(includeMyProc, neighProcList, numNeigh)
  use Grid_data, ONLY : gr_meshMe
  use Grid_interface, ONLY : Grid_getListOfBlocks
  use ut_qsortInterface, ONLY : ut_qsort
  use gr_interface, ONLY : gr_findAllNeghID
  use gr_interfaceTypeDecl, ONLY: AllBlockRegions_t
  implicit none
  logical, intent(IN) :: includeMyProc
  integer, dimension(:), pointer :: neighProcList
  integer, intent(OUT) :: numNeigh
  integer, allocatable, dimension(:) :: procs
  integer, dimension(MAXBLOCKS) :: blockList
  integer :: i, j, k, b, p, n, ENDI, ENDJ, ENDK, allCenters, numNegh, &
       procCount, procSize, blockCount, blockID, procID, currentProc
  type(AllBlockRegions_t) :: surrBlksSummary

  nullify(neighProcList)
  numNeigh = 0

  ENDI = 3
  ENDJ = max(1,K2D*3)
  ENDK = max(1,K3D*3)
  allCenters = 2**NDIM
  call Grid_getListOfBlocks(LEAF, blockList, blockCount)

  !Up to "2**(NDIM-1)" neighbors for each guard cell region
  !Exactly "((ENDI*ENDJ*ENDK)-1)" guard cell regions in each block
  !Exactly "blockCount" LEAF blocks in this MPI rank.
  procSize = (2**(NDIM-1) * ((ENDI*ENDJ*ENDK)-1) * blockCount)
  allocate(procs(procSize))

  procCount = 0
  do b = 1, blockCount
     blockID = blockList(b)
     call gr_findAllNeghID(blockID, surrBlksSummary)
     do k = 1, ENDK
        do j = 1, ENDJ
           do i = 1, ENDI
              if ((i*j*k) /= allCenters) then
                 numNegh = surrBlksSummary % regionInfo(i,j,k) % numNegh
                 do n = 1, numNegh
                    procID = surrBlksSummary % regionInfo(i,j,k) % &
                         details(PROCNO,n)
                    if (procID /= gr_meshMe .or. includeMyProc) then
                       procCount = procCount + 1
                       procs(procCount) = procID
                    end if
                 end do
              end if
           end do
        end do
     end do
  end do

  if (procCount > 0) then
     call ut_qsort(procs, procCount)

     currentProc = procs(1)
     numNeigh = 1
     do i = 2, procCount !Replace duplicated processors with a -1
        if (currentProc == procs(i)) then
           procs(i) = -1
        else
           currentProc = procs(i)
           numNeigh = numNeigh + 1
        end if
     end do

     allocate(neighProcList(numNeigh))
     p = 1
     do i = 1, procCount !Copy the unique processors
        if (procs(i) /= -1) then
           neighProcList(p) = procs(i)
           p = p + 1
        end if
     end do
  end if
  deallocate(procs)
end subroutine Grid_getNeighProcList
