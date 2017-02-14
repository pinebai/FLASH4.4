!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/gr_pfftGridPointTable
!!
!! NAME 
!!
!! gr_pfftGridPointTable
!!
!! SYNOPSIS
!!
!! gr_pfftGridPointTable(integer(IN) :: pfft_inLen)
!!
!! DESCRIPTION 
!!  
!! Creates a lookup table which contains the grid points belonging 
!! to each processor.
!!
!! ARGUMENTS
!!
!!  pfft_inLen - the lengths along the dimensions
!!
!!
!!***
subroutine gr_pfftGridPointTable(pfft_inLen)
#include "constants.h"
#include "Flash.h"

  use gr_pfftData, ONLY : pfft_globalLen, pfft_ndim, pfft_procGrid
  use gr_pfftReconfigData, ONLY : pfft_procLookup
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, dimension(1:MDIM), intent(IN) :: pfft_inLen
  integer :: i, j, numPfftProcs, currentGridPoint, remainingGridPoints, error


  !The amount of data residing on each processor is given by pfft_inLen.  
  !However, if the global grid is not a multiple of pfft_inLen we have 
  !more PFFT grid points than FLASH grid points.  The 
  !data structure "pfft_procLookup" records which FLASH grid points 
  !actually reside on each PFFT processor.
  !-----------------------------------------------------------------  
  do i = 1, pfft_ndim
     numPfftProcs = pfft_procGrid(i)
     allocate(pfft_procLookup(i) % procInfo(0:numPfftProcs-1), STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("Severe error: Memory cannot be allocated!")
     end if

     currentGridPoint = 1
     j = 0

     !Iterate over every single PFFT processor along this dimension.
     do while (j < numPfftProcs)

        remainingGridPoints = pfft_globalLen(i) - currentGridPoint
        if (remainingGridPoints >= pfft_inLen(i)) then

           !We can assign a complete chunk to processor j.
           pfft_procLookup(i) % procInfo(j) % globalStartGridPoint = &
                currentGridPoint
           pfft_procLookup(i) % procInfo(j) % globalEndGridPoint = &
                currentGridPoint + pfft_inLen(i) - 1

        else

           !We do not have enough grid points to assign a whole chunk to processor j.  
           !Therefore, we will just give processor j whatever is left.         
           pfft_procLookup(i) % procInfo(j) % globalStartGridPoint = &
                currentGridPoint
           pfft_procLookup(i) % procInfo(j) % globalEndGridPoint = &
                currentGridPoint + remainingGridPoints

           !If there are any more processors, we just give them zero sized chunks.
           do while (j < numPfftProcs - 1)
              j = j + 1
              pfft_procLookup(i) % procInfo(j) % globalStartGridPoint = -1
              pfft_procLookup(i) % procInfo(j) % globalEndGridPoint = -1
           end do

        end if

        currentGridPoint = currentGridPoint + pfft_inLen(i)
        j = j + 1
     end do

  end do
  !call PrintPfftProcGrid()  !For checking.

end subroutine gr_pfftGridPointTable



!Print information about how the grid points are distributed amongst PFFT processors.
subroutine PrintPfftProcGrid()

  use gr_pfftData, ONLY : pfft_ndim, pfft_procGrid, pfft_myPE
  use gr_pfftReconfigData, ONLY : pfft_procLookup
  implicit none
  integer :: i, j

  if (pfft_myPE == 0) then
     print *, "Pfft processor grid shape:", pfft_procGrid
  end if

  do i = 1, pfft_ndim
     if (pfft_myPE == 0) then
        print *, "----- Dimension:", i
     end if
     do j = 0, ubound(pfft_procLookup(i) % procInfo,1)
        if (pfft_myPE == 0) then
           print *, "startPoint:", pfft_procLookup(i) % procInfo(j) % globalStartGridPoint, &                
                "endPoint:", pfft_procLookup(i) % procInfo(j) % globalEndGridPoint
        end if
     end do
  end do

end subroutine PrintPfftProcGrid
