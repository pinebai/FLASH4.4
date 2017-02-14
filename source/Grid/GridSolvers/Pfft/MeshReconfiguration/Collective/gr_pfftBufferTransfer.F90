!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Collective/gr_pfftBufferTransfer
!!
!! NAME
!!
!!  gr_pfftBufferTransfer
!!
!! SYNOPSIS
!!
!!  gr_pfftBufferTransfer(integer(IN)   :: direction,
!!                        integer(IN) ::  axis,
!!                        real, dimension(:,:) :: buffer,
!!                        real, dimension(:) :: pfftArray)
!! 
!! DESCRIPTION
!!  
!!  This routine is used in two different ways depending on how 
!!  it is called:
!!    1.  It copies data from buffer (pfft_recvBuf in calling subroutine) to 
!!        pfftArray (pfft_inArray in calling subroutine).
!!    2.  It copies data from pfftArray (pfft_outArray in calling subroutine) to 
!!        buffer (pfft_sendBuf in calling subroutine).
!!  
!!  Whichever way the data is being copied, the recvMap (which points 
!!  to either pfft_recvJMap or pfft_recvKMap) describes the coordinates
!!  where the data should be moved.
!!
!!  The data inside buffer is laid out as buffer(pointID, proc),
!!  where, pointID is a plain integer identifier.  The data can be associated with both
!!  a FLASH and a PFFT grid point with the help of the module array, pfft_recv[J|K]Map.
!!  To obtain the grid point in FLASH space we use pointID as a counter in order to  
!!  restore the grid points in the same way that we packed the grid points.
!!  To obtain the grid point in PFFT space we convert the associated NDIM PFFT 
!!  coordinates in pfft_recv[J|K]Map into a 1D coordinate in pfftArray.
!!
!!
!! ARGUMENTS
!!
!!  direction - Specifies whether we copy from the buffer into pfftArray
!!              or vice versa.
!!  axis -      Specifies the map which should be used to describe the 
!!              location of the data.  Can be JAXIS or KAXIS.
!!  buffer -    2D array that is either a receive buffer (recvBuf) or a 
!!              send buffer (sendBuf).
!!  pfftArray - 1D array that is either the pfft input or output array.
!!
!!  
!!***
subroutine gr_pfftBufferTransfer(direction, axis, buffer, pfftArray)

#include "constants.h"
#include "Pfft.h"
#include "Flash.h"

  use gr_pfftData, ONLY : pfft_inLen, pfft_ndim, pfft_me, pfft_myPE
  use gr_pfftReconfigData, ONLY : pfft_recvJMap, pfft_recvKMap, & 
       pfft_maxProcs, pfft_procLookup
  use ut_conversionInterface, ONLY : ut_convertToMemoryOffset
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, intent(IN) :: direction, axis
  real, dimension(:,:), intent(INOUT) :: buffer
  real, dimension(:), intent(INOUT) :: pfftArray

  integer, dimension(:,:), pointer :: recvMap
  integer, dimension(MDIM) :: elementCoord, arrayLBound
  integer :: i, j, k, eachDim
  integer :: np, nblks, ptr, bc, ri
  integer :: ib, ie, jb, je, kb, ke
  integer :: ioff, joff, koff, nprocs, oneDimensionalCoord, memoryOffset
  integer, dimension(1:MDIM) :: offset

  !Input options specify how to use the subroutine.
  !-----------------------------------------------------------------------
  if ( (direction /= TO_PFFT) .and. (direction /= FROM_PFFT) ) then
     call Driver_abortFlash("[gr_pfftBufferTransfer]: Direction not recognised!")
  end if
   
  !We always query the receive map regardless of the direction.
  if (axis == JAXIS) then
     recvMap => pfft_recvJMap
  else if (axis == KAXIS) then
     recvMap => pfft_recvKMap
  else
     call Driver_abortFlash("[gr_pfftBufferTransfer]: Axis must be JAXIS or KAXIS!")
  end if
  !-----------------------------------------------------------------------

  arrayLBound(:) = 1  !Required for ut_convertToMemoryOffset subroutine.
  nprocs = pfft_maxProcs

  offset(1:MDIM) = 0
  do eachDim = 1, pfft_ndim     
     offset(eachDim) = pfft_procLookup(eachDim) % procInfo(pfft_me(eachDim)) % globalStartGridPoint - 1
  end do


  ioff = offset(IAXIS)
  joff = offset(JAXIS)
  koff = offset(KAXIS)


  do np = 1, nProcs
     ri = 1
     ptr = PFFT_MAPEND
     nblks = recvMap(PFFT_NUMBLKS,np)

     do bc = 1, nblks
        ib = recvMap(ptr+PFFT_SPOSI,np) - ioff
        ie = recvMap(ptr+PFFT_EPOSI,np) - ioff
        jb = recvMap(ptr+PFFT_SPOSJ,np) - joff
        je = recvMap(ptr+PFFT_EPOSJ,np) - joff
        kb = recvMap(ptr+PFFT_SPOSK,np) - koff
        ke = recvMap(ptr+PFFT_EPOSK,np) - koff

        do k = kb, ke
           do j = jb, je
              do i = ib, ie

                 elementCoord(IAXIS) = i; elementCoord(JAXIS) = j; elementCoord(KAXIS) = k;
                 if (any(elementCoord < 0)) then 
                    print *, "ERROR!!!! Processor:", pfft_myPE, "elementCoord:", elementCoord
                    call Driver_abortFlash("Mistake with global coordinates!")
                 end if

                 call ut_convertToMemoryOffset(MDIM, elementCoord, arrayLBound, pfft_inLen, memoryOffset)              
                 oneDimensionalCoord = memoryOffset + 1  !Fortran arrays begin at element 1.


                 if (direction == TO_PFFT) then
                    pfftArray(oneDimensionalCoord) = buffer(ri,np)
                 else if (direction == FROM_PFFT) then
                    buffer(ri,np) = pfftArray(oneDimensionalCoord)
                 end if

                 ri = ri + 1
              end do
           end do
        end do
        ptr = ptr + PFFT_MAPEND
     end do
  end do

  nullify(recvMap)
  return

end subroutine gr_pfftBufferTransfer
