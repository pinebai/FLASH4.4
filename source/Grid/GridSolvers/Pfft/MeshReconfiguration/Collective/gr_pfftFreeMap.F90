!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/Collective/gr_pfftFreeMap
!!
!! NAME
!!
!!  gr_pfftFreeMap
!!
!! SYNOPSIS
!!
!!  gr_pfftFreeMap()
!!
!! DESCRIPTION
!!
!!  Finalises any variables, arrays used in Collective implementation.
!!
!! ARGUMENTS
!!
!!***

subroutine gr_pfftFreeMap()
  use gr_pfftData, ONLY : pfft_ndim
  use gr_pfftReconfigData, ONLY : pfft_sendJMap, &
       pfft_recvJMap, pfft_sendKMap, pfft_recvKMap, pfft_procLookup, &
       pfft_maxProcData, pfft_sendBuf, pfft_recvBuf
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer :: error, i

  if (pfft_ndim > 1) then
     deallocate(pfft_sendJMap, pfft_recvJMap, STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("[gr_pfftReconfigFinalise]:" // &
             "Severe error. JAXIS maps cannot be deallocated!")
     end if
  end if

  if (pfft_ndim > 2) then
     deallocate(pfft_sendKMap, pfft_recvKMap, STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("[gr_pfftReconfigFinalise]:" // &
             "Severe error. KAXIS maps cannot be deallocated!")
     end if
  end if

  do i = 1, pfft_ndim
     deallocate(pfft_procLookup(i) % procInfo, STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("[gr_pfftReconfigFinalise]:" // &
             "Severe error. pfft_procLookup cannot be deallocated!")
     end if
  end do

  if (pfft_maxProcData > 0) then
     deallocate(pfft_sendBuf, pfft_recvBuf, STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("[gr_pfftReconfigFinalise]:" // &
             "Severe error. send/recv buffers cannot be deallocated!")
     end if
  end if
end subroutine gr_pfftFreeMap
