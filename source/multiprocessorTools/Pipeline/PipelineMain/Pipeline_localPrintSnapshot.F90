!!****if* source/multiprocessorTools/Pipeline/PipelineMain/Pipeline_localPrintSnapshot
!!
!! NAME
!!  
!!  Pipeline_localPrintSnapshot
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localPrintSnapshot (subroutine, intent (in) :: printItemDetails)
!!
!! DESCRIPTION
!!
!!  The routine will print a snapshot of the current local state of the pipeline.
!!  It prints all items in each of the buffer, sending and receiving channels in
!!  the pipeline on the local processor.
!!
!! ARGUMENTS
!!
!!  printItemDetails : the routine specifying the print statements for the items
!!
!! NOTES
!!
!!  The routine 'printItemDetails' must be defined outside the pipeline unit.
!!  It will contain the specifics of the items for which the pipeline was used.
!!
!!***

subroutine Pipeline_localPrintSnapshot (printItemDetails)

  use Pipeline_data,     ONLY : pl_itemBuf,     &
                                pl_itemCount,   &
                                pl_numChannels, &
                                pl_recvBuf,     &
                                pl_recvCount,   &
                                pl_recvRequest, &
                                pl_sendBuf,     &
                                pl_sendCount

  implicit none

  include "Flash_mpi.h"

  interface
    subroutine printItemDetails (item, itemDescription)
      real,              intent (in) :: item (:)
      character (len=*), intent (in) :: itemDescription
    end subroutine printItemDetails
  end interface

  character (len=100) :: itemDescription

  integer :: i, n
  integer :: sendCount, recvCount
!
!
!     ...Show all items in local item buffer.
!
!
  do i = 1, pl_itemCount
     write (itemDescription,'(a,i10)') 'itemBuf: item ',i
     call printItemDetails (pl_itemBuf (:,i), trim (itemDescription))
  end do
!
!
!     ...Show all items in local sending buffers for each channel. It should always be
!        OK to read what is in the sending buffers, even if the send is not completed,
!        because:
!
!         1) "The sender should not modify any part of the send buffer after a
!             nonblocking send operation is called, until the send completes."
!             [MPI-3 3.7.2]
!
!         2) "(the send operation itself leaves the content of the send buffer
!             unchanged)" [MPI-3 3.7.3]
!
!
  do n = 1, pl_numChannels
     sendCount = pl_sendCount (n)
     if (sendCount > 0) then
         do i = 1, sendCount
            write (itemDescription,'(2(a,i10))') 'sendBuf: channel ',n,' , item ',i
            call printItemDetails (pl_sendBuf (:,i,n), trim (itemDescription))
         end do
     end if
  end do
!
!
!     ...Show all items in local receiving buffers for each channel. It is only OK to
!        read what is in the receiving buffers after the receive completed,
!        because:
!
!         1) "The receiver should not access any part of the receive buffer
!             after a nonblocking receive operation is called, until the
!             receive completes." [MPI-3 3.7.2]
!
!
  do n = 1, pl_numChannels
     recvCount = pl_recvCount (n)
     if (recvCount > 0 .and. pl_recvRequest (n) == MPI_REQUEST_NULL) then
         do i = 1, recvCount
            write (itemDescription,'(2(a,i10))') 'recvBuf: channel ',n,' , item ',i
            call printItemDetails (pl_recvBuf (:,i,n), trim (itemDescription))
         end do
     end if
  end do
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_localPrintSnapshot
