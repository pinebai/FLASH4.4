!!****if* source/Grid/GridParticles/gr_ptFinalize
!!
!! NAME
!!
!!  gr_ptFinalize
!!
!! SYNOPSIS
!!
!!  gr_ptFinalize()
!!
!! DESCRIPTION
!!
!!  deallocate the scratch buffers used in moving particles data
!!
!! ARGUMENTS
!!
!!
!!***


subroutine gr_ptFinalize()
  use Logfile_interface, ONLY : Logfile_stamp
  use gr_ptData, ONLY :   gr_ptDestBuf,gr_ptSourceBuf
  use Grid_data, ONLY :   gr_meshMe, gr_useParticles
  implicit none

  integer :: istat

  if (allocated(gr_ptDestBuf) .OR. allocated(gr_ptSourceBuf) .OR. gr_useParticles) then
     ! These are allocated in gr_ptInit
     deallocate(gr_ptDestBuf,stat=istat)
     if (istat .NE. 0) then
        call Logfile_stamp('WARNING - failed to deallocate gr_ptDestBuf','[gr_ptFinalize]')
     end if
     deallocate(gr_ptSourceBuf,stat=istat)
     if (istat .NE. 0) then
        call Logfile_stamp('WARNING - failed to deallocate gr_ptSourceBuf','[gr_ptFinalize]')
     end if
  end if

  return
end subroutine gr_ptFinalize
