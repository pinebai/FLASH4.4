!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonUtilities/pi_closeMonitorFile
!!
!! NAME
!!
!!  pi_closeMonitorFile
!!
!! SYNOPSIS
!!
!!  call pi_closeMonitorFile ()
!!
!! DESCRIPTION
!!
!!  Closes the monitor file for recording proton imaging progress.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!          
!!  Only the master processor closes the file.
!!
!!***

subroutine pi_closeMonitorFile ()

  use ProtonImaging_data, ONLY : pi_globalMe,       &
                                 pi_monitorFileUnit

  implicit none

#include "constants.h"
#include "Flash.h"
!
!
!   ...Close the file only on the master processor.
!
!
  if (pi_globalMe == MASTER_PE) then
      close (pi_monitorFileUnit)
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine pi_closeMonitorFile
