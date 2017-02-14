!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonUtilities/pi_openMonitorFile
!!
!! NAME
!!
!!  pi_openMonitorFile
!!
!! SYNOPSIS
!!
!!  call pi_openMonitorFile ()
!!
!! DESCRIPTION
!!
!!  Opens the monitor file for recording proton imaging progress. The file name
!!  has already been set during initialization.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!          
!!  Only the master processor opens the file.
!!
!!***

subroutine pi_openMonitorFile ()

  use ProtonImaging_data, ONLY : pi_globalMe,        &
                                 pi_monitorFileName, &
                                 pi_monitorFileUnit

  implicit none

#include "constants.h"
#include "Flash.h"

  logical, save :: firstCall = .true.

  integer :: ut_getFreeFileUnit
!
!
!   ...Open the file only on the master processor.
!
!
  if (pi_globalMe == MASTER_PE) then

      pi_monitorFileUnit = ut_getFreeFileUnit ()

      if (firstCall) then
          open  (unit   = pi_monitorFileUnit, &
                 file   = pi_monitorFileName, &
                 form   = 'FORMATTED',        &
                 status = 'UNKNOWN'           )
          firstCall = .false.
      else
          open  (unit     = pi_monitorFileUnit, &
                 file     = pi_monitorFileName, &
                 form     = 'FORMATTED',        &
                 status   = 'OLD',              &
                 position = 'APPEND'            )
      end if
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine pi_openMonitorFile
