!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDisk/pi_closeDiskProtonFile
!!
!! NAME
!!
!!  pi_closeDiskProtonFile
!!
!! SYNOPSIS
!!
!!  call pi_closeDiskProtonFile (character (len=3), intent (in) :: whichOne)
!!
!! DESCRIPTION
!!
!!  Closes the old or new disk proton files.
!!
!! ARGUMENTS
!!
!!  whichOne : controls which (old or new) disk proton file is to be closed
!!
!! NOTES
!!          
!!  The old disk proton file is closed by all processors. The new disk proton file can only
!!  be closed by the master processor.
!!
!!***

subroutine pi_closeDiskProtonFile (whichOne)

  use ProtonImaging_data,  ONLY : pi_diskProtonNewFileID, &
                                  pi_diskProtonOldFileID, &
                                  pi_globalMe

  implicit none

#include "constants.h"
#include "Flash.h"

  character (len=3), intent (in) :: whichOne

  logical :: fileOpen
!
!
!   ...Close the wanted disk proton file.
!
!
  if (whichOne == 'old') then

      inquire (unit = pi_diskProtonOldFileID, opened = fileOpen)

      if (fileOpen) then
          close (pi_diskProtonOldFileID)
      end if

  else if (whichOne == 'new') then

      if (pi_globalMe == MASTER_PE) then

          inquire (unit = pi_diskProtonNewFileID, opened = fileOpen)

          if (fileOpen) then
              close (pi_diskProtonNewFileID, status = 'DELETE')
          end if
      end if

  else
      call Driver_abortFlash ("pi_closeDiskProtonFile: No old/new disk proton file specified!")
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine pi_closeDiskProtonFile
