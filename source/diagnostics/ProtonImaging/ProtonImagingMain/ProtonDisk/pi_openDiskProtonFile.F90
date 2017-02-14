!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDisk/pi_openDiskProtonFile
!!
!! NAME
!!
!!  pi_openDiskProtonFile
!!
!! SYNOPSIS
!!
!!  call pi_openDiskProtonFile (character (len=3), intent (in) :: whichOne,
!!                              logical, optional, intent (in) :: saveOldRecords)
!!
!! DESCRIPTION
!!
!!  Opens unformatted old or new disk proton files for accumulating/reading disk protons.
!!  Both files are given a corresponding file ID, which will serve for future writing
!!  and reading reference to the file.
!!
!! ARGUMENTS
!!
!!  whichOne       : controls which (old or new) disk proton file is to be opened
!!  saveOldRecords : if present and true, the old disk proton file is to be appended
!!
!! NOTES
!!          
!!  The old disk proton file can be opened by all processors, which need to read from
!!  that file. The new disk proton file can only be opened by the master processor.
!!  If no old disk proton file is present the routine simply does nothing, unless
!!  the processor is the master processor, in which case it creates a new old disk
!!  proton file.
!!
!!  If the optional keyword 'saveOldRecords' is passed and it is set to true, then
!!  the old disk proton file is opened in append mode, i.e. all previous records
!!  are saved.
!!
!!***

subroutine pi_openDiskProtonFile (whichOne, saveOldRecords)

  use Driver_interface,    ONLY : Driver_abortFlash

  use ProtonImaging_data,  ONLY : pi_diskProtonNewFileID,       &
                                  pi_diskProtonNewFileName,     &
                                  pi_diskProtonNewFileNbatches, &
                                  pi_diskProtonOldFileID,       &
                                  pi_diskProtonOldFileName,     &
                                  pi_globalMe

  implicit none

#include "constants.h"
#include "Flash.h"

  character (len=3), intent (in) :: whichOne
  logical, optional, intent (in) :: saveOldRecords

  logical :: appendFile
  logical :: fileExists
  integer :: ut_getFreeFileUnit
!
!
!   ...Open the wanted disk proton file.
!
!
  if (whichOne == 'old') then

      inquire (file = pi_diskProtonOldFileName, exist = fileExists)

      if (fileExists) then
          pi_diskProtonOldFileID = ut_getFreeFileUnit ()

          appendFile = .false.
          if (present (saveOldRecords)) then
              appendFile = saveOldRecords
          end if

          if (appendFile) then
              open (unit     = pi_diskProtonOldFileID,   &
                    file     = pi_diskProtonOldFileName, &
                    form     = 'UNFORMATTED',            &
                    position = 'APPEND',                 &
                    status   = 'OLD')
          else
              open (unit   = pi_diskProtonOldFileID,   &
                    file   = pi_diskProtonOldFileName, &
                    form   = 'UNFORMATTED',            &
                    status = 'OLD')
          end if
      else
          if (pi_globalMe == MASTER_PE) then
              pi_diskProtonOldFileID = ut_getFreeFileUnit ()
              open (unit   = pi_diskProtonOldFileID,   &
                    file   = pi_diskProtonOldFileName, &
                    form   = 'UNFORMATTED',            &
                    status = 'NEW')
          end if
      end if

  else if (whichOne == 'new') then

      if (pi_globalMe == MASTER_PE) then
          pi_diskProtonNewFileID = ut_getFreeFileUnit ()
          open (unit   = pi_diskProtonNewFileID,   &
                file   = pi_diskProtonNewFileName, &
                form   = 'UNFORMATTED',            &
                status = 'NEW')
          pi_diskProtonNewFileNbatches = 0
      end if

  else
      call Driver_abortFlash ("pi_openDiskProtonFile: No old/new disk proton file specified!")
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine pi_openDiskProtonFile
