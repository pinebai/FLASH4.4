!!****if* source/diagnostics/ProtonImaging/localAPI/pi_openDiskProtonFile
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

  implicit none

  character (len=3), intent (in) :: whichOne
  logical, optional, intent (in) :: saveOldRecords

  return
end subroutine pi_openDiskProtonFile
