!!****if* source/diagnostics/ProtonImaging/localAPI/pi_closeDiskProtonFile
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

  implicit none

  character (len=3), intent (in) :: whichOne

  return
end subroutine pi_closeDiskProtonFile
