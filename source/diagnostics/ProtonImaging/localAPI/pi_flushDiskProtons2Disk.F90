!!****if* source/diagnostics/ProtonImaging/localAPI/pi_flushDiskProtons2Disk
!!
!! NAME
!!
!!  pi_flushDiskProtons2Disk
!!
!! SYNOPSIS
!!
!!  call pi_flushDiskProtons2Disk ()
!!
!! DESCRIPTION
!!
!!  When calling this routine, the disk protons accumulated by each processor are all
!!  send to the master processor and processed for writing out to disk. The following
!!  steps are performed:
!!
!!     1) Gather at the master processor the info of how many disk protons each
!!        processor currently has. Calculate the offsets to prepare for storage
!!        of the disk protons on the master processor.
!!
!!     2) Gather all the disk protons on the master processor. The disk protons
!!        on the master remain in position, the others are appended according to the
!!        offsets calculated previously.
!!
!!     3) On the master processor, write out the complete set of disk protons to disk
!!        either by appending them to an already existing disk proton file or by
!!        creating a new disk proton file.
!!
!!     4) Once all disk protons have been written to file, reset the disk proton counter
!!        to 0 on all processors.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pi_flushDiskProtons2Disk ()

  implicit none

  return
end subroutine pi_flushDiskProtons2Disk
