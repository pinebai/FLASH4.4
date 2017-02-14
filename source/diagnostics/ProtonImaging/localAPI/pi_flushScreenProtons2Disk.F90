!!****if* source/diagnostics/ProtonImaging/localAPI/pi_flushScreenProtons2Disk
!!
!! NAME
!!
!!  pi_flushScreenProtons2Disk
!!
!! SYNOPSIS
!!
!!  call pi_flushScreenProtons2Disk ()
!!
!! DESCRIPTION
!!
!!  When calling this routine, the screen protons accumulated by each processor are all
!!  send to the master processor and processed for writing out to disk. The following
!!  steps are performed:
!!
!!     1) Gather at the master processor the info of how many screen protons each
!!        processor currently has. Calculate the offsets to prepare for storage
!!        of the screen protons on the master processor.
!!
!!     2) Gather all the screen protons on the master processor. The screen protons
!!        on the master remain in position, the others are appended according to the
!!        offsets calculated previously.
!!
!!     3) On the master processor, sort the collection of all screen protons into
!!        buckets according to their detector number. There are as many buckets
!!        as are screen detectors. If a bucket is full, trigger writing to disk to
!!        the assigned detector printout file.
!!
!!     4) Once all screen protons have been processed and all buckets have been emptied,
!!        reset the screen proton counter to 0 on all processors and reset the bucket
!!        counters on the master processor to 0 as well.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!  The routine has been designed in such a way that it can be called multiple times
!!  during a proton imaging simulation. The screen protons written out will be simply
!!  appended to the existing printout detector files.
!!
!!***

subroutine pi_flushScreenProtons2Disk ()

  implicit none

  return
end subroutine pi_flushScreenProtons2Disk
