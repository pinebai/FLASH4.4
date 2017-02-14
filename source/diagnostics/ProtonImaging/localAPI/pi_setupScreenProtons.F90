!!****if* source/diagnostics/ProtonImaging/localAPI/pi_setupScreenProtons
!!
!! NAME
!!
!!  pi_setupScreenProtons
!!
!! SYNOPSIS
!!
!!  call pi_setupScreenProtons ()
!!
!! DESCRIPTION
!!
!!  Sets up the screen protons, which means to allocate the needed screen proton array.
!!  The screen protons are used to accumulate protons on the detector screen(s).
!!  The routine also sets up and commits a mpi type structure corresponding to the
!!  screen protons, which will be used when sending screen protons between processors.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pi_setupScreenProtons ()

  implicit none

  return
end subroutine pi_setupScreenProtons
