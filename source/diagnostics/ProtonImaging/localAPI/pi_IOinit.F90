!!****if* source/diagnostics/ProtonImaging/localAPI/pi_IOinit
!!
!! NAME
!!
!!  pi_IOinit
!!
!! SYNOPSIS
!!
!!  call pi_IOinit ()
!!
!! DESCRIPTION
!!
!!  This routine initializes variables for plotting proton data. As this needs
!!  the exact number of total protons that will be emitted, this routine is only
!!  successful, if all the beams have been set up. The number of emitted protons
!!  might change slightly during the beams setup, so we cannot use the requested
!!  proton emission number values from the flash.par.
!!
!!***

subroutine pi_IOinit ()
  
  implicit none

  return
end subroutine pi_IOinit
