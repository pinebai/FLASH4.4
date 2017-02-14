!!****if* source/diagnostics/ProtonImaging/localAPI/pi_IOdetectorScreens
!!
!! NAME
!!
!!  pi_IOdetectorScreens
!!
!! SYNOPSIS
!!
!!  call pi_IOdetectorScreens ()
!!
!! DESCRIPTION
!!
!!  This routine stores IO protons that travel along the edges of each detector
!!  screen into the corresponding IO proton arrays. Five x,y,z positions need to
!!  be recorded in order for each detector IO proton to complete the screen
!!  perimeter. Each such detector IO proton is going to be assigned a tag that
!!  exceeds the current maximum tag value.
!!
!!***

subroutine pi_IOdetectorScreens ()
  
  implicit none

  return
end subroutine pi_IOdetectorScreens
