!!****if* source/diagnostics/ProtonImaging/localAPI/pi_traceProtons3DRec
!!
!! NAME
!!
!!  pi_traceProtons3DRec
!!
!! SYNOPSIS
!!
!!  call pi_traceProtons3DRec ()
!!
!! DESCRIPTION
!!
!!  Processes all protons on the collection of blocks on the current processor for
!!  those geometries consisting formally of 3D rectangular grids (cartesian).
!!  On exit, each proton has either:
!!
!!            i)  reached a different (yet unknown) block
!!           ii)  has reached the domain boundary and exited.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!  The use of threading is possible for tracing the protons through each block. 
!!
!!***

subroutine pi_traceProtons3DRec ()

  implicit none

  return
end subroutine pi_traceProtons3DRec
