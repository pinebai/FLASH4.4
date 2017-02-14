!!****if* source/diagnostics/ProtonImaging/localAPI/pi_updateProtons
!!
!! NAME
!!
!!  pi_updateProtons
!!
!! SYNOPSIS
!!
!!  call pi_updateProtons (logical (in) :: doMove)
!!
!! DESCRIPTION
!!
!!  This routine updates the block ID info of the current set of protons.
!!  It consists of four stages:
!!
!!                   i) sort proton array on old block ID
!!                  ii) move proton array to new block ID
!!                 iii) sort proton array on new block ID
!!                  iv) extract new block ID info
!!
!!  When it is certain, that the old blockID's are still valid, then
!!  there is the option of calling this routine and having it execute
!!  only steps i) and iv). All nonexistent protons will be removed (i.e.
!!  placed at the very end of the protons array).
!!
!! ARGUMENTS
!!
!!  doMove : logical keyword to invoke movement of the protons (if false -> no move)
!!
!! NOTES
!!
!!***

subroutine pi_updateProtons (doMove)

  implicit none

  logical, intent (in) :: doMove

  return
end subroutine pi_updateProtons
