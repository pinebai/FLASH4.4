!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_updateRays
!!
!! NAME
!!
!!  ed_updateRays
!!
!! SYNOPSIS
!!
!!  call ed_updateRays (logical (in) :: doMove)
!!
!! DESCRIPTION
!!
!!  This routine updates the block ID info of the current set of rays.
!!  It consists of four stages:
!!
!!                   i) sort ray array on old block ID
!!                  ii) move ray array to new block ID
!!                 iii) sort ray array on new block ID
!!                  iv) extract new block ID info
!!
!!  When it is certain, that the old blockID's are still valid, then
!!  there is the option of calling this routine and having it execute
!!  only steps i) and iv). All nonexistent rays will be removed (i.e.
!!  placed at the very end of the rays array).
!!
!! ARGUMENTS
!!
!!  doMove : logical keyword to invoke movement of the rays (if false -> no move)
!!
!! NOTES
!!
!!***

subroutine ed_updateRays (doMove)

  implicit none

  logical, intent (in) :: doMove

  return
end subroutine ed_updateRays
