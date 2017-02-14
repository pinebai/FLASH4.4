!!****if* source/physics/Diffuse/DiffuseMain/diff_advanceTherm
!!
!!  NAME 
!!
!!  diff_advanceTherm
!!
!!  SYNOPSIS
!!
!!  call diff_advanceTherm(integer(IN)                  :: blockCount,
!!                         integer(IN)                  :: blockList(blockCount),
!!                         real(IN)                     :: dt,
!!                         integer, OPTIONAL, intent(IN):: pass)
!!
!!  DESCRIPTION 
!!      This routine advances the heat diffusion equation (i.e., heat conduction).
!!      An implicit scheme is used. 
!!
!!      This is a stub.
!!
!! ARGUMENTS
!!
!!  blockCount   - The number of blocks in the list
!!  blockList(:) - The list of blocks on which the solution must be updated
!!   dt           : The time step
!!   pass         : pass=1 directional order of solution sweep X-Y-Z, 
!!                  pass=2 directional order of solution sweep Z-Y-X.
!!  dt           - The time step
!!
!!
!!
!! NOTES
!!
!!  The interface of this subroutine must be explicitly known to code that
!!  calls it.  The simplest way to make it so is to have something like
!!    use diff_interface,ONLY: diff_advanceTherm
!!  in the calling routine.
!!***

subroutine diff_advanceTherm(blockCount,blockList,dt,pass)

  implicit none

  integer,intent(IN)                       :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  real,intent(in)                          :: dt
  integer, OPTIONAL, intent(IN)            :: pass
  
  return
end subroutine diff_advanceTherm
