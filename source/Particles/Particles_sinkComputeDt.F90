!!****f* source/Particles/Particles_sinkComputeDt
!!
!! NAME
!!
!!  Particles_sinkComputeDt
!!
!! SYNOPSIS
!!
!!  call Particles_sinkComputeDt(integer, INTENT(in) :: blockid,
!!                               real, INTENT(inout) :: dt_sink,
!!                            integer, INTENT(inout) :: dt_minloc)
!!
!! DESCRIPTION
!!
!!  Constrains the global timestep based on sink particle velocities and
!!  gravitational accelerations.
!!
!! ARGUMENTS
!!
!!   blockid - ID of block in current processor
!!
!!   dt_sink - dt constrained by sinks
!!
!!   dt_minloc - location of the cell that constrains the sink time step
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2008-2015
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   added acceleration timestep constraint; Christoph Federrath, 2013
!!   added sinks off-domain support; Christoph Federrath, 2015
!!
!!***

subroutine Particles_sinkComputeDt(blockID,dt_sink,dt_minloc)
  implicit none
  integer, INTENT(in)    :: blockID
  real, INTENT(inout)    :: dt_sink
  integer, INTENT(inout) :: dt_minloc(5)
end subroutine Particles_sinkComputeDt
