!!****f* source/Simulation/Simulation_initParticleAttrib
!!
!! NAME
!!
!!  Simulation_initParticleAttrib
!!
!! SYNOPSIS
!!
!!  call Simulation_initParticleAttrib(logical(IN) :: restart)
!!
!! DESCRIPTION
!!
!!  A stub for initializing additional particle attributes if necessary.
!!
!!  This interface is called during initialization after
!!   o  particles positions,
!!   o  particle velocities and other properties needed to advance
!!      particle trajectories,
!!   o  mass properties for active particles,
!!   o  and particle tags
!!  have all been initialized (or possibly read from a checkpoint, if
!!  restart is true).
!!  
!!
!!  By default this does nothing. A typical use would be to initialize
!!  particle properties that were defined in a simulation directory.
!!  It is likely that an implementation only needs to take action
!!  if the restart flag is false.
!!
!!
!! ARGUMENTS
!!
!!  restart - true if restarting from a checkpoint, false otherwise.
!!
!! SEE ALSO
!!
!!  Driver_intFlash
!!  Particles_initPositions
!!
!!***

subroutine Simulation_initParticleAttrib(restart)
  implicit none
  logical,intent(in) :: restart

end subroutine Simulation_initParticleAttrib
