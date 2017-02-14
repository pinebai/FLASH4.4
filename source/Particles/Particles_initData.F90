!!****f* source/Particles/Particles_initData
!!
!! NAME
!!
!!  Particles_initData
!!
!! SYNOPSIS
!!
!!  call Particles_initData(logical(IN)    :: restart,
!!                          logical(INOUT) :: partPosInitialized)
!!
!! DESCRIPTION
!!
!!  Initialization of particle attributes after they have been
!!  dropped into their positions.
!!
!!  This interface is called during Flash initialization after
!!   o  particles positions,
!!   o  mass properties for active particles,  (NOTE: it is true in
!!                                              most cases, but not a 
!!                                              requirement)
!!   o  and particle tags
!!  have been initialized (or possibly read from a checkpoint, if
!!  restart is true).
!!
!!  This is where velocity components are initialized.
!!  For passive particles, this is done by mapping from
!!  the fluid velocity field. The sequence of initialization
!!  is such that fluid fields are guaranteed to have been
!!  initialized to their final starting values, and the
!!  configuration of the grid has been finalized including
!!  creation of blocks up to the maximum refinement level,
!!  when the Particles_initData interface is called.
!!
!!  The default implementation in the ParticlesInitialization
!!  subunit calls Simulation_initParticleAttrib before returning,
!!  so if custom initialization of additional particle properties
!!  is needed, that's the place to implement it.
!!
!! ARGUMENTS
!!
!!  restart - true if restarting from a checkpoint, false otherwise.
!!
!! SIDE EFFECTS
!!
!!  Updates particles data that is private to the Particles unit.
!!
!!  The defaul implementation normally sets an internal flag
!!  (pt_velInitialized) that keeps track of whether initialization of
!!  particle velocities is complete.
!!
!! SEE ALSO
!!
!!  Driver_intFlash
!!  Simulation_initParticleAttrib
!!
!!***

subroutine Particles_initData(restart,partPosInitialized)
  implicit none
  logical, intent(IN) :: restart
  logical,intent(INOUT) :: partPosInitialized
  
  partPosInitialized=.false.

end subroutine Particles_initData
