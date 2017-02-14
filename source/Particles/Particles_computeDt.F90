!!****f* source/Particles/Particles_computeDt
!!
!! NAME
!!  Particles_computeDt
!!
!! SYNOPSIS
!!  Particles_computeDt( integer(in) :: blockID,
!!                       real(inout) :: dt_part,
!!                       integer(5)(inout) :: dt_minloc)
!!
!! DESCRIPTION
!!   Timestep computation routine for the particle unit.
!!   This routine sets the timestep by requiring that
!!   particles travel no more than some fraction pt_dtFactor
!!   during a single step.  Even for time integration
!!   schemes that can deal with pt_dtFactor > 1, you should
!!   not set pt_dtFactor to be larger than half the number
!!   of guard cells.  Otherwise particles that leave the block
!!   that "owns" them may overshoot the block's immediate
!!   neighbors, causing problems when particle data are transmitted
!!   to the neighbors.
!!
!!
!! ARGUMENTS
!!
!!    blockID:         local block ID
!!    
!!    dt_part:         variable to hold timestep constraint
!!    dt_minloc(5):    array to hold limiting zone info:  index[1-3] are zone
!!                     indices (i,j,k); index[4]=block ID; index[5]=processor 
!!                     number. The zone indices indicate the zone containing
!!                     the particle whose velocity restricted the timestep.
!!
!! PARAMETERS
!!
!!   pt_dtFactor:      REAL [default 0.5] A factor multiplying dx/|v| to limit
!!                       the movement of particles outside a single block. 
!!
!!
!! NOTES
!!
!!   Note that this routine does NOT calculate a timestep limitation
!!   based on the theoretical stability limits of the time integration
!!   routines.  It SHOULD, but these limits were deemed to be unlikely
!!   to be restrictive based on the high velocities involved in normal
!!   astronomical simulations.
!!
!!***

 subroutine Particles_computeDt (blockID, dt_part, dt_minloc)

!===============================================================================

  implicit none
    
  
  real, INTENT(inout)    :: dt_part
  integer, INTENT(inout) :: dt_minloc(5)
  integer, INTENT(in)    :: blockID

  
  return
end subroutine Particles_computeDt

