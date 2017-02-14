!!****f* source/Particles/Particles_updateRefinement
!!
!! NAME
!!
!!  Particles_updateRefinement
!!
!! SYNOPSIS
!!
!!  Particles_updateRefinement(real(inout) :: oldLocalNumBlocks)
!!
!! DESCRIPTION
!!   This routine provides a hook into the particle data structure
!!   for the Grid. It is called during update Refinement process
!!   by the Grid. The routine passes the control right back
!!   grid, with Particles specific data structures in the argument
!!   list, so that Grid can operate on them.
!!
!! ARGUMENTS
!!
!!    oldLocalNumBlocks :   number of blocks on a processor before 
!!                          refinement. 
!!
!! PARAMETERS
!!  
!!
!!***

subroutine Particles_updateRefinement(oldLocalNumBlocks)
  implicit none 
  integer,intent(INOUT) :: oldLocalNumBlocks

  return
end subroutine Particles_updateRefinement
