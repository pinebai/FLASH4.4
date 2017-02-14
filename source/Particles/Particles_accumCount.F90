!!****f* source/Particles/Particles_accumCount
!!
!! NAME
!!    Particles_accumCount
!!
!! SYNOPSIS
!!
!!    Particles_accumCount(integer, intent(IN) :: var)
!!
!! DESCRIPTION
!!    
!!   This routine is to be used in the refinement based upon the number of
!!   particles in a block. It adds a weight to the cell of the grid in the
!!   specified grid variable if a particle is found to be contained in the 
!!   cell 
!!
!! ARGUMENTS
!!
!!  var:        the grid variable to add the weights if particle found in the cell
!!
!!
!!***


subroutine Particles_accumCount(var)

  implicit none

  integer, intent(IN) :: var
!----------------------------------------------------------------------

  
  return
  
  !----------------------------------------------------------------------
  
end subroutine Particles_accumCount


