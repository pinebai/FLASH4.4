!!****f* source/Particles/Particles_addNew
!!
!! NAME
!!    Particles_addNew
!!
!! SYNOPSIS
!!    call Particles_addNew( integer(in)  :: count,
!!                  optional,real(in)     :: pos(MDIM,count),
!!                           logical(out) :: success)
!!
!! DESCRIPTION
!!
!!    This routine allows particles to be added during evolution.
!!    In the particles data structure it always initializes the tag and
!!    processor ID. If the optional argument "pos" is present then it
!!    will also initialize the position and block ID attributes in the
!!    particles. It returns the value FALSE in success if there isn't
!!    enough space left in the data structure for the requested number
!!    of particles.
!!
!! ARGUMENTS
!!
!!     count   :: the count of particles to be added
!!     pos     :: optional, contains the coordinates of the particles
!!     success :: This arg returns TRUE if there was enough space 
!!                in the particles data structure to add the requested
!!                number of particles, FALSE otherwise.
!!
!!    
!!  NOTES
!!
!! The constant MDIM is defined in constants.h .
!!
!!***

!!#define DEBUG_PARTICLES

subroutine Particles_addNew (count, pos, success)
  
  implicit none
  
#include "constants.h"

  integer, INTENT(in) :: count
  real, optional, dimension(MDIM,count), intent(IN)::pos
  logical, intent(OUT) :: success

  success = .FALSE.
  return
  
end subroutine Particles_addNew
