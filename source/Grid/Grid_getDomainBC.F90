!!****f* source/Grid/Grid_getDomainBC
!!
!! NAME
!!  Grid_getDomainBC
!!
!! SYNOPSIS
!!
!!  Grid_getDomainBC(integer(OUT) :: boundary(LOW:HIGH,MDIM),
!!                    
!! DESCRIPTION 
!!  Returns the boundary condition for each face of the domain
!!  
!!  The boundary conditions are defined in the header file
!!  constants.h, ie. OUTFLOW, REFLECTING, PERIODIC etc
!!   
!! ARGUMENTS 
!!
!!  boundary   - array returned holding boundary conditions, except when periodic,
!!            if any of the faces of the block are on a physical boundary. 
!!                
!!            the first index of the array can take on values LOW or
!!            HIGH, and the second index can be IAXIS, JAXIS or KAXIS
!!
!! NOTES
!!
!!   The #define constants LOW, HIGH, IAXIS, JAXIS and KAXIS
!!   are defined in constants.h and are
!!   meant to ease the readability of the code.       
!!   instead of boundary(2,3) = PERIODIC the code reads
!!   boundary(HIGH, KAXIS) = PERIODIC
!!
!!***


subroutine Grid_getDomainBC(boundary)

#include "constants.h"
  implicit none
  integer, dimension(LOW:HIGH,MDIM),intent(out):: boundary

  boundary(LOW:HIGH,1:MDIM) = NOT_BOUNDARY
 
  return
end subroutine Grid_getDomainBC





