!!****f* source/Simulation/Simulation_getRenormGroup
!!
!! NAME
!!
!!  Simulation_getRenormGroup
!!
!! SYNOPSIS
!!
!!  Simulation_getRenormGroup(integer, intent(in)  :: mscalar,
!!                            integer, intent(out)  :: group)
!!
!! DESCRIPTION
!!  Takes the index of a Mass Scalar and returns the renormalization group
!!  in which it appears.
!!
!!  For Mass Scalars which do not need to be normalized, the group is 0.
!!
!! ARGUMENTS
!!
!!   mscalar : The target mass scalar
!!
!!   group : The renormalization group to which mscalar belongs.
!!           0 if mscalar does not need to be renormalized.
!!
!!
!! NOTES
!!  An implementation of this subroutine is usually automatically
!!  generated at setup time and should not be edited.
!!
!!***


!!  call Simulation_getRenormGroup(mscalar,group)
!!
!!  Takes the index of a Mass Scalar and returns the renormalization group
!!  in which it appears. 
!!   
!!  For MassScalars which don't need to be normalized, the group is 0
!!


subroutine Simulation_getRenormGroup(mscalar,group)
implicit none 

#include "constants.h"
#include "Flash.h"

   integer, intent(out) ::group 
   integer, intent(in) :: mscalar

   group = 0 !! Dont renorm group 0, norm everything else

end subroutine Simulation_getRenormGroup

