!!****f* source/Grid/Grid_pfftMapFromOutput
!!
!! NAME 
!!
!!   Grid_pfftMapFromOutput
!!
!! SYNOPSIS
!!
!!   Grid_pfftMapFromOutput(integer(IN) :: gridVar,
!!                 real, dimension(:),(IN)   :: pfft_outArray)
!!
!! DESCRIPTION 
!!
!!  Takes the data from the Pfft grid and redistribute it back
!!  to the the variable gridVar in the mesh data structure.
!!  based upon thh map determined by the routine Grid_pfftInit.
!!
!! ARGUMENTS
!!
!!  gridVar          - variable on the mesh on which pfft is to be applies
!!  pfft_outArray     - array containing the output of Pfft
!!
!! NOTE 
!!
!!  Users must make sure that Grid_pfftInput has been called correctly
!!  before calling this routine.
!!
!!***

subroutine Grid_pfftMapFromOutput(gridVar,pfft_outArray)


  implicit none

  integer,intent(IN) :: gridVar
  real, dimension(:),intent(out) :: pfft_outArray  
  pfft_outArray = 0.0

  return
end subroutine Grid_pfftMapFromOutput
