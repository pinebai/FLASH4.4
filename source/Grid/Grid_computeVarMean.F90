!!****f* source/Grid/Grid_computeVarMean
!!
!! NAME
!!
!!  Grid_computeVarMean
!!
!! SYNOPSIS
!!
!!  call Grid_computeVarMean(integer(in) :: iunk,
!!                           real(out) :: mean)
!!
!! DESCRIPTION
!!
!!  Calculates the mean of a variable in UNK
!!
!! ARGUMENTS
!!
!!   iunk : the variable (index into the UNK array)
!!
!!   mean : the mean value returned
!!
!!
!!
!!***

subroutine Grid_computeVarMean(iUnk, mean)
  implicit none
  
  integer, intent(in) :: iUnk
  real, intent(out) :: mean

  mean = 0.0

end subroutine Grid_computeVarMean
