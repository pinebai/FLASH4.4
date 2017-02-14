!!****if* source/Grid/GridMain/Grid_computeVarMean
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
  use gr_interface, ONLY: gr_findMean
  implicit none
  
  integer, intent(in) :: iUnk
  real, intent(out) :: mean

  call gr_findMean(iUnk, iType=2, bGuardCell=.FALSE., mean=mean)

end subroutine Grid_computeVarMean
