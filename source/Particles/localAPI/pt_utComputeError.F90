!!****if* source/Particles/localAPI/pt_utComputeError
!!
!! NAME
!!
!!  pt_utComputeError
!!
!! SYNOPSIS
!!
!!  pt_utComputeError(real(in) :: dtOld,
!!                           real(in) :: dtNew,
!!                           real(in) :: t)
!!
!! DESCRIPTION
!!
!!  Compute error, that is, compare actual to analytic solution.
!!
!!  Compares particles' POS{X,Y,Z}_PART_PROP and POSANAL{X,Y,Z}_PART_PROP
!!  properties.
!!
!! ARGUMENTS
!!
!!   dtOld -- previous time increment
!!   dtNew -- current time increment
!!   t     -- time at which solutions are compared
!!  
!!***

!===============================================================================

subroutine pt_utComputeError (dtOld,dtNew,t)
    
  
  implicit none

  
  real, INTENT(in)  :: dtOld, dtNew, t
     
  return
!!------------------------------------------------------------------------------
  
end subroutine pt_utComputeError


