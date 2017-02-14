!!****f* source/Grid/Grid_initDomain
!!
!! NAME
!!  Grid_initDomain
!!
!! SYNOPSIS
!!
!!  call Grid_initDomain(logical(IN)  :: restart,
!!                       logical(INOUT) :: particlesInitialized)
!!
!!
!! DESCRIPTION
!!  Initialize the Grid's data stuctures for the discretized mesh and
!!  apply the initial conditions.
!!
!! ARGUMENTS
!!  restart  : whether starting from initial conditions or from checkpoint
!!  particlesInitialized : is true if particle positions were initialized before returning
!!                         from this routine
!!
!!***

subroutine Grid_initDomain( restart,particlesInitialized)

implicit none
  logical, intent(IN) :: restart
  logical, intent(INOUT) :: particlesInitialized
  
end subroutine Grid_initDomain
