!!****f* source/Driver/Driver_mpiThreadSupport
!!
!! NAME
!!  Driver_mpiThreadSupport
!!
!! SYNOPSIS
!!
!!  Driver_mpiThreadSupport(logical(OUT) :: mpiThreadSupport)
!!               
!!  
!! DESCRIPTION 
!!
!!  Confirms whether the MPI library has provided us with a safe
!!  thread support level.
!!
!! ARGUMENTS
!!
!!  mpiThreadSupport - Safe level? (TRUE/FALSE)
!!
!!***

subroutine Driver_mpiThreadSupport(mpiThreadSupport)
  implicit none
  logical, intent(OUT) :: mpiThreadSupport
  mpiThreadSupport = .false.
end subroutine Driver_mpiThreadSupport
