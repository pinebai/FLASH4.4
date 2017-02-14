!!****if* source/Driver/DriverMain/Driver_mpiThreadSupport
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
  use Driver_data, ONLY : dr_mpiThreadSupport
  implicit none
  logical, intent(OUT) :: mpiThreadSupport
  mpiThreadSupport = dr_mpiThreadSupport
end subroutine Driver_mpiThreadSupport
