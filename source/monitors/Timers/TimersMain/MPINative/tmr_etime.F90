!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_etime
!!
!! NAME
!!   
!! tmr_etime
!!
!! SYNOPSIS
!!
!!  call tmr_etime(real(OUT) :: time)
!!
!! DESCRIPTION
!!
!!  Return the elapsed time.
!!
!! ARGUMENTS
!!
!!  time - elapsed time
!!
!!***
subroutine tmr_etime (time)
 
  implicit none 
  include 'mpif.h'
  real, intent(out) :: time

  time = MPI_WTime()
  return
end subroutine tmr_etime
