!!****if* source/IO/IOParticles/io_ptCorrectNextPartTime
!!
!! NAME
!!
!!  io_ptCorrectNextPartTime
!!
!! SYNOPSIS
!!
!!  io_ptCorrectNextPartTime(real, intent(IN)  :: simtime,
!!                           real, intent(OUT)  :: oldtime)
!!
!! DESCRIPTION
!!
!!  This is used in correcting the particle file output time during 
!!  checkpointing.  Without correcting this time in the checkpoint file,
!!  on restart, the count for the next particle file output will be
!!  reset to zero rather than being in the same place that the simulation was
!!  when the checkpoint file was output.
!!
!! ARGUMENTS
!!
!!   simtime : the current simulation time 
!!
!!   oldtime : the next time to output a file by
!!
!!
!!
!!***

subroutine io_ptCorrectNextPartTime(simTime, oldTime)

  use IOParticles_data, ONLY : io_nextParticleFileTime, io_particleFileIntervalTime

  implicit none

  real, intent(IN) :: simTime
  real, intent(OUT) :: oldTime
  real :: newTime

  oldTime = io_nextParticleFileTime

  if( simTime >= io_nextParticleFileTime .and. &
     io_particleFileIntervalTime > 0.e0 ) then
    io_nextParticleFileTime = io_nextParticleFileTime + &
      (int((simTime-io_nextParticleFileTime)/io_particleFileIntervalTime)+1) &
       * io_particleFileIntervalTime
  end if

  

end subroutine io_ptCorrectNextPartTime
