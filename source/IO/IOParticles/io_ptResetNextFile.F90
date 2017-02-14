!!****if* source/IO/IOParticles/io_ptResetNextFile
!!
!! NAME
!!
!!  io_ptResetNextFile
!!
!! SYNOPSIS
!!
!!  io_ptResetNextFile(real, intent(IN)  :: savednext)
!!
!! DESCRIPTION
!!
!!  This allows you to reset the next time to write a particle file to the 
!!  value sacedNext.  This is necessary to preserve the next output time
!!  across restarts.
!!
!! ARGUMENTS
!!
!!   savednext : the externally stored next time to output a particle file
!!
!!
!!
!!***

subroutine io_ptResetNextFile(savedNext)

  use IOParticles_data, ONLY : io_nextParticleFileTime
 
  implicit none
  
  real, intent(IN) :: savedNext
  
  io_nextParticleFileTime = savedNext

end subroutine
