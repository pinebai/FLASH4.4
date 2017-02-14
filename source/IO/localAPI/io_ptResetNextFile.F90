!!****if* source/IO/localAPI/io_ptResetNextFile
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
!!  savednext : the externally stored next time to output a particle fil
!!
!!
!!
!!***

subroutine io_ptResetNextFile(savedNext)


  implicit none
  real, intent(IN) :: savedNext


end subroutine
