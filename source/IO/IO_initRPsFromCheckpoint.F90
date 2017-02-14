!!****f* source/IO/IO_initRPsFromCheckpoint
!!
!! NAME
!!
!!  IO_initRPsFromCheckpoint
!!
!! SYNOPSIS
!!
!!  call IO_initRPsFromCheckpoint(character(len=*)(IN) :: filename,
!!                                integer(OUT)         :: ierr)
!!
!! DESCRIPTION
!!
!!   Reads runtime parameter information and some integer scalars
!!   from a checkpoint file, then sets the runtime parameters in
!!   the linked list for the current run accordingly.
!!
!!   
!! ARGUMENTS
!!
!!
!!   filename : name of checkpoint file
!!
!!   ierr : if not 0, indicates that something went wrong
!!
!! NOTES
!!
!!   This subroutine is called early during initialization.
!!   The implementation, and subroutines called from it, must
!!   operate correctly before IO_init has been called.
!!
!!   This is not meant to be called by users.
!!   
!!***

subroutine IO_initRPsFromCheckpoint( filename, ierr)
  implicit none
  character(len=*),intent(IN) :: filename
  integer,intent(OUT) :: ierr

  ierr = 0

end subroutine IO_initRPsFromCheckpoint
