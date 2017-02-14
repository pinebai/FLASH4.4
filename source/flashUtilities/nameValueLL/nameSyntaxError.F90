!!****if* source/flashUtilities/nameValueLL/nameSyntaxError
!!
!! NAME
!!  nameSyntaxError
!!
!! SYNOPSIS
!!
!!  nameSyntaxError(character(len=*)(in) :: inline)
!!
!! DESCRIPTION
!!
!!  Signals a syntax error in the runtime parameter file.
!!  Currently just echoes the offending line
!!
!!
!! ARGUMENTS
!!
!!  inline - string with syntax error that gets printed to stdout
!!
!!***


subroutine nameSyntaxError (inline)
  
  implicit none
  
  character(len=*),intent(in) :: inline
  
  write (*,*) 'read :  skipping the following', & 
       &              ' line due to a syntax error:'
  write (*,*) inline
  
  return
end subroutine nameSyntaxError
