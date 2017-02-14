!!****if* source/RuntimeParameters/RuntimeParametersMain/rp_getArgument
!!
!! NAME 
!!    rp_getArgument
!!
!! SYNOPSIS
!!    call rp_getArgument(integer(in) :: pos,
!!                        character(len=MAX_STRING_LENGTH)(out) :: arg)
!!
!! DESCRIPTION
!!
!!     Plug-in routine to return a command-line argument.
!!     Specify the position of the argument on the command
!!     line (beginning with 1) and a string variable to
!!     receive the argument value.
!!
!!     This is a wrapper which can be replaced when running
!!     on systems that do not have this extension to the standard.
!!
!! ARGUMENTS
!!
!!     pos :   integer argument position
!!     arg :  the returned command-line argument
!!
!!
!!
!!***
 
 
subroutine rp_getArgument (pos, arg)

#ifdef NAGF95
  use f90_unix_env, ONLY: getarg
#endif

  implicit none
 
#include "constants.h"

  integer, intent(in) ::          pos
  character(len=MAX_STRING_LENGTH), intent(out) :: arg
 
  call getarg (pos, arg)
 
  return
end subroutine rp_getArgument
