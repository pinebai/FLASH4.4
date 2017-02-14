!!****if* source/RuntimeParameters/RuntimeParametersMain/rp_getOpt
!!
!! NAME
!!
!!  rp_getOpt
!!
!!
!! SYNOPSIS
!!
!!  call rp_getOpt(character(len=20)(in) :: optString, 
!!                 integer(in)           :: hasValue, 
!!                 integer(out)          :: exists, 
!!                 character(len=MAX_STRING_LENGTH)(out) :: optValue)
!!
!!
!! DESCRIPTION
!!
!!  Look for character string command line opt and optionally its
!!  associated value.  If opt is there, exists is set to 1.  If
!!  has_value is set to 1, then also try to parse out a value
!!  and return the value in the character variable opt_value.
!!
!! ARGUMENTS
!!
!!   optString :   input: opt character string to look for
!!   hasValue  :   input: if 1, look for value
!!   exists     :   output: set to 1 if opt found
!!   optValue  :   output: the associated value if the opt has one
!!
!! EXAMPLE
!!
!!  To get the name of a par  file from a command line like this:
!!   mpirun -np 1 flash2 -parfile flash.par
!!
!!   call rp_getOpt('-parfile', 1, exists, parameter_file_string)
!!
!!   if (exists == 1) then
!!     ... do stuff with parameter_file_string
!!
!!
!! NOTES
!!
!!  This does not support -- style arguments.
!!  
!!
!!***

subroutine rp_getOpt(optString, hasValue, exists, optValue)

#ifdef NAGF95
  use f90_unix_env, ONLY: iargc
#endif

  implicit none

#include "constants.h"

  character(len=20), INTENT(IN) :: optString
  integer, INTENT(IN) :: hasValue
  integer, INTENT(OUT) :: exists
  character (len=MAX_STRING_LENGTH), INTENT(OUT) :: optValue

  character (len=MAX_STRING_LENGTH) :: currentArg
  integer  nargs, i, argLength,strlen
  logical found

#ifndef NAGF95
  integer iargc
#endif

  nargs = iargc()
  strlen = len_trim(optString)
  found = .false.
  i = 1

  exists = 0

  do while (.not. found .and. i <= nargs)
     call rp_getArgument(i, currentArg)
     argLength = len_trim(currentArg)
#ifdef DEBUG_RP
     print*,optString             ! these are useful for DEBUG
     print*,strlen,argLength      ! DEBUG
#endif
     if ((argLength==(strlen)).and.(currentArg==optString)) then
        exists = 1
        found = .true.
        if (hasValue == 1) then
           call rp_getArgument(i+1, optValue)
        end if
     end if
     i = i + 1
  end do


end subroutine rp_getOpt








