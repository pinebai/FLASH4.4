!!****f* source/RuntimeParameters/RuntimeParameters_read
!!
!! NAME
!!  RuntimeParameters_read
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_read( character(in), (len=MAX_STRING_LEN) :: parmfile)
!!
!!
!! DESCRIPTION
!!
!!
!!         Parses the parameter file parmfile, which contains
!!         job-dependent parameter definitions.  Syntax of
!!         parameter file lines is:
!!
!!                       # comment               rest of line is a comment
!!                       variable = value        set variable to value
!!                       strvar   = "word"     set string variable strvar
!!                                                to word
!!
!!         Some error checking is performed:  if a variable is
!!         unrecognized, it is ignored.  Syntax errors are signaled
!!         along with the offending lines.  However, type mismatch
!!         errors force termination.  Also, no checking is done to
!!         determine whether all of the variables have received some
!!         value.  In case of repeated definitions of the same
!!         variable, the last definition overrides the others.
!!
!! ARGUMENTS
!!
!!   parmfile :       the name of the parameter file to read
!! 
!! NOTES
!!
!!   This routine is called during FLASH initialization to read the flash.par file.
!!   In general it would not be used by users.  Instead, use the restart = .true.
!!   capability within a flash.par to restart a run from a checkpoint file,
!!   using a different flash.par configuration.
!!
!!
!!
!!***



subroutine RuntimeParameters_read (parmfile)

#include "constants.h"

  implicit none
  

  character(len=MAX_STRING_LENGTH), intent(in)    :: parmfile


  return
end subroutine RuntimeParameters_read



