!!****if* source/flashUtilities/general/concatStringWithInt
!!
!! NAME
!!    concatStringWithInt
!!
!! SYNOPSIS
!!
!!    concatStringWithInt()
!!
!! DESCRIPTION
!!
!!  This routine is used for generating a new string name where
!!  the base name is appended with some integer value at the end
!!  For example if stringBaseName is "refine_var_", and the value
!!  of the integer is "1", then the routine will return "refine_var_1"
!!
!!
!!***

subroutine concatStringWithInt(inString, ind, outString )
  
#include "constants.h"

  implicit none
  character(len=*),intent(IN) :: inString
  integer,intent(IN) :: ind
  character(len=MAX_STRING_LENGTH),intent(OUT) :: outString
  
  
  character(len=MAX_STRING_LENGTH) :: parname
  integer :: pos


  pos=len_trim(inString)+1
  parname= ""
  write(parname, *) ind
  parname = ADJUSTL(parname)
  outString=inString(:pos-1)//parname
  

end subroutine concatStringWithInt
