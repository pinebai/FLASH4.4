!!****if* source/flashUtilities/general/removeNullChar
!!
!! NAME
!!
!!  removeNullChar
!!
!!
!! SYNOPSIS
!!
!!  removeNullChar(intent(inout) :: str1) 
!!
!!
!! DESCRIPTION
!!
!!  This is a small little routine that takes a string and
!!  looks for null chars at the end.  If it finds one it replaces
!!  it with a blank character.  Additionally, if one NullChar is 
!!  found, any following chars are set to a blank char.
!!   This routine is necessary because
!!  fortran and C handle strings differently.  C attaches the null
!!  character to the end of strings.  Since most of our work is
!!  done in fortran this rarely matters, but in the case of I/O
!!  particularly reading in data from a checkpoint file.  C strings
!!  can be different from Fortran strings
!!  
!!  
!!  
!!
!! ARGUMENTS
!! 
!!  str1 - string with null chars returned with null chars removed
!!  
!!
!!***

subroutine removeNullChar(str1)

  implicit none
 
  character(len=*), intent(inout) :: str1
  integer :: i
  logical :: nullCharFound
  integer :: temp_int


  nullCharFound = .false.

  do i=1, len(str1)
     if(nullCharFound) then
        str1(i:i) = " "
     else if(str1(i:i) == "\0") then
        str1(i:i) = " "
        nullCharFound = .true.
     else 
        temp_int = ichar(str1(i:i))
        if(temp_int == 0) then
           str1(i:i) = " "
           nullCharFound = .true.
        end if
     end if
     
  end do


end subroutine removeNullChar
