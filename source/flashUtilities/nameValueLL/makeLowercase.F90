!!****if* source/flashUtilities/nameValueLL/makeLowercase
!! 
!!  NAME     
!!    makeLowercase
!! 
!!  SYNOPSIS
!!    makeLowercase(character(len=*)(INOUT) :: str)
!!
!!  DESCRIPTION 
!!    Convert a string to lowercase.
!!
!!  ARGUMENTS
!!    str:    string to be converted 
!!
!!***

subroutine makeLowercase (str)
  
implicit none
  character(len=*),intent(INOUT) :: str
  integer         :: i
  
  do i = 1, len_trim(str)
     if (lge(str(i:i), 'A') .and. lle(str(i:i), 'Z')) & 
          &      str(i:i) = achar( iachar(str(i:i)) + 32 )
  enddo
  
  return
end subroutine makeLowercase
