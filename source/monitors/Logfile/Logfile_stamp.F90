!!****f* source/monitors/Logfile/Logfile_stamp
!! NAME
!!   
!! Logfile_stamp
!!
!! SYNOPSIS
!!  Logfile_stamp(int/real/str/log(in) :: val,
!!                char*(in)            :: tag,
!!                char*(in)            :: attrib)
!!
!! DESCRIPTION
!!
!!   Logfile_stamp is an overloaded subroutine and includes
!!   Logfile_stampInt
!!   Logfile_stampReal
!!   Logfile_stampStr
!!   Logfile_stampIntArray
!!   Logfile_stampRealArray
!!   Logfile_stampStrArray
!!   Logfile_stampStrPair
!!   
!!   Each subroutine stamps the date and time along with a value and
!!   message into the logfile.  Typical uses are when the grid package 
!!   stamps the logfile
!!   if it refines or derefines the grid.  The IO unit will call 
!!   Logfile_stamp when it reads or writes data to a checkpoint file.
!!
!! ARGUMENTS
!!
!!  Arguments vary slightly for the routines but in general
!!  val  - int_val, real_val, intArrayVal, a value to put into the logfile
!!  tag  - string identifier, often of the routine stamping the logfile.
!!         (see example below)
!!  attrib - !!DEV! optional argument, not currently implemented.
!!
!! NOTES
!!  
!!  Because Logfile_stamp is an overloaded subroutine _and_ also has
!!  optional arguments, most compilers require that any routine calling
!!  Logfile_stamp must include the header file Logfile.h.
!!
!!
!!  Variables that begin with "log_" are defined in the fortran 
!!  module Logfile_data.  The prefix "log_" is meant to indicate
!!  that these variables have Logfile unit scope.  Other variables
!!  are local to the individual subroutines
!!
!! EXAMPLE
!!   example call for stamping a string
!!   call Logfile_stamp( 'read 159 blocks from file', '[IO_readCheckpoint]')
!!
!! 
!!***




#include "constants.h"

subroutine Logfile_stampInt( intVal, tag, attrib)

implicit none
  integer, intent(in)                   ::  intVal  
  character(len=*), intent(in)           :: tag
  character(len=*), intent(in), OPTIONAL :: attrib

  return
  
end subroutine Logfile_stampInt

subroutine Logfile_stampIntArray( intArr, len, tag, attrib)

implicit none

  character(len=*),intent(in)           :: tag
  character(len=*), intent(in),OPTIONAL :: attrib
  integer,intent(in)                    :: len
  integer, dimension(len), intent(in)   :: intArr

  return
  
end subroutine Logfile_stampIntArray

subroutine Logfile_stampReal( realVal, tag, attrib)

  implicit none
  real, intent(in)                       :: realVal
  character(len=*), intent(in)           :: tag
  character(len=*), intent(in), OPTIONAL  :: attrib
  
  return
  
end subroutine Logfile_stampReal

subroutine Logfile_stampRealArray( realArr, len, tag, attrib)

  implicit none
  character(len=*),intent(in)            :: tag
  character(len=*), intent(in), OPTIONAL :: attrib
  integer,intent(in)                     :: len
  real, dimension(len),intent(in)        :: realArr

  return
  
end subroutine Logfile_stampRealArray

subroutine Logfile_stampStr( string, tag, attrib)

implicit none
  character(len=*), intent(in)           :: string
  character(len=*), intent(in), OPTIONAL :: tag, attrib
  
  return
  
end subroutine Logfile_stampStr

subroutine Logfile_stampStrArray( strArr, len, tag, attrib)

implicit none
  character(len=*), intent(in)             :: tag
  character(len=*), intent(in), OPTIONAL   :: attrib
  integer,intent(in)                       :: len
  character(len=*), dimension(len),intent(in) :: strArr
  
  return    
  
end subroutine Logfile_stampStrArray


subroutine Logfile_stampStrPair( strArr, len, dimen, tag, attrib)

implicit none

  character(len=*),intent(in)           :: tag
  character(len=*),intent(in), OPTIONAL :: attrib
  integer,intent(in)                    :: len, dimen
  character(len=*),dimension(len,dimen),intent(in) :: strArr
  
  return    
  
end subroutine Logfile_stampStrPair
