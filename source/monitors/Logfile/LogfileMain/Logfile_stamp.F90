!!****if* source/monitors/Logfile/LogfileMain/Logfile_stamp
!! NAME
!!   
!! Logfile_stamp
!!
!! SYNOPSIS
!!
!!  Logfile_stamp(int/real/str/log(in) :: val,
!!                char*(in) :: tag,
!!                char*(in) :: attrib)
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
!!   stamps the logfile if it refines or derefines the grid.  The IO
!!   unit will call Logfile_stamp when it reads or writes data to a
!!   checkpoint file.
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
!!  Variables that begin with "log_" are defined in the fortran module
!!  Logfile_data.  The prefix "log_" is meant to indicate that these
!!  variables have Logfile unit scope.  Other variables are local to
!!  the individual subroutines
!!
!! EXAMPLE
!!   example call for stamping a string
!!   Logfile_stamp( 'read 159 blocks from file', '[IO_readCheckpoint]')
!!
!! 
!!***

#include "Flash.h"

#ifndef FLASH_FLUSH
! Do not actually flush if FLASH_FLUSH is not already defined
#define FLASH_FLUSH(u)
#endif

subroutine Logfile_stampInt( intVal, tag, attrib)

  use Logfile_data, ONLY : log_globalMe,  log_fileOpen, log_lun, &
       log_keepOpenAfterStamp, log_flushLevel
  use Logfile_interface, ONLY : Logfile_close, Logfile_open

  implicit none

#include "constants.h"

  integer, intent(in)                   ::  intVal  
  character(len=*), intent(in)          :: tag
  character(len=*), intent(in), OPTIONAL  :: attrib

  character(len=32)                     :: numToStr
  character(len=40), save               :: dateStr

  integer :: logUnit
  logical :: logUnitLocal=.false.

  if (log_globalMe .ne. MASTER_PE) return  ! only master processor writes to logfile

  if (.not. log_fileOpen) then
     call Logfile_open(logUnit,logUnitLocal)
  else
     logUnit = log_lun
  end if

  call current_date_time (dateStr)
  write(numToStr, "(I12)") intVal
  if (log_fileOpen) &
       &    write (log_lun, '(" [ ", A," ] ", A, ": ", A)') dateStr(1:len_trim(dateStr)), trim(tag), &
       &    trim(adjustl(numToStr))


  if (log_fileOpen .AND. .NOT. log_keepOpenAfterStamp) then
     call Logfile_close()
  else if (log_fileOpen .AND. log_flushLevel > 0) then
     FLASH_FLUSH(logUnit)
  end if

  return
  
end subroutine Logfile_stampInt

subroutine Logfile_stampIntArray( intArr, len, tag, attrib)

  use Logfile_data, ONLY : log_globalMe,  log_fileOpen, log_lun, &
       log_keepOpenAfterStamp, log_flushLevel
  use Logfile_interface, ONLY : Logfile_close, Logfile_open

  implicit none

#include "constants.h"

  character(len=*), intent(in)           :: tag
  character(len=*), OPTIONAL, intent(in) :: attrib
  integer, intent(in)                    :: len
  integer, dimension(len),intent(in)     :: intArr

  character(len=32)          :: numToStr
  character(len=40), save    :: dateStr    
  integer                    :: i
  integer :: logUnit
  logical :: logUnitLocal=.false.


  if (log_globalMe .ne. MASTER_PE) return  ! only master processor writes to logfile
  
  if (.not. log_fileOpen) then
     call Logfile_open(logUnit,logUnitLocal)
  else
     logUnit = log_lun
  end if
  
  call current_date_time (dateStr)
  
  do i=1, len
     write(numToStr, "(I12)") intArr(i)
     if (log_fileOpen) &
          &  write (log_lun, '(" [ ", A," ] ", A, ": ", A)') dateStr(1:len_trim(dateStr)), trim(tag), &
          &  trim(adjustl(numToStr))
  end do


  if (log_fileOpen .AND. .NOT. log_keepOpenAfterStamp) then
     call Logfile_close()
  else if (log_fileOpen .AND. log_flushLevel > 0) then
     FLASH_FLUSH(logUnit)
  end if
  
  return
  
end subroutine Logfile_stampIntArray

subroutine Logfile_stampReal( realVal, tag, attrib)

  use Logfile_data, ONLY : log_globalMe,  log_fileOpen, log_lun, &
       log_keepOpenAfterStamp, log_flushLevel
  use Logfile_interface, ONLY : Logfile_close, Logfile_open

  implicit none

#include "constants.h"
  real, intent(in)                       :: realVal
  character(len=*), intent(in)           :: tag
  character(len=*), OPTIONAL, intent(in) :: attrib
  character(len=32)          :: numToStr
  character(len=40), save    :: dateStr    
  integer :: logUnit
  logical :: logUnitLocal=.false.


  if (log_globalMe .ne. MASTER_PE) return  ! only master processor writes to logfile

  if (.not. log_fileOpen) then
     call Logfile_open(logUnit,logUnitLocal)
  else
     logUnit = log_lun
  end if

  call current_date_time (dateStr)

  write(numToStr, "(ES20.13)") realVal
  if (log_fileOpen) &
       &    write (log_lun, '(" [ ", A," ] ", A, ": ", A)') dateStr(1:len_trim(dateStr)), trim(tag), &
       &    trim(adjustl(numToStr))


  if (log_fileOpen .AND. .NOT. log_keepOpenAfterStamp) then
     call Logfile_close()
  else if (log_fileOpen .AND. log_flushLevel > 0) then
     FLASH_FLUSH(logUnit)
  end if
  
  return
  
end subroutine Logfile_stampReal

subroutine Logfile_stampRealArray( realArr, len, tag, attrib)

  use Logfile_data, ONLY : log_globalMe,  log_fileOpen, log_lun, &
       log_keepOpenAfterStamp, log_flushLevel
  use Logfile_interface, ONLY : Logfile_close, Logfile_open

  implicit none

#include "constants.h"
  character(len=*),intent(in)            :: tag
  character(len=*), OPTIONAL, intent(in) :: attrib
  integer, intent(in)                    :: len 
  real, dimension(len), intent(in)       :: realArr
  integer                    :: i
  integer :: logUnit
  logical :: logUnitLocal=.false.

  character(len=32)          :: numToStr
  character(len=40), save    :: dateStr    


  if (log_globalMe .ne. MASTER_PE) return  ! only master processor writes to logfile
  
  if (.not. log_fileOpen) then
     call Logfile_open(logUnit,logUnitLocal)
  else
     logUnit = log_lun
  end if
  
  call current_date_time (dateStr)

  do i=1, len
     write(numToStr, "(ES20.13)") realArr(i)
     if (log_fileOpen) &
          &  write (log_lun, '(" [ ", A," ] ", A, ": ", A)') dateStr(1:len_trim(dateStr)), trim(tag), &
          &  trim(adjustl(numToStr))
  end do



  if (log_fileOpen .AND. .NOT. log_keepOpenAfterStamp) then
     call Logfile_close()
  else if (log_fileOpen .AND. log_flushLevel > 0) then
     FLASH_FLUSH(logUnit)
  end if
  return
  
end subroutine Logfile_stampRealArray

subroutine Logfile_stampStr( string, tag, attrib)

  use Logfile_data, ONLY : log_globalMe,  log_fileOpen, log_lun, &
       log_keepOpenAfterStamp, log_flushLevel
  use Logfile_interface, ONLY : Logfile_close, Logfile_open

  implicit none

#include "constants.h"
  character(len=*), intent(in)           :: string
  character(len=*), OPTIONAL, intent(in) :: tag, attrib
  character(len=40), save                :: dateStr    
  character(len=MAX_STRING_LENGTH)       :: myTag
  integer :: logUnit
  logical :: logUnitLocal=.false.




  if (log_globalMe .ne. MASTER_PE) return  ! only master processor writes to logfile

  if (.not. log_fileOpen) then
     call Logfile_open(logUnit,logUnitLocal)
  else
     logUnit = log_lun
  end if
  
  call current_date_time (dateStr)

  if (.NOT. present(tag)) then
     myTag = "message"
  else
     myTag = tag
  end if


  if (log_fileOpen) &
       &     write (log_lun, '(" [ ", A, " ] ", A, ": ", A)') dateStr(1:len_trim(dateStr)), trim(myTag), trim(string)


  if (log_fileOpen .AND. .NOT. log_keepOpenAfterStamp) then
     call Logfile_close()
  else if (log_fileOpen .AND. log_flushLevel > 0) then
     FLASH_FLUSH(logUnit)
  end if

  return
  
end subroutine Logfile_stampStr













subroutine Logfile_stampStrArray( strArr, len, tag, attrib)

  use Logfile_data, ONLY : log_globalMe,  log_fileOpen, log_lun, &
       log_keepOpenAfterStamp, log_flushLevel
  use Logfile_interface, ONLY : Logfile_close, Logfile_open

  implicit none

#include "constants.h"
  character(len=*), intent(in)                      :: tag
  character(len=*), intent(in), OPTIONAL             :: attrib
  integer, intent(in)                               :: len
  character(len=*), dimension(len),intent(in)       :: strArr
  character(len=256)                     :: temp
  integer                                :: i
  character(len=40), save                :: dateStr  
  integer :: logUnit
  logical :: logUnitLocal=.false.

  if (log_globalMe .ne. MASTER_PE) return  ! only master processor writes to logfile
  
  if (.not. log_fileOpen) then
     call Logfile_open(logUnit,logUnitLocal)
  else
     logUnit = log_lun
  end if
  
  call current_date_time (dateStr)
  
  do i=1, len
     write (temp(1:), "(A, A, A, A, A, A, A)") ' [ ', dateStr(1:len_trim(dateStr)), ' ] ', tag , ": ", &
          &  trim(strArr(i))
     if (log_fileOpen) write (log_lun, '(A)') trim(temp)
  end do


  if (log_fileOpen .AND. .NOT. log_keepOpenAfterStamp) then
     call Logfile_close()
  else if (log_fileOpen .AND. log_flushLevel > 0) then
     FLASH_FLUSH(logUnit)
  end if

  return    
  
end subroutine Logfile_stampStrArray

subroutine Logfile_stampStrPair( strArr, len, dimen, tag, attrib)

  use Logfile_data, ONLY : log_globalMe,  log_fileOpen, log_lun, &
       log_keepOpenAfterStamp, log_flushLevel
  use Logfile_interface, ONLY : Logfile_close, Logfile_open

implicit none
#include "constants.h"

  character(len=*), intent(in)             :: tag
  character(len=*), OPTIONAL, intent(in)   :: attrib
  integer, intent(in)                      :: len, dimen
  character(len=*),dimension(len,dimen),intent(in)  :: strArr

  character(len=40), save                  :: dateStr
  integer, parameter                       :: temp_len=256
  character(len=temp_len)                  :: temp
  integer                                  :: i, j, cStrLen
  integer :: logUnit
  logical :: logUnitLocal=.false.

  if (log_globalMe .ne. MASTER_PE) return  ! only master processor writes to logfile
  

  if (.not. log_fileOpen) then
     call Logfile_open(logUnit,logUnitLocal)
  else
     logUnit = log_lun
  end if
  
  call current_date_time (dateStr)
  
  ! calculate compound string length from array
  cStrLen = 0
  do i=1, len
     do j=1, dimen
        if (j>2 .AND. j.GE.dimen-2) cycle
        cStrLen = cStrLen + len_trim(adjustl(strArr(i,j)))
     end do
  end do

  cStrLen = cStrLen + (len * 2) - 1 ! calculate space for data delimiter)

  ! it is assumed that the array has two dimensions, i.e. name and value        
  if (cStrLen > MAX_STRING_LENGTH) then
     do i=1, len
        if (len_trim(strArr(i,1)) > 0 .and. len_trim(strArr(i,2)) > 0) then
           write (temp(1:), "(A)") ' [ '//dateStr(1:len_trim(dateStr))//' ] '//trim(adjustl(tag))//": "//&
                &trim(strArr(i,1))//"="//trim(strArr(i,2))
        else if (len_trim(strArr(i,1)) > 0) then
           write (temp(1:), "(A)") ' [ '//dateStr(1:len_trim(dateStr))//' ] '//trim(adjustl(tag))&
                &//": "//trim(strArr(i,1))
        else if  (len_trim(strArr(i,2)) > 0) then
           write (temp(1:), "(A)") ' [ '//dateStr(1:len_trim(dateStr))//' ] '//trim(adjustl(tag))&
                &//": "//trim(strArr(i,2))
        else
           temp = ""
        end if
        if (len_trim(adjustl(temp)) > 0) then
           if (log_fileOpen) write (log_lun, '(A)') trim(temp)
        end if
     end do
  else
     write (temp(1:), "(A)") ' [ '//dateStr(1:len_trim(dateStr))//' ] '//trim(adjustl(tag))//": "
     ! it is assumed that the array has at least two dimensions, i.e. name and value
     do i=1, len
        if (len_trim(temp) + 1 + len_trim(strArr(i,1)) + len_trim(strArr(i,2)) + &
             min(1,len_trim(strArr(i,1)))*min(1,len_trim(strArr(i,2))) > temp_len) then
           if (log_fileOpen) write (log_lun, '(A)') trim(temp) !flush temp
           write (temp(1:), "(A)") ' [ '//dateStr(1:len_trim(dateStr))//' ] '//trim(adjustl(tag))//"... "
        end if
        if (len_trim(strArr(i,1)) > 0 .and. len_trim(strArr(i,2)) > 0) then
           write (temp(len_trim(temp)+2:), "(A, A, A)") & ! +2 is needed to keep white space at end
                & trim(strArr(i,1)), "=", trim(strArr(i,2))
        else if (len_trim(strArr(i,1)) > 0) then
           write (temp(len_trim(temp)+2:), "(A)") & ! +2 is needed to keep white space at end
                & trim(strArr(i,1))
        else if  (len_trim(strArr(i,2)) > 0) then
           write (temp(len_trim(temp)+2:), "(A)") & ! +2 is needed to keep white space at end
                & trim(strArr(i,2))
        else
           temp = ""
        end if
        do j=3,dimen
           if (len_trim(temp)+len_trim(strArr(i,j)) .LE. temp_len) &
                write (temp(len_trim(temp)+1:), "(A)") strArr(i,j) !addtl. delimiters or similar
        end do
        if ( i .lt. len .AND. len_trim(temp) < temp_len) write (temp(len_trim(temp)+1:), "(A)") " "
     end do
     if (len_trim(adjustl(temp)) > 0) then
        if (log_fileOpen) write (log_lun, '(A)') trim(temp)
     end if
  end if


  if (log_fileOpen .AND. .NOT. log_keepOpenAfterStamp) then
     call Logfile_close()
  else if (log_fileOpen .AND. log_flushLevel > 0) then
     FLASH_FLUSH(logUnit)
  end if
  
  return    
  
end subroutine Logfile_stampStrPair
