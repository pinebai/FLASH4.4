!!****f* source/monitors/Logfile/Logfile_getDateTimeStr
!!
!! NAME
!!  Logfile_getDateTimeStr
!!
!! SYNOPSIS
!!  call Logfile_getDateTimeStr(character(len=28)(OUT) :: dateTimeStr)
!!
!!
!! DESCRIPTION
!!  Get a string that contains the current date and time in the same
!!  format used in Logfie_stamp messages.
!!
!!  This interface is meant for use by code outside the Logfile unit
!!  that wants to generate messages in the format of the Logfile unit
!!  but circumventing  the Logfile_stamp interface.
!!
!!
!! ARGUMENTS
!!  dateTimeStr - the date and time string is returned here.
!!
!! NOTES
!!  Currently, the caller could just call current_date_time directly
!!  to the same effect as calling this subroutine, with the following
!!  difference:
!!  o The string returned by this implementation is enclosed in
!!    [ square brackets ].
!!
!!***

subroutine Logfile_getDateTimeStr(dateTimeStr)

  implicit none

  character(len=28), intent(OUT)        :: dateTimeStr


  dateTimeStr = ' '
  

end subroutine Logfile_getDateTimeStr
