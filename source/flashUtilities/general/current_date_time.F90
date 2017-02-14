!!****if* source/flashUtilities/general/current_date_time
!!
!! NAME
!!
!!  current_date_time
!!
!!
!! SYNOPSIS
!!
!!  current_date_time(character(len=40),(OUT) :: date_string)
!!
!!
!!
!! DESCRIPTION
!!
!!  return the current date and time nicely formatted in a single character
!!  string.  The argument date_string is a 40 character string containing
!!  the date and time
!!
!!
!! ARGUMENTS
!!
!!  date_string - returned string with date and time
!!
!!
!!
!!***

subroutine current_date_time(date_string)

  implicit none
  
  character(len=*),intent(out) :: date_string
  
  character (len=8)  :: current_date
  character (len=10) :: current_time
  character (len=5)  :: time_zone

  integer   :: date_values(8)


  character (len=2) :: month_string, day_string, hour_string, &
       minute_string, seconds_string

  character (len=4) :: year_string
  
  call date_and_time(current_date, current_time, time_zone, date_values)

  write (year_string, '(i4.4)') date_values(1)
  write (month_string, '(i2.2)') date_values(2)
  write (day_string, '(i2.2)') date_values(3)
  write (hour_string, '(i2.2)') date_values(5)
  write (minute_string, '(i2.2)') date_values(6)
  write (seconds_string, '(i2.2)') date_values(7)

  date_string = month_string // '-' // day_string // '-' //          &
       year_string // '  ' // hour_string // ':' // minute_string // &
       ':' // seconds_string // '.' // current_time(8:10)


  return
end subroutine current_date_time
