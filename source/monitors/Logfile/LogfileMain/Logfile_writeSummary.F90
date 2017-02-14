!!****if* source/monitors/Logfile/LogfileMain/Logfile_writeSummary
!!
!! NAME
!!   Logfile_writeSummary
!!
!! SYNOPSIS
!!
!!   Logfile_writeSummary(character(len=*):: strArr(length,dim),
!!                        integer(in)     :: length,
!!                        integer(in)     :: dim,
!!                        integer(in)     :: strLen,
!!                        integer(in)     :: numHeaders,
!!                   logical(in),optional :: reduced,
!!                   logical(in),optional :: separateFiles)
!!
!! DESCRIPTION
!!
!!   Logfile_writeSummary writes the data handed to it by the Timers unit.
!!   Logfile_writeSummary does not do any of the performance calculations.  Its only
!!   role is to format the data handed to it by Timers and write the data neatly to the
!!   Logfile.
!!
!!   The summary implemented here has two sections: a header section
!!   and a timers summary listing section.  The header section has two columns
!!   in each row, a name and a value, and there are numHeaders rows of
!!   these.  strArr(1:numHeaders,1) hold the names, and strArr(1:numHeaders,2) holds the
!!   values.  strArr(numHeaders+1, :) holds strings that are the column names 
!!   for the timers summary listing section.  strArr(numHeaders+2:, 1) holds the 
!!   amount of indentation for each timer, and strArr(numHeadears+2:, 2:) holds 
!!   the data for each of the columns in timers summary. 
!!
!!   There are two kinds of summaries that may be generated:
!!    - the "traditional" summary, which contains timing information collected only
!!      for the processor writing to the log file, normally processor 0; and
!!    - a summary with timing information "reduced" across all processors.
!!
!! ARGUMENTS
!!
!!   strArr  - array holding all the run summary information (like evolved zones, seconds
!!             in monitoring period, etc).  The first rows hold header information,
!!             the next hold timer summary information. 
!!   length  - first dimension of strArr; the number of headers + the number 
!!             of lines in the summary, + 1 for the summary column names
!!   dim     - the number of columns in the timer summary + 1 for the indentation 
!!             of the timers
!!   strLen  - length of each string entry (likely MAX_STRING_LENGTH)
!!   numHeaders - the number of name/value pairs in the header of the summary
!!   reduced - if present and .TRUE., generate summary of reduced timer data;
!!             otherwise generate a normal local-processor summary.
!!   separateFiles - if true, every processor writes its summary to its own file named
!!                   timer_summary_<process id>
!!
!! EXAMPLE
!!
!!  A typical one-processor summary will look like this:
!!
!! ==============================================================================
!! perf_summary: code performance summary
!!                      beginning : 03-29-2006  19:23.18
!!                         ending : 03-29-2006  19:23.25
!!   seconds in monitoring period :                6.870
!!         number of subintervals :                   11
!!        number of evolved zones :                15680
!!               zones per second :             2282.359
!! ------------------------------------------------------------------------------
!! accounting unit                       time sec  num calls   secs avg  time pct
!! ------------------------------------------------------------------------------
!! initialization                          0.469      1           0.469     6.828
!!  guardcell internal                     0.104      8           0.013     1.515
!! evolution                               6.400      1           6.400    93.158
!!  hydro                                  5.943     20           0.297    86.507
!! ...
!!
!! NOTES
!!
!!  The user will likely never call this routine directly.  Developers will need to
!!  understand it if output format to logfile is changed.
!!
!!***

subroutine Logfile_writeSummary(strArr, length, dim, strLen, numHeaders, reduced, separateFiles)

  use Logfile_data, ONLY : log_globalMe,  log_lun, log_fileOpen   
  use Logfile_interface, ONLY : Logfile_break, Logfile_close, &
    Logfile_open

  implicit none

#include "constants.h"

  integer, intent(in)                                  :: length, dim, strLen, numHeaders
  character(len=MAX_STRING_LENGTH), intent(in), dimension(length,dim)  :: strArr
  logical, optional, intent(IN)                        :: reduced
  logical, optional, intent(IN)                        :: separateFiles
  character(len=MAX_STRING_LENGTH)                     :: indentStr, tag
  
  integer  :: i, tmpLen
  logical  :: doreduced, doseparate
  integer  :: summary_lun
  integer :: logUnit
  logical :: logUnitLocal=.false.

  if (present(reduced)) then
     doreduced = reduced
  else
     doreduced = .FALSE.
  end if

  if (present(separateFiles)) then
     doseparate = separateFiles
  else
     doseparate = .FALSE.
  end if

  ! if separate, everyone writes to his own file

  indentStr(:) = " "

  if (doseparate) then
     summary_lun = 931
     call log_openSummaryFile(log_globalMe, summary_lun)
     call log_summaryBreak(log_globalMe, summary_lun, "=")
  else if (log_globalMe .eq. MASTER_PE) then
     if(.not. log_fileOpen) then
        call Logfile_open(logUnit,logUnitLocal)
     end if

     call Logfile_break( "=")
     summary_lun = log_lun
  end if

  if (doseparate .or. log_globalMe .eq. MASTER_PE) then
     if (doreduced .and. log_globalMe .eq. MASTER_PE) then
        write (summary_lun, *) 'perf_summary: code performance summary statistics'
     else
        write (summary_lun, *) 'perf_summary: code performance summary for process ', log_globalMe
     end if
        
     do i=1, length
        if (i <= numHeaders) then
           
           write (summary_lun, "(1X, A30, ' : ', A20)") trim(adjustl(strArr(i,1))), trim(adjustl(strArr(i,2)))
        else if (i == numHeaders + 1) then
           call log_summaryBreak(log_globalMe, summary_lun, "-")
           if (.NOT. doreduced) then
              write (summary_lun, "(1X, A15, 23X, A8, 2X, A9, 3X, A8, 2X, A8)") trim(adjustl(strArr(i,2))),&
                   & trim(adjustl(strArr(i,3))), &
                   & trim(adjustl(strArr(i,4))), trim(adjustl(strArr(i,5))), trim(adjustl(strArr(i,6)))
           else if (log_globalMe .eq. MASTER_PE) then
              write (summary_lun, "(1X, A15, 18X, A12, 2X, A12, 1X, A12, 2X, A10)") trim(adjustl(strArr(i,2))),&
                   & trim(adjustl(strArr(i,3))), &
                   & trim(adjustl(strArr(i,4))), trim(adjustl(strArr(i,5))), trim(adjustl(strArr(i,6)))
           end if
           call log_summaryBreak(log_globalMe, summary_lun, "-")
        else
           tmpLen = len_trim(adjustl(strArr(i,1)))
           write(tag,*) indentStr(1:tmpLen) // trim(adjustl(strArr(i,2)))
           if (.NOT. doreduced) then
              write (summary_lun, "(A28, 6X, A12, 1X, A6, 1X, A15, 2X, A8)") tag, &
                   & trim(adjustl(strArr(i,3))), &
                   & trim(adjustl(strArr(i,4))), trim(adjustl(strArr(i,5))), trim(adjustl(strArr(i,6)))
           else if (log_globalMe .eq. MASTER_PE) then
              write (summary_lun, "(A28, 6X, A12, 2X, A12, 1X, A12, 2X, A10)") tag, &
                   & trim(adjustl(strArr(i,3))), &
                   & trim(adjustl(strArr(i,4))), trim(adjustl(strArr(i,5))), trim(adjustl(strArr(i,6)))
           end if
        end if
     end do
        
     call log_summaryBreak(log_globalMe, summary_lun, "=")    
     if (doseparate) then
        close(summary_lun)
     else
        call Logfile_close()
     end if
  end if

  
end subroutine Logfile_writeSummary


subroutine log_openSummaryFile(log_globalMe, summary_lun)

  implicit none

  integer, intent(in) :: log_globalMe
  integer, intent(in) :: summary_lun
  
  character(len=MAX_STRING_LENGTH), save :: fileName       ! Name of log file
  integer :: attempts
  integer :: io_status
  
  write (fileName, '(A,I6.6)') 'timer_summary_', log_globalMe 

  attempts = 0
  io_status = 1
  do while ((io_status /= 0) .and. (attempts < 10))
     open (summary_lun, file=fileName, position='append', iostat=io_status)
     attempts = attempts + 1
  enddo
  if (io_status /= 0) then
     write (*,*) 'log_openSummaryFile:  Error: cannot open file; ', fileName, &
          'io_status = ', io_status
  end if

end subroutine log_openSummaryFile


subroutine log_summaryBreak (log_globalMe, summary_lun, char)

  implicit none

#include "constants.h"
  
  integer, intent(in)                :: log_globalMe
  integer, intent(in)                :: summary_lun
  character(len=1), intent(in)                   :: char
  character(len=MAX_STRING_LENGTH-2) :: buff
  
  character(len=1)                   :: lchar
  !logical                            :: do_open = .false.

  if (ichar(char) < 32) then
     lchar = ' '
  else
     lchar = char
  end if
  buff = repeat(lchar, MAX_STRING_LENGTH-2) ! 78 characters

  write(summary_lun, *) buff

  return
end subroutine log_summaryBreak
