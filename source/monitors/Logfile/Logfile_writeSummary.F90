!!****f* source/monitors/Logfile/Logfile_writeSummary
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



subroutine Logfile_writeSummary( strArr, length, dim, strLen, numHeaders, reduced, separateFiles)

#include "constants.h"

  implicit none
  integer, intent(in)                      :: length, dim, strLen, numHeaders
  character(len=MAX_STRING_LENGTH), dimension(length,dim), intent(in)  :: strArr
  logical, optional, intent(IN)                        :: reduced
  logical, optional, intent(IN)                        :: separateFiles
 

end subroutine Logfile_writeSummary

