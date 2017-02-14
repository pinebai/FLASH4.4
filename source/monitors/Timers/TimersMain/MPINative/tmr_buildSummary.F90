!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_buildSummary
!!
!! NAME
!!
!!   tmr_buildSummary
!!
!! SYNOPSIS
!!
!!  call tmr_buildSummary(character, len=MAX_STRING_LENGTH, dimension(length,columns,intent(inout) ::summaryArray, 
!!                        integer, intent(in)         :: length, 
!!                        integer, intent(in)         :: columns, 
!!                        integer, intent(inout)      :: currentIndex, 
!!                        integer, intent(in)         :: indentation, 
!!                        type(tmr_stack), intent(in) :: rootStack,
!!                        logical, intent(in)         :: reduced)
!!
!!
!! DESCRIPTION
!!
!!   A recursive routine that takes as input the summary output string
!!   data, a root timer stack and goes through the Timers runtime
!!   timer stacks and builds the summary output.  This is an
!!   inherently recursive process-- what you do at each level of a
!!   stack is what you want to do after you add a level to that stack.
!!   Each time this routine is called with a stack of timer names, it
!!   goes through the fill list of timer names (in the order in which
!!   they were created) and checks if adding the current timer name to
!!   the current stack matches an entry in the list of all timer
!!   stacks that occured during execution of Flash, and if it does, it
!!   adds that entry into the summary, and calls itself with the new
!!   root stack as the one it just found by adding the current timer
!!   name.
!!
!!   The output is kindof strange, but it is the summary data as
!!   strings in summaryArray.  The first index of summaryArray
!!   represents a row of output, and the second, a column, more or
!!   less.  The first column is the exception-- summaryArray(i, 1) is
!!   a number, n, of a non-blank characters, in this case, 'c', which
!!   represents that the output of row i should be indented n
!!   characters.  This is how the logfile knows how to indent the
!!   timer name at the beginning of each row of output to create the
!!   effect of showing the timer in the summary.  The second column,
!!   summaryArray(i, 2) is the name of the timer, summaryArray(i,3) is
!!   the amount of time spent in the timer, summaryArray(i,4) is the
!!   number of times it was started, summaryArray(i,5) is the time
!!   spent / times called, and summaryArray(i,6) is the percentage of
!!   total runtime spent in this timer.  If you want to change the
!!   nature of the summary, you would need to change this file, and
!!   the logfile code that labels the columns in the logfile.
!!
!!   If reduced is .TRUE., this subroutine does the MPI reduce calls for
!!   generating min, max, and average times across all processors and
!!   formats results differently.
!!
!! ARGUMENTS
!!
!! summaryArray -- see above
!! length       -- first index size of array summaryArray 
!! columns      -- second index size of array summaryArray
!! currentIndex -- the row in summaryArray to add to next
!! indentation  -- how much we're already indented when we enter this routine
!! rootStack    -- which timer stack we're starting at upon entry
!! reduced      -- If true, produce timer summary information reduced across
!!                 processors, instead of the default one-processor summary.
!!
!! NOTES
!!  When the reduced flag is used, this subroutine may hang, if the set of
!!  Timers_start and Timers_stop calls executed on each processor is not the same.
!!  This may happen in particular when there are fewer (leaf) blocks than processors.
!!***

recursive subroutine tmr_buildSummary(summaryArray, length, columns, currentIndex, indentation, rootStack,reduced)

  use Timers_data, ONLY: tmr_stack, tmr_numSegments, tmr_acctSegs, tmr_numSegments,&
       tmr_globalNumProcs, tmr_globalComm
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "constants.h"
#include "Flash_mpi.h"
  
  integer, intent(in) :: length
  integer, intent(in) :: columns
  character(len=MAX_STRING_LENGTH), dimension(length, columns), intent(inout) :: summaryArray
  integer, intent(inout)  :: currentIndex
  integer, intent(in)     :: indentation
  type(tmr_stack), intent(in) :: rootStack
  logical, intent(IN)     :: reduced
  
  integer i, pushResult
  integer :: indicies(tmr_numSegments)
  real :: maxTime, minTime, avgTime
  integer ::  ierr

  type(tmr_stack) :: currentStack
  logical :: doreduced
  doreduced = reduced

  do i = 1, tmr_numSegments
     call tmr_stackListIndex(tmr_acctSegs(i)%stacks, rootStack, indicies(i))
     if (indicies(i) > 0) then
        if (doreduced) then

           !Initializing avgTime to 0.0 stops FPE on myPE /= MASTER_PE 
           !when evaluating avgTime = avgTime / tmr_globalNumProcs expression.
           maxTime = 0.0; minTime = 0.0; avgTime = 0.0

           call MPI_Reduce(tmr_acctSegs(i)%time(indicies(i)), maxTime, 1, &
                FLASH_REAL, MPI_MAX,     &
                MASTER_PE, tmr_globalComm, ierr)
           call MPI_Reduce(tmr_acctSegs(i)%time(indicies(i)), minTime, 1, &
                FLASH_REAL, MPI_MIN,     &
                MASTER_PE, tmr_globalComm, ierr)
           call MPI_Reduce(tmr_acctSegs(i)%time(indicies(i)), avgTime, 1, &
                FLASH_REAL, MPI_SUM,     &
                MASTER_PE, tmr_globalComm, ierr)
           avgTime = avgTime / tmr_globalNumProcs
        end if

        write (summaryArray(currentIndex, 1), *) repeat("c", indentation)
        write (summaryArray(currentIndex, 2), "(A)") trim(tmr_acctSegs(i)%name(:min(len_trim(tmr_acctSegs(i)%name),25)))


        if (doreduced) then
           write (summaryArray(currentIndex,3), "(F15.3)") maxTime
           write (summaryArray(currentIndex,4), "(F15.3)") minTime
           write (summaryArray(currentIndex,5), "(F15.3)") avgTime
           write (summaryArray(currentIndex,6), "(I10)") tmr_acctSegs(i)%timesCalled(indicies(i))
        else
           write (summaryArray(currentIndex,3), "(F15.3)") tmr_acctSegs(i)%time(indicies(i))
           write (summaryArray(currentIndex,4), "(I6)") tmr_acctSegs(i)%timesCalled(indicies(i))
           write (summaryArray(currentIndex,5), "(F15.3)") tmr_acctSegs(i)%time(indicies(i)) &
                / tmr_acctSegs(i)%timesCalled(indicies(i))
           write (summaryArray(currentIndex,6), "(F8.3)") tmr_acctSegs(i)%pctTime(indicies(i))
        end if
        
        currentIndex = currentIndex + 1
        call tmr_stackAssign(currentStack, tmr_acctSegs(i)%stacks%stacks(indicies(i)))
        call tmr_stackPush(currentStack, i, pushResult)
        if (pushResult < 0) then
           call Driver_abortFlash("tmr_buildSummary: ran out of space building summary")
        end if
        call tmr_buildSummary(summaryArray, length, columns, currentIndex, indentation+1, currentStack, doreduced)
     end if
  end do
end subroutine tmr_buildSummary
