!!****if* source/Grid/GridMain/paramesh/gr_setMaxRefineByTime
!!
!! NAME
!!
!!  gr_setMaxRefineByTime
!!
!! SYNOPSIS
!!
!!  call gr_setMaxRefineByTime()
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!
!! SEE ALSO
!!   gr_markRefineDerefine
!!   Grid_markRefineDerefine
!!   
!!***
subroutine gr_setMaxRefineByTime()
#include "constants.h"
#include "Flash.h"
  use Grid_data, ONLY : gr_maxRefine, &
       gr_meshMe, &
       gr_lrefinemaxByTime, &
       gr_lrefmaxTimes, &
       gr_lrefmaxTimeValues
  use tree, ONLY : lrefine_max, lrefine_min
  use Driver_interface, ONLY : Driver_getSimTime
  use logfile_interface, ONLY : Logfile_stamp
  implicit none

  integer :: i, prevMaxRefine
  real :: time

  call Driver_getSimTime(time)

  prevMaxRefine = gr_maxRefine
  do i = 1, GR_LREFMAXTIMES
     if(time >= gr_lrefmaxTimes(i) .and. gr_lrefmaxTimes(i) > 0.0) then
        gr_maxRefine = max(lrefine_min, gr_lrefmaxTimeValues(i))
     end if
  end do

  if (prevMaxRefine /= gr_maxRefine) then
     if(gr_meshMe == MASTER_PE) then
        print*,'[gr_setMaxRefineByTime] now the maximum refine is ', &
             gr_maxRefine, ' was ', prevMaxRefine
     end if
     call Logfile_stamp(gr_maxRefine,'[gr_markDerefineByTime] set gr_maxRefine')        
  end if

end subroutine gr_setMaxRefineByTime
