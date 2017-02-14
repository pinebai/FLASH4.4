!!****if* source/Grid/GridMain/paramesh/gr_markDerefineByTime
!!
!! NAME
!!
!!  gr_markDerefineByTime
!!
!! SYNOPSIS
!!
!!  call gr_markDerefineByTime()
!!
!! DESCRIPTION
!!
!!  Modifies the global gr_maxRefine which is meant to be a user-
!!  changeable version of lrefine_max.
!!
!!  Some PARAMESH versions don't take kindly to user code modifying
!!  lrefine_max itself.  So enter gr_maxRefine.  It is supposed to
!!  be < lrefine_max in order to be meaningful.  The AMR refinement
!!  criteria have to cooperate, by checking gr_maxRefine and setting
!!  derefine flags / cancelling refine flags accordingly.
!!
!! ARGUMENTS
!!
!!
!! SEE ALSO
!!   gr_markRefineDerefine
!!   Grid_markRefineDerefine
!!   
!! HISTORY
!!   June 2009   created  KW, based on logic by Sean Couch
!!***

subroutine gr_markDerefineByTime()
#include "constants.h"
  use Grid_data, ONLY : gr_maxRefine,&
                        gr_lrefineMaxRedTimeScale,&
                        gr_lrefineMaxRedTRef,&
                        gr_lrefineMaxRedLogBase, gr_meshMe
  use tree, ONLY : lrefine_max, lrefine_min
  use Driver_interface, ONLY : Driver_getSimTime
  use logfile_interface, ONLY : Logfile_stamp
  implicit none

  real :: time
  integer       :: refine_steps, prevMaxRefine

  call Driver_getSimTime(time)

  if (time > gr_lrefineMaxRedTRef) then
     refine_steps = max(0,int( &
          log10((time - gr_lrefineMaxRedTRef)/gr_lrefineMaxRedTimeScale) &
          /log10(gr_lrefineMaxRedLogBase) &
          + 1))
     prevMaxRefine = gr_maxRefine
     gr_maxRefine = max(lrefine_min,lrefine_max - refine_steps)
     if (gr_maxRefine .NE. prevMaxRefine) then
        if (gr_meshMe == MASTER_PE) &
             print*,'[gr_markDerefineByTime] now the maximum refine is ',gr_maxRefine
        call Logfile_stamp(gr_maxRefine,'[gr_markDerefineByTime] set gr_maxRefine')
     end if

  end if

end subroutine gr_markDerefineByTime
