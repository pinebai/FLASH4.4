!!****if* source/monitors/Debugger/DebuggerMain/Debugger_init
!!
!! NAME
!!  Debugger_init
!!
!! SYNOPSIS
!!
!!  Debugger_init()
!!                   
!!  
!! DESCRIPTION 
!!  
!!  Initialize the debugger unit
!!  
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Debugger_init()
  use Debugger_data, ONLY : dbg_globalMe, dbg_doHeapCheck
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype
  implicit none

  call Driver_getMype(GLOBAL_COMM,dbg_globalMe)

  call RuntimeParameters_get("doHeapCheck", dbg_doHeapCheck)
#if defined(DMALLOC) || defined(DDT)
  if (dbg_doHeapCheck) then
     call dbg_enable_heap_check()
  else
     call dbg_disable_heap_check()
  end if
#else
  if (dbg_globalMe == MASTER_PE) then
     print *, "*** DMALLOC / DDT support not compiled into this application ***"
  end if
#endif
end subroutine Debugger_init
