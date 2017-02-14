/*
 source/monitors/Profiler/ProfilerMain/mpihpm/pr_prof_control

 NAME
  pr_prof_control

 SYNOPSIS
  pr_prof_control(const int mode)

 DESCRIPTION
  Turns gprof or vprof data collection on and off.
  If you intend to use gprof then you must link FLASH using -pg
  for this subroutine to do anything.

  Both BG/P and BG/Q provide libmpihpm which includes MPI profiling
  and hardware counter profiling.  It can also include CPU profiling
  through gprof (BG/P) or vprof (BG/Q) and will dump CPU profiling data
  for MPI task 0, the MPI task with the maximum communication time, the
  MPI task with the minimum communication time and the MPI task with
  the median communication time.  I check whether I am on BG/P or
  BG/Q using the __PPC64__ predefined macro.

 ARGUMENTS
  mode - turns profiling on or off
*/

#include "Profiler.h"

void pr_prof_control(const int mode)
{
/* The __PPC64__ macro is defined by the xlc compiler
   but not the xlf compiler. */
#ifdef __PPC64__
  switch (mode) {
  case PRF_START_PROFILER:
    vprof_start(); /* void vprof_start(void) */
    break;
  case PRF_STOP_PROFILER:
    vprof_stop(); /* void vprof_stop(void) */
    break;
  default:
    break;
  }
#else
  switch (mode) {
  case PRF_START_PROFILER:
    moncontrol(1);
    break;
  case PRF_STOP_PROFILER:
  case PRF_DISABLE_PROFILER:
    moncontrol(0);
    break;
  default:
    break;
  }
#endif
}
