!!****if* source/Particles/ParticlesMain/Particles_getGlobalNum
!!
!! NAME
!!  Particles_getGlobalNum
!!
!! SYNOPSIS
!!
!!  Particles_getGlobalNum(integer(OUT)  :: globalNumParticles)
!!                
!! DESCRIPTION 
!!
!!  Returns the global (total) number of particles in the 
!!  simulation after computing it by calling an MPI reduction routine.
!!
!!  Needed for input/output across all processors.
!!
!! ARGUMENTS
!!
!!  globalNumParticles:     integer        number of particles across all processors
!!
!! SIDE EFFECTS
!!
!!  Writes to the log file a line with the current global number of particles
!!  when first called, and then whenever the number has changed since the previous
!!  time.
!!
!!***

subroutine Particles_getGlobalNum(globalNumParticles)

  use Logfile_interface, ONLY: Logfile_stamp
  use Particles_data, ONLY : pt_numLocal, pt_meshMe, useParticles, pt_meshComm
#define LOG_STATS .TRUE.

  implicit none
#include "mpif.h"

  integer, intent(out)  :: globalNumParticles

  integer :: ierr
  integer,save :: lastLoggedNum = -1

  if (useParticles) then
     call MPI_AllReduce(pt_numLocal, globalNumParticles, 1, MPI_INTEGER, &
          MPI_SUM, pt_meshComm, ierr)
  else
     globalNumParticles = 0 !NOTE: pt_meshMe must be set for Logfile_stamp.
  end if

  if (LOG_STATS) then
     if (globalNumParticles .NE. lastLoggedNum) then
        call Logfile_stamp(globalNumParticles, "[Particles_getGlobalNum]: Number of particles now")
        lastLoggedNum = globalNumParticles
     end if
  end if
  return
end subroutine Particles_getGlobalNum
