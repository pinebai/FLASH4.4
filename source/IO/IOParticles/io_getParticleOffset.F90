!!****if* source/IO/IOParticles/io_getParticleOffset
!!
!! NAME
!!  io_getParticleOffset
!!
!! SYNOPSIS
!!
!!  io_getParticleOffset( integer(in)  :: localNumParticles,
!!                        integer(out) :: globalNumParticles,
!!                        integer(out) :: particleOffset)
!!               
!!  
!! DESCRIPTION 
!!
!!  Returns the global(total) number of particles in the 
!!  simulation as well as the particle offset.  DEV: This routine
!!  is redundant with Particles_getGlobalNum....
!!
!! ARGUMENTS
!!
!!  localNumParticles - number of particles per proc
!!  globalNumParticles - total number of particles in the simulation
!!  particleOffset - offset value of current processor's first particle
!!                   (if all processors are writing particles to one big array
!!                    proc1's particles get the first slots in the array,
!!                    proc2's particles get the next spots and so on)
!!
!! NOTES
!!
!!  !DEV: Old comment said: currently only implemented with the uniform grid!
!! 
!!
!!***

subroutine io_getParticleOffset( localNumParticles, globalNumParticles, particleOffset)

  use IO_data, ONLY : io_outputSplitNum, io_splitFileNum, io_globalMe, io_globalNumProcs

  implicit none

#include "Flash_mpi.h"
#include "constants.h"


  integer, intent(in)   ::  localNumParticles
  integer, intent(out)  :: globalNumParticles, particleOffset
  integer           :: partsToLeft(0:io_globalNumProcs-1)
  integer :: ierr, i
  



  !gather all local particles from all procs and put them into an array partsToLeft
  call MPI_Allgather(localNumParticles, 1, MPI_INTEGER, &
       partsToLeft, 1, MPI_INTEGER, &
       MPI_COMM_WORLD, ierr)
  


  globalNumParticles = 0
  

  do i = 0, io_globalNumProcs -1
     globalNumParticles = globalNumParticles + partsToLeft(i)
  end do
  
  do i = io_globalNumProcs-1,1,-1
     partsToLeft(i) = partsToLeft(i-1)
  end do
  
  partsToLeft(0) = 0

  do i = 2,io_globalNumProcs-1
     partsToLeft(i) = partsToLeft(i) + partsToLeft(i-1)
  end do
  
  particleOffset = partsToLeft(io_globalMe)
 

  return
end subroutine io_getParticleOffset
