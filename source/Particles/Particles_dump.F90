!!****f* source/Particles/Particles_dump
!!
!! NAME
!!    Particles_dump
!!
!! SYNOPSIS
!!    Particles_dump()
!!                   integer(IN) :: blockCount,
!!                   integer(IN), dimension(:)  :: blockList,
!!                   integer(IN) :: nstep,
!!                   real(IN)    :: time)
!!
!! DESCRIPTION
!!
!!   Dump particle information to a plain old file.  The output file is called
!!      test_ParticlesDump_00##  where the ## relates to the processor number.
!!   Quick and dirty output implemented for StirTurb BG/L testing.  You can see an
!!      example of the output by running the simulation unitTest/Particles.
!!   On the positive side, there is an associated fidlr3/IDL routine particles_dump.pro which will
!!      read the input from this file.
!!
!! ARGUMENTS
!!
!!  
!!  blockCount:             integer(IN)     number of blocks within this processor
!!  blockList(blockCount):  integer(IN)     block IDs within this processor
!!  nstep:                  integer(IN)     current time step index
!!  time:                   real(IN)        current simulation time
!!
!! PARAMETERS
!!
!! NOTES
!!
!!  The source for the real routine is in the ParticlesMain/unitTest subdirectory
!!
!!***


!-------------------------------------------------------------------
subroutine Particles_dump(blockCount,blockList,nstep,time,dt)

#include "Flash.h"
#include "constants.h"

  implicit none

  
  integer, intent(IN) :: blockCount
  integer, intent(IN) :: blockList(blockCount)
  integer, intent(IN) :: nstep
  real, intent(IN)    :: time, dt
 
  
end subroutine Particles_dump

