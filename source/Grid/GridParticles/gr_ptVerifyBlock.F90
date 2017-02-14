!!****if* source/Grid/GridParticles/gr_ptVerifyBlock
!!
!! NAME
!!  gr_ptVerifyBlock
!!
!! SYNOPSIS
!!
!!  gr_ptVerifyBlock(real(IN)    :: particles(:,:),
!!                   integer(IN) :: propCount,
!!                   integer(IN) :: localNumParticles,
!!                   integer(IN) :: maxParticlesPerProc)
!!                    
!!  
!! DESCRIPTION 
!!  
!! Used to check whether each particle is associated with:
!!   1.  A valid block (i.e. not NONEXISTENT)
!!   2.  The correct block (i.e. it exists within this block's bounding box).
!!
!! ARGUMENTS 
!!
!!  particles : List of particles. It is a 2D real array, the first dimension
!!              represents particle's properties, and second dimension is 
!!              index to particles.
!!
!!  propCount : The count of fields in the particles data structure
!!
!!  localNumParticles : The current number of particles mapped to this processor.
!!
!!  maxParticlesPerProc : This is parameter determined at runtime, 
!!                        and is the maximum number of local
!!                        particles that a simulation expects to have. 
!!                        All the arrays in the particles
!!                        unit are allocated based on this number
!!***

subroutine gr_ptVerifyBlock(particles,propCount,localNumParticles,maxParticlesPerProc)

  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_outsideBoundBox
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_meshMe
  use gr_ptData, ONLY : gr_ptBlk, gr_ptTag, gr_ptPosx, gr_ptPosy, gr_ptPosz
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: propCount
  integer,intent(IN) :: localNumParticles
  integer,intent(IN) :: maxParticlesPerProc
  real, dimension(propCount, maxParticlesPerProc),intent(IN) :: particles
  real, dimension(LOW:HIGH, MDIM) :: bndBox
  integer, dimension(MDIM) :: Negh
  integer :: i, blockID, logUnit
  logical :: outside, inErrorState, logUnitLocal

  inErrorState = .false.
  logUnitLocal = .true.

  !Does each particle exist in the correct bounding block.
  do i = 1, localNumParticles
     blockID = int(particles(gr_ptBlk, i))

     if (blockID == NONEXISTENT) then

        inErrorState = .true.
        call Logfile_open(logUnit,logUnitLocal)
        write(logUnit, *) "[gr_ptVerifyBlock]: particle", i, "tag", particles(gr_ptTag,i), &
             "has NONEXISTENT block" 
        write(logUnit,*)'and position is',particles(gr_ptPosx:gr_ptPosz,i)
        call Logfile_close(logUnitLocal)

     else

        call Grid_getBlkBoundBox(blockID,bndBox)
        if(particles(gr_ptPosx,i) < bndBox(LOW,IAXIS) .or. &
             particles(gr_ptPosx,i) >= bndBox(HIGH,IAXIS) .or. &
             particles(gr_ptPosy,i) < bndBox(LOW,JAXIS) .or. &
             particles(gr_ptPosy,i) >= bndBox(HIGH,JAXIS) .or. &
             particles(gr_ptPosz,i) < bndBox(LOW,KAXIS) .or. &
             particles(gr_ptPosz,i) >= bndBox(HIGH,KAXIS)) then

           inErrorState = .true.
           call Grid_outsideBoundBox(particles(gr_ptPosx:gr_ptPosz,i),bndBox,outside,Negh)

           call Logfile_open(logUnit,logUnitLocal)
           write(logUnit,*) "[gr_ptVerifyBlock]: particle", i, "tag", particles(gr_ptTag,i), &
                "at position:", particles(gr_ptPosx:gr_ptPosz,i), &
                "but block encloses bndBox(IAXIS):", bndBox(LOW:HIGH,IAXIS), &
                "bndBox(JAXIS):", bndBox(LOW:HIGH,JAXIS), &
                "bndBox(KAXIS):", bndBox(LOW:HIGH,KAXIS), &
                "Outside attribute:", outside, "Negh:", Negh
           call Logfile_close(logUnitLocal)
        end if

     end if
  end do

  if (inErrorState) then
     call Driver_abortFlash("Particle on wrong block.  See local error logfiles.")
  end if

end Subroutine gr_ptVerifyBlock
