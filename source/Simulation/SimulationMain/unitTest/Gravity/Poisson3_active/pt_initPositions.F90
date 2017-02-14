!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3_active/pt_initPositions
!!
!! NAME
!!
!!   pt_initPositions
!!
!! SYNOPSIS
!!
!!   pt_initPositions(integer, INTENT(in) :: blockId, 
!!                    logical, INTENT(out) :: success)
!!
!! DESCRIPTION
!!   An override for the modified poisson3 problem.
!!
!!
!! ARGUMENTS
!!
!!   blockId : Id of block in current processor
!!   success: logical argument indicating whether particles initalised 
!!            correctly on this block.
!!
!!
!!***

subroutine pt_initPositions(blockId, success)

  use Simulation_data, ONLY : sim_threshold
  use Particles_data, ONLY : pt_numLocal, pt_maxPerProc, pt_numProcs, pt_meshMe
  use pt_interface, ONLY : pt_createTag
  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, &
       Grid_getBlkCenterCoords, Grid_getPointData, Grid_getBlkIndexLimits, &
       Grid_getSingleCellCoords, Grid_putPointData, Grid_getSingleCellVol, Grid_getDeltas
  use Grid_data, ONLY : gr_meshMe

#include "Flash.h"
#include "constants.h"

  implicit none

  integer, intent(in) :: blockId
  logical, intent(out) :: success

  integer       :: i, j, k
  real          :: xr, yr, zr, vxr, vyr, vzr
  integer, save :: p
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real, dimension(MDIM) :: cellCoords
  integer, dimension(MDIM) :: point
  real :: densityToLose, cellDensity, massToGiveParticle, cellVolume
  real :: del(MDIM)
  real    :: urp  ! random number generator
  integer :: seed_size
  logical, save :: first_call = .true.
  print *, "Processor", gr_meshMe, "entering the customised pt_initPositions: pt_maxPerProc=", pt_maxPerProc


  ! Particle slot number (incremented and saved between calls)
  p = pt_numLocal

  !-------------------------------------------------------------------------
  ! Initialise the random number generator.
  if(first_call.eqv..true.) then     
     !  returned value seed_size gives the number of integers the processor uses for the
     !    starting value
     call random_seed(SIZE=seed_size)

     !  generates a large (from pt_pRand) and unique integer for each processor

     i = int(0.5 * gr_meshMe) + pt_numProcs

     !  initializes the random number vector with a fixed seed (from i)
     call random_seed(PUT=(/(i, j = 1, seed_size)/))

     !  now call the random_number generator the same number of times as has been seeded.
     !  WHY we call the random_number generator so many times this is unknown -- perhaps to
     !    really mix things up?  If so, then why a fixed seed?
     do j = 1,2*i
        !  returns a single uniform random number urp between zero and
        call random_number (harvest=urp)
     end do
     first_call = .false.
  end if
  !-------------------------------------------------------------------------


  !We need to loop over each cell in each block.  
  !If the density is above a threshold, then we reduce the 
  !density in that cell to the threshold.  The lost density 
  !is then multiplied by the volume of the cell to give the 
  !lost mass.  This mass is then used as the mass for a
  !newly initilialised particle.

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getDeltas(blockId, del)

  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     point(3) = k
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        point(2) = j
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           point(1) = i

           !Extract the density that has been placed in each cell.
           call Grid_getPointData(blockID, CENTER, DENS_VAR, EXTERIOR, point, cellDensity)


           !Create a particle if cell density is above a threshold.
           if (cellDensity > sim_threshold) then
              densityToLose = cellDensity - sim_threshold

              call Grid_getSingleCellVol(blockID, EXTERIOR, point, cellVolume)
              massToGiveParticle = densityToLose * cellVolume


              !Store reduced density in the grid
              call Grid_putPointData(blockID, CENTER, DENS_VAR, EXTERIOR, point, sim_threshold)


              !Now we initialise a single particle, which is given the lost mass.

              p = p + 1  !Add one to total particle number.

              !Get the left edge position of the fluid. 
              xr = 0; yr = 0; zr = 0;
              call Grid_getSingleCellCoords(point, blockId, LEFT_EDGE, EXTERIOR, cellCoords)

              !I can't place the particle smack bang on the central point, because in 
              !this case the weighting factor specifies all the particle mass is 
              !mapped to only the single central cell.  Part of the reason for this 
              !test is to check prolongation and restriction in particle mapping.
              call random_number (harvest=urp)
              xr = cellCoords(IAXIS)+(urp*del(IAXIS))

              call random_number (harvest=urp)
              yr = cellCoords(JAXIS)+(urp*del(JAXIS))

              call random_number (harvest=urp)
              zr = cellCoords(KAXIS)+(urp*del(KAXIS))


              !Get the cell centered current velocity of the fluid.
              vxr = 0; vyr = 0; vzr = 0
              call Grid_getPointData(blockID, CENTER, VELX_VAR, EXTERIOR, point, vxr)
              call Grid_getPointData(blockID, CENTER, VELY_VAR, EXTERIOR, point, vyr)
              call Grid_getPointData(blockID, CENTER, VELZ_VAR, EXTERIOR, point, vzr)


              call InitSingleParticle(p, xr, yr, zr, vxr, vyr, vzr, massToGiveParticle, blockId)

           end if

        end do
     end do
  end do


  ! Set the particle database local number of particles.
  pt_numLocal = p
  print *, "Processor", gr_meshMe, "intialised", pt_numLocal, "particles"
  success = .true.

  return
end subroutine pt_initPositions


!**********************************************************************
!  Routine:     InitSingleParticle

!  Description: Initialize a single particle with given characteristics.


subroutine InitSingleParticle (p, xpos, ypos, zpos, xvel, yvel, zvel, &
                               mass, block)

  use Particles_data, ONLY : pt_maxPerProc, particles
  use Driver_interface, ONLY : Driver_abortFlash

#include "Flash.h"
#include "constants.h"

  implicit none

  real, intent(in)    :: xpos, ypos, zpos, xvel, yvel, zvel, mass
  integer, intent(in) :: p, block

  if (p > pt_maxPerProc) &
    call Driver_abortFlash("InitSingleParticle:  Exceeded max # of particles!")

  particles(BLK_PART_PROP,p)  = real(block)
  particles(PROC_PART_PROP,p) = real(pt_meshMe)
  particles(MASS_PART_PROP,p) = mass
  particles(POSX_PART_PROP,p) = xpos
  particles(VELX_PART_PROP,p) = xvel
  particles(POSY_PART_PROP,p) = ypos
  particles(VELY_PART_PROP,p) = yvel
  particles(POSZ_PART_PROP,p) = zpos
  particles(VELZ_PART_PROP,p) = zvel

return
end subroutine InitSingleParticle
