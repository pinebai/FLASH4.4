!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRaysCreate/ed_createRayTags
!!
!! NAME
!!
!!  ed_createRayTags
!!
!! SYNOPSIS
!!
!!  call ed_createRayTags (integer (in), optional :: passSplitDriver)
!!
!! DESCRIPTION
!!
!!  Creates a unique global tag number for each ray. This is achieved by
!!  performing the following three steps:
!!
!!     1) Block unique tag assignment --> Assign consecutive tag values
!!        starting from 1 in each block. After this step, the tags of the
!!        rays in each LEAF block will have duplicate values.
!!
!!     2) Processor unique tag assignment --> Assign consecutive tag values
!!        starting from 1 for all rays on the current processor. After this
!!        step, the tags are unique at processor level but duplicated among
!!        processors. This step is done by looping over all blocks on the
!!        processor adding up all number of rays on the blocks lower than
!!        the block being analyzed.
!!
!!     3) Global unique tag assignment --> Assign consecutive tag values
!!        starting from 1 for all rays present on all processors at this
!!        moment. Using MPI_Allgather, the same is done as for step 2) but
!!        at processor level, i.e adding up all number of rays on the
!!        processors lower than the processor being analyzed.
!!
!!    When the split driver is used during a simulation, the tags are
!!    offset in the second half of the time step to ensure that each
!!    ray launched over the entire time-step has a unique tag. Thus,
!!    every ray can be uniquely identified by a cycle-tag pair.
!!
!! ARGUMENTS
!!
!!  passSplitDriver : indicates first/second half of time step for split driver
!!
!! NOTES
!!
!!  Since the tags consist of consecutive integers, This method of tag
!!  generation will work for N rays in a simulation, where N is the maximum
!!  integer representable on the machine.
!!
!!***

subroutine ed_createRayTags (passSplitDriver)

  use EnergyDeposition_data,  ONLY : ed_globalComm,     &
                                     ed_globalMe,       &
                                     ed_globalNumProcs, &
                                     ed_maxRayCount,    &
                                     ed_rayCount,       &
                                     ed_rays,           &
                                     ed_tagStart
  
  use Driver_interface,       ONLY : Driver_abortFlash
  use Grid_interface,         ONLY : Grid_getLocalNumBlks

  implicit none

#include "Flash.h"
#include "constants.h"
#include "EnergyDeposition.h"

  include "Flash_mpi.h"
  
  integer, optional, intent (in) :: passSplitDriver

  integer :: block
  integer :: blockID
  integer :: error
  integer :: localNumBlocks
  integer :: myBlockTagOffset
  integer :: myProcBlockTagOffset
  integer :: myProcTagOffset
  integer :: proc
  integer :: ray

  integer, allocatable :: localNumRays (:)     ! Number of rays on each process
  integer, allocatable :: myBlockRays  (:)     ! Number of rays on each block
!
!
!     ...Get the number of blocks on the current processor and allocate the
!        necessary arrays.
!
!
  call Grid_getLocalNumBlks (localNumBlocks)

  allocate (myBlockRays  (0:localNumBlocks ))
  allocate (localNumRays (1:ed_globalNumProcs))
!
!
!     ...Step 1: Assign each ray a unique contiguous tag on this process.
!
!
  myBlockRays (0:localNumBlocks) = 0

  do ray = 1, ed_rayCount
     blockID = int (ed_rays (RAY_BLCK,ray))
     myBlockRays (blockID)  = myBlockRays (blockID) + 1
     ed_rays (RAY_TAGS,ray) = real (myBlockRays (blockID))
  end do
!
!
!     ...Step 2: Calculate array that will determine block offsets for each block.
!
!  
  do block = 2,localNumBlocks
     myBlockRays (block) = myBlockRays (block) + myBlockRays (block - 1)
  end do
!
!
!     ...Step 3: Calculate processor offset for current processor.
!
!
  call MPI_Allgather (ed_rayCount,   &       ! this is being sent from i-th process
                      1,             &       ! how many elements ?
                      MPI_INTEGER,   &       ! type of elements sent
                      localNumRays,  &       ! stored into i-th array element on all processes
                      1,             &       ! number of elements stored
                      MPI_INTEGER,   &       ! type of elements stored
                      ed_globalComm, &       ! communicator handle
                      error          )       ! error handle

  myProcTagOffset = 0

  do proc = 1, ed_globalMe
     myProcTagOffset = myProcTagOffset + localNumRays (proc)
  end do
!
!
!     ...Combine steps 1,2 and 3 to calculate the globally unique tag. 
!
!
  do ray = 1, ed_rayCount

     blockID = int (ed_rays (RAY_BLCK,ray))

     myBlockTagOffset     = myBlockRays (blockID - 1)
     myProcBlockTagOffset = myProcTagOffset + myBlockTagOffset

     ed_rays (RAY_TAGS,ray) = ed_rays (RAY_TAGS,ray) + real (ed_tagStart + myProcBlockTagOffset)

  end do
!
!
!     ...Set tag offset for second half of time step. 
!
!
  ed_tagStart = 0

  if (present (passSplitDriver)) then
      if (passSplitDriver == 1) then
          do proc = 1, ed_globalNumProcs
             ed_tagStart = ed_tagStart + localNumRays (proc)
          end do
      end if
  end if
!
!
!     ...Deallocate arrays used. 
!
!
  deallocate (myBlockRays )
  deallocate (localNumRays)
!
!
!     ...Ready! 
!
!
  return
end subroutine ed_createRayTags
