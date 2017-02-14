!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonsCreate/pi_createProtonTags
!!
!! NAME
!!
!!  pi_createProtonTags
!!
!! SYNOPSIS
!!
!!  call pi_createProtonTags ()
!!
!! DESCRIPTION
!!
!!  Creates a unique global tag number for each proton created during the entire
!!  time step cycle. This is achieved by performing the following four steps:
!!
!!     1) Block unique tag assignment --> Assign consecutive tag values
!!        starting from 1 in each block. After this step, the tags of the
!!        protons in each LEAF block will have duplicate values.
!!
!!     2) Processor unique tag assignment --> Assign consecutive tag values
!!        starting from 1 for all protons on the current processor. After this
!!        step, the tags are unique at processor level but duplicated among
!!        processors. This step is done by looping over all blocks on the
!!        processor adding up all number of protons on the blocks lower than
!!        the block being analyzed.
!!
!!     3) Batch global unique tag assignment --> Assign consecutive tag values
!!        starting from 1 for all protons present on all processors for the
!!        current batch. Using MPI_Allgather, the same is done as for step 2)
!!        but at processor level, i.e. adding up all number of protons on the
!!        processors lower than the processor being analyzed.
!!
!!     4) Overall global unique tag adjustment --> Adjust consecutive tag values
!!        such that all protons created during the entire simulation will have
!!        a unique tag. This is simply done by adding the value of the previous
!!        maximum tag value. After this step has been done, the maximum tag value is
!!        updated by the number of protons generated for the current batch on
!!        all processors.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!  Since the tags consist of consecutive integers, this method of tag
!!  generation will work for N protons in a simulation, where N is the
!!  maximum integer representable on the machine.
!!
!!***

subroutine pi_createProtonTags ()

  use ProtonImaging_data,  ONLY : pi_globalComm,     &
                                  pi_globalMe,       &
                                  pi_globalNumProcs, &
                                  pi_maxProtonCount, &
                                  pi_protonCount,    &
                                  pi_protons,        &
                                  pi_tagMax
  
  use Driver_interface,    ONLY : Driver_abortFlash
  use Grid_interface,      ONLY : Grid_getLocalNumBlks

  implicit none

#include "Flash.h"
#include "constants.h"
#include "ProtonImaging.h"

  include "Flash_mpi.h"

  integer :: block
  integer :: blockID
  integer :: error
  integer :: localNumBlocks
  integer :: myBlockTagOffset
  integer :: myProcBlockTagOffset
  integer :: myProcTagOffset
  integer :: proc
  integer :: proton

  integer, allocatable :: localNumProtons (:)     ! Number of protons on each process
  integer, allocatable :: myBlockProtons  (:)     ! Number of protons on each block
!
!
!     ...Get the number of blocks on the current processor and allocate the
!        necessary arrays.
!
!
  call Grid_getLocalNumBlks (localNumBlocks)

  allocate (myBlockProtons  (0:localNumBlocks ))
  allocate (localNumProtons (1:pi_globalNumProcs))
!
!
!     ...Step 1: Assign each proton a unique contiguous tag on this process.
!
!
  myBlockProtons (0:localNumBlocks) = 0

  do proton = 1, pi_protonCount
     blockID = int (pi_protons (PROTON_BLCK,proton))
     myBlockProtons (blockID) = myBlockProtons (blockID) + 1
     pi_protons (PROTON_TAGS,proton) = real (myBlockProtons (blockID))
  end do
!
!
!     ...Step 2: Calculate array that will determine block offsets for each block.
!
!  
  do block = 2,localNumBlocks
     myBlockProtons (block) = myBlockProtons (block) + myBlockProtons (block - 1)
  end do
!
!
!     ...Step 3: Calculate processor offset for current processor.
!
!
  call MPI_Allgather (pi_protonCount,  &       ! this is being sent from i-th process
                      1,               &       ! how many elements ?
                      MPI_INTEGER,     &       ! type of elements sent
                      localNumProtons, &       ! stored into i-th array element on all processes
                      1,               &       ! number of elements stored
                      MPI_INTEGER,     &       ! type of elements stored
                      pi_globalComm,   &       ! communicator handle
                      error            )       ! error handle

  myProcTagOffset = 0

  do proc = 1, pi_globalMe
     myProcTagOffset = myProcTagOffset + localNumProtons (proc)
  end do
!
!
!     ...Combine steps 1,2,3 and 4 to calculate the time step globally unique tag. 
!
!
  do proton = 1, pi_protonCount

     blockID = int (pi_protons (PROTON_BLCK,proton))

     myBlockTagOffset     = myBlockProtons (blockID - 1)
     myProcBlockTagOffset = myProcTagOffset + myBlockTagOffset

     pi_protons (PROTON_TAGS,proton) = pi_protons (PROTON_TAGS,proton) + real (pi_tagMax + myProcBlockTagOffset)

  end do
!
!
!     ...For next batch step 4: Update the maximum tag value by the number of protons
!        generated during the current batch.
!
!
  pi_tagMax = pi_tagMax + sum (localNumProtons (1:pi_globalNumProcs))
!
!
!     ...Deallocate arrays used. 
!
!
  deallocate (myBlockProtons )
  deallocate (localNumProtons)
!
!
!     ...Ready! 
!
!
  return
end subroutine pi_createProtonTags
