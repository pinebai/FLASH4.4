!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_createRayTags
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

  implicit none

  integer, optional, intent (in) :: passSplitDriver

  return
end subroutine ed_createRayTags
