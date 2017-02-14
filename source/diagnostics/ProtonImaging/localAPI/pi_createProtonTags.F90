!!****if* source/diagnostics/ProtonImaging/localAPI/pi_createProtonTags
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
!!  Creates a unique global tag number for each proton. This is achieved by
!!  performing the following three steps:
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
!!     3) Global unique tag assignment --> Assign consecutive tag values
!!        starting from 1 for all protons present on all processors at this
!!        moment. Using MPI_Allgather, the same is done as for step 2) but
!!        at processor level, i.e. adding up all number of protons on the
!!        processors lower than the processor being analyzed.
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

  implicit none

  return
end subroutine pi_createProtonTags
