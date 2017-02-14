!!****f* source/physics/sourceTerms/Flame/Flame_step
!!
!! NAME
!!
!!  Flame_step
!!
!! SYNOPSIS
!!
!!  call Flame_step ( integer(in)                       :: num_blocks,
!!                    integer(in),dimension(num_blocks) :: blockList(:),
!!                    real(in)                          :: dt )
!!
!! DESCRIPTION
!!
!!  Performs a flame step.
!!
!!  Pricipal public function for Flame unit.
!!  Evolve flame forward for all blocks in blockList by one step
!!  of size dt.
!!  Flame speed and flame effects sub-units are called within this
!!  subroutine as they are required
!!  for ADR Applies unsplit reaction-diffusion operater to FLAM_MSCALAR
!!  May or may not deposit energy, depending on which
!!  Flame_Effects module had been included
!!
!! ARGUMENTS
!!
!!   num_blocks - number of blocks.
!!    blockList - the block list.
!!           dt - the time step.
!!
!! SEE ALSO
!!
!!  see Flame_interface.F90 for possible updates
!!
!!***

! this is a stub for when the Flame Unit is not included
!
! Dean Townsley 2008
!
subroutine Flame_step( num_blocks, blockList, dt  )    
       
  implicit none
  integer, INTENT(in)                        :: num_blocks
  integer, INTENT(in), DIMENSION(num_blocks) :: blockList
  real,    INTENT(in)                        :: dt

  return
end subroutine Flame_step
