!!****f* source/physics/sourceTerms/Turb/Turb_calc
!!
!! NAME
!!
!!  Turb_calc
!!
!! SYNOPSIS
!!
!!  call Turb_init( integer(in)                       :: num_blocks,
!!                  integer(in),dimension(num_blocks) :: blockList(:) )
!!
!! DESCRIPTION
!!
!!  Calculate the laplacian of the velocity field and store
!!  in TURB_VAR, LAPY_VAR, and LAPZ_VAR for each block, then
!!  fill guard cells.
!!  Calculate the curl of the laplacian for each block, then
!!  fill guard cells again.
!!  This procedure gives OP2 in Colin et al. (2000)
!!
!! ARGUMENTS
!!
!!   num_blocks - number of blocks.
!!    blockList - the block list.
!!
!! SEE ALSO
!!
!!  See Turb_interface.F90 for possible updates
!!
!!***

! Stub for calculation of turbulent intensity for TFI implementation
!
! Aaron Jackson 2010
subroutine Turb_calc(num_blocks, blockList)
  implicit none
  integer, intent(in)                        :: num_blocks
  integer, intent(in), dimension(num_blocks) :: blockList

  return
end subroutine Turb_calc
