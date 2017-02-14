!!****f* source/Driver/Driver_driftBlock
!!
!! NAME
!!  Driver_driftBlock
!!
!! DESCRIPTION
!!
!!  Register the content hash of the provided block data, if it is different
!!  from the last time this block was hashed then log it to file.
!!
!! ARGUMENTS
!!
!!  src_file: source file location to log in case of changed hash
!!  src_line: source line location to log in case of changed hash
!!  blk: processor local block number
!!  ptr: block data to hash, all values will be hashed so this should not include guard cells
!!  gds: grid data struct type of data
!!
!!***
subroutine Driver_driftBlock(src_file, src_line, blk, ptr, gds)
  implicit none
  character(len=*), intent(in) :: src_file
  integer, intent(in) :: src_line
  integer, intent(in) :: blk
  real, intent(in) :: ptr(:,:,:,:)
  integer, intent(in) :: gds
end subroutine Driver_driftBlock
