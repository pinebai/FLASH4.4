!!****f* source/Driver/Driver_driftUnk
!!
!! NAME
!!  Driver_driftUnk
!!
!! DESCRIPTION
!!
!!  Compute one hash per unk variable over all registered block hashes.  Vars with
!!  hashes that have changed since the previous call will be logged.
!!
!! ARGUMENTS
!!
!!  src_file: source file location to log in case of changed hash
!!  src_line: source line location to log in case of changed hash
!!  flags: bitmask of options (for now there is only one flag)
!!    - DRIFT_NO_PARENTS: if present then only leaf blocks are included
!!
!!***
subroutine Driver_driftUnk(src_file, src_line, flags)
  implicit none
  character(len=*), intent(in) :: src_file
  integer, intent(in) :: src_line
  integer, intent(in), optional :: flags
end subroutine Driver_driftUnk
