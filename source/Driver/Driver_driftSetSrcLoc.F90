!!****f* source/Driver/Driver_driftSetSrcLoc
!!
!! NAME
!!  Driver_driftSetSrcLoc
!!
!! DESCRIPTION
!!
!!  Saves the current source file & line to module vars in Drift_data.
!!  It is recommended that you call this method like so:
!!    call Driver_driftSetSrcLoc(__FILE__,__LINE__)
!!
!! ARGUMENTS
!!
!!  filename: source file location to log in case of changed hash
!!  line: source line location to log in case of changed hash
!!
!!***
subroutine Driver_driftSetSrcLoc(filename, line)
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: line
end subroutine Driver_driftSetSrcLoc
