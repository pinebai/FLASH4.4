module Driver_interface
  implicit none
contains
  subroutine Driver_abortFlash (errorMessage)
    implicit none
    character(len=*), intent(in) :: errorMessage
    print *, errorMessage
    stop     
    return
  end subroutine Driver_abortFlash
end module Driver_interface
