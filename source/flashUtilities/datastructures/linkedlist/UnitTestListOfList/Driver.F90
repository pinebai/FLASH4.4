module Driver
  implicit none
  contains

  subroutine Driver_abortFlash(msg)
    implicit none
    character (len=*), intent(IN) :: msg
    print *, "ERROR!!!", msg
    stop
  end subroutine Driver_abortFlash
end module Driver
