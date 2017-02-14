program test
  implicit none  
  real, parameter :: val1 = -0.0000001, val2 = 0.0
  real :: ret

  !NaN caught by -qflttrap=enable:invalid
  print *, achar(10), " ******** Testing sqrt of a negative number ********"
  call sqrt_error(val1, ret)

  !NaN caught by -qflttrap=enable:invalid and not linking
  !against libmass!
  print *, achar(10), " ******** Testing log of a negative number ********"
  call log_error(val1, ret)

  !Division by zero caught by -qflttrap=enable:zerodivide
  print *, achar(10), " ******** Testing divide by zero ********"
  call divide(val2, ret)

  !Quiet NaN caught by -qflttrap=enable:nanq
  !Signaling NaN caught by -qflttrap=enable:invalid
  print *, achar(10), " ******** Testing uninitialized data in automatic variable ********"
  call automatic_variable_error(ret)

  print *, achar(10), " ******** Testing uninitialized data in allocatable array ********"
  call allocatable_error(ret)
end program test

subroutine sqrt_error(x,y)
  real, intent(IN) :: x
  real, intent(OUT) :: y
  y = sqrt(x)
end subroutine sqrt_error

subroutine log_error(x,y)
  real, intent(IN) :: x
  real, intent(OUT) :: y
  y = log(x)
end subroutine log_error

!I create a subroutine which divides to prevent the compiler
!from replacing an exception-producing calculation with a
!constant infinity value.  (If I wanted to catch a compile
!time y = 1.0 / 0.0 then I would need to compile with
!-qfloat=nofold)
subroutine divide(x,y)
  real, intent(IN) :: x
  real, intent(OUT) :: y
  y = 1.0 / x
end subroutine divide

subroutine automatic_variable_error(y)
  implicit none
  real, intent(OUT) :: y
  real :: x
  y = x * 2.0
  print *, "variable after assignment from automatic variable = ", y
end subroutine automatic_variable_error

subroutine allocatable_error(y)
  implicit none
  real, intent(OUT) :: y
  real, allocatable, dimension(:) :: x
  allocate(x(10))
  y = x(3) * 2.0
  deallocate(x)
  print *, "variable after assignment from allocatable array = ", y
end subroutine allocatable_error
