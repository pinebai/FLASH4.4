pure real function sqrt3(val)
  implicit none
  real, intent(IN) :: val
  real, parameter :: third = 1.0e0/3.0e0
  sqrt3 = val ** third
end function sqrt3
