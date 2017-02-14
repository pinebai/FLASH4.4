subroutine ut_qsortInt(inOutArray, arrayLength, ascOrderArg)
  implicit none
  integer, dimension(:), intent(INOUT) :: inOutArray
  integer, intent(IN) :: arrayLength
  logical, optional, intent(IN) :: ascOrderArg
  logical :: ascOrder

  if (present(ascOrderArg)) then
     ascOrder = ascOrderArg
  else
     ascOrder = .true.
  end if

  if (ascOrder .eqv. .true.) then
     call ut_qsort_int_asc(inOutArray, arrayLength)
  else
     call ut_qsort_int_Desc(inOutArray, arrayLength)  
  end if

end subroutine ut_qsortInt


subroutine ut_qsortFloat(inOutArray, arrayLength, ascOrderArg)
  implicit none
  real*4, dimension(:), intent(INOUT) :: inOutArray
  integer, intent(IN) :: arrayLength
  logical, optional, intent(IN) :: ascOrderArg
  logical :: ascOrder

  if (present(ascOrderArg)) then
     ascOrder = ascOrderArg
  else
     ascOrder = .true.
  end if

  if (ascOrder .eqv. .true.) then
     call ut_qsort_float_asc(inOutArray, arrayLength)
  else
     call ut_qsort_float_desc(inOutArray, arrayLength)  
  end if

end subroutine ut_qsortFloat


subroutine ut_qsortDouble(inOutArray, arrayLength, ascOrderArg)
  implicit none
  real*8, dimension(:), intent(INOUT) :: inOutArray
  integer, intent(IN) :: arrayLength
  logical, optional, intent(IN) :: ascOrderArg
  logical :: ascOrder

  if (present(ascOrderArg)) then
     ascOrder = ascOrderArg
  else
     ascOrder = .true.
  end if

  if (ascOrder .eqv. .true.) then
     call ut_qsort_double_asc(inOutArray, arrayLength)
  else
     call ut_qsort_double_desc(inOutArray, arrayLength)  
  end if

end subroutine ut_qsortDouble
