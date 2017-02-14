Module ut_qsortInterface
  implicit none  
  interface ut_qsort
    subroutine ut_qsortInt(inOutArray, arrayLength, ascOrderArg)
      implicit none
      integer, dimension(:), intent(INOUT) :: inOutArray
      integer, intent(IN) :: arrayLength
      logical, optional, intent(IN) :: ascOrderArg
    end subroutine ut_qsortInt

    subroutine ut_qsortFloat(inOutArray, arrayLength, ascOrderArg)
      implicit none
      real*4, dimension(:), intent(INOUT) :: inOutArray
      integer, intent(IN) :: arrayLength
      logical, optional, intent(IN) :: ascOrderArg
    end subroutine ut_qsortFloat

    subroutine ut_qsortDouble(inOutArray, arrayLength, ascOrderArg)
      implicit none
      real*8, dimension(:), intent(INOUT) :: inOutArray
      integer, intent(IN) :: arrayLength
      logical, optional, intent(IN) :: ascOrderArg
    end subroutine ut_qsortDouble
 end interface ut_qsort
end Module ut_qsortInterface
