!!****ih* source/flashUtilities/contiguousConversion/ut_conversionInterface
!!
!! This is the public interface file for the array conversion  
!! subroutines.
!!***
module ut_conversionInterface

  implicit none      

  interface
     subroutine ut_convertToMemoryOffset(dims, elementCoord, arrayLBound, arrayUBound, memoryOffset)
       implicit none
       integer, intent(IN) :: dims
       integer, dimension(1:dims), intent(IN) :: elementCoord, arrayUBound, arrayLBound
       integer, intent(OUT) :: memoryOffset
     end subroutine ut_convertToMemoryOffset
  end interface

  interface
     subroutine ut_convertToArrayIndicies(dims, memoryOffset, arrayLBound, arrayUBound, elementCoord)
       implicit none
       integer, intent(IN) :: dims, memoryOffset
       integer, dimension(1:dims), intent(IN) :: arrayLBound, arrayUBound
       integer, dimension(1:dims), intent(OUT) :: elementCoord
     end subroutine ut_convertToArrayIndicies
  end interface

end module ut_conversionInterface
