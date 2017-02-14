!!****if* source/flashUtilities/contiguousConversion/ut_convertToMemoryOffset
!!
!! NAME
!!
!!  ut_convertToMemoryOffset
!!
!! SYNOPSIS
!!
!!  ut_convertToMemoryOffset(integer, intent(IN) :: dims, 
!!                           integer, dimension(1:dims), intent(IN) :: elementCoord),
!!                           integer, dimension(1:dims), intent(IN) :: arrayLBound, 
!!                           integer, dimension(1:dims), intent(IN) :: arrayUBound,
!!                           integer, intent(OUT) :: memoryOffset)
!!
!! DESCRIPTION
!!
!!  Routine to convert integer indicies of a hypothetical N-dimensional 
!!  contiguous array into a memory offset in terms of the underlying datatype.  
!!  Here, the memory offset is an input parameter which holds 
!!  the offset relative to the start of the N-dimensional array.
!!  e.g. suppose we have the following two dimensional array: myArray(1:2,1:4), and 
!!  we want to know the memory offset corresponding to the indicies (2,3).  In this 
!!  case the subroutine returns the memory offset 5 because the contiguous
!!  indicy coordinates are (1,1), (2,1), (1,2), (2,2), (1,3), (2,3).  Here, 
!!  (1,1) is the starting coordinate, so it corresponds to memory offset 0.
!!
!! ARGUMENTS
!!     
!!  dims : Dimensionality of the hypothetical array.
!!  elementCoord : The input integer indicies.
!!  arrayLBound : Array containing the lower bounds of the hypothetical array.
!!  arrayUBound : Array containing the upper bounds of the hypothetical array.
!!  memoryOffset : Memory offset from the start of the hypothetical array. 
!!                 The memory offset is given in terms of the datatype of 
!!                 the array, as opposed to bytes.
!!  
!!***
!*******************************************************************************

subroutine ut_convertToMemoryOffset(dims, elementCoord, arrayLBound, arrayUBound, memoryOffset)

implicit none
integer, intent(IN) :: dims
integer, dimension(1:dims), intent(IN) :: elementCoord, arrayUBound, arrayLBound
integer, intent(OUT) :: memoryOffset
integer, dimension(1:dims) :: elementsToNextIndicy
integer :: i

#ifdef DEBUG_MODE
do i = 1, dims, 1
   if (elementCoord(i) < arrayLBound(i) .or. elementCoord(i) > arrayUBound(i)) then
      print *, "[ConvertToMemoryOffset]: Invalid data passed!"
      return
   end if
end do
#endif

!This is a count of the number of memory elements between successive array elements.
elementsToNextIndicy(1) = 1
do i = 2, dims, 1
   elementsToNextIndicy(i) = elementsToNextIndicy(i-1) * (arrayUbound(i-1) - arrayLbound(i-1) + 1)
end do

!Make a cumulative count of the number of memory elements to reach the user-specified 
!element in each dimension.  Notice that the lower limit of the target array is 
!taken into account so that the caller can pass a non-unit based array.  
memoryOffset = sum( (elementCoord(1:dims) - arrayLBound(1:dims)) * elementsToNextIndicy(1:dims) )

end subroutine ut_convertToMemoryOffset
