!!****if* source/flashUtilities/contiguousConversion/ut_convertToArrayIndicies
!!
!! NAME
!!
!!  ut_convertToArrayIndicies
!!
!! SYNOPSIS
!!
!!  ut_convertToArrayIndicies(integer, intent(IN) :: dims, 
!!                            integer, intent(IN) :: memoryOffset,
!!                            integer, dimension(1:dims), intent(IN) :: arrayLBound, 
!!                            integer, dimension(1:dims), intent(IN) :: arrayUBound,
!!                            integer, dimension(1:dims), intent(OUT) :: elementCoord)
!!
!! DESCRIPTION
!!
!!  Routine to convert a memory offset into integer indicies 
!!  of a hypothetical N-dimensional contiguous array.  Here,
!!  the memory offset is an input parameter which holds 
!!  the offset relative to the start of the N-dimensional array.
!!  e.g. suppose we have the following two dimensional array: myArray(1:2,1:4), and 
!!  we want to know the indicies corresponding to a memory offset of 5.  In this 
!!  case the subroutine returns the indicies (2,3) because the contiguous
!!  indicy coordinates are (1,1), (2,1), (1,2), (2,2), (1,3), (2,3).  Here, 
!!  (1,1) is the starting coordinate, so it corresponds to memory offset 0.
!!
!! ARGUMENTS
!!     
!!  dims : Dimensionality of the hypothetical array.
!!  memoryOffset : Memory offset from the start of the hypothetical array.  This
!!                 is then converted into integer indicies in the subroutine.  
!!                 The memory offset is required in terms of the datatype of 
!!                 the array, as opposed to bytes.
!!  arrayLBound : Array containing the lower bounds of the hypothetical array.
!!  arrayUBound : Array containing the upper bounds of the hypothetical array.
!!  elementCoord : The returned integer indicies.
!!
!!  
!!***
!*******************************************************************************

subroutine ut_convertToArrayIndicies(dims, memoryOffset, arrayLBound, arrayUBound, elementCoord)

  implicit none
  integer, intent(IN) :: dims, memoryOffset
  integer, dimension(1:dims), intent(IN) :: arrayLBound, arrayUBound
  integer, dimension(1:dims), intent(OUT) :: elementCoord

  integer, dimension(1:dims) :: elementCount, elementsToNextIndicy
  integer :: i, j, memoryOffsetRemainder

  !This is a count of the number of memory elements between successive array elements.
  elementsToNextIndicy(1) = 1
  do i = 2, dims, 1
     elementsToNextIndicy(i) = elementsToNextIndicy(i-1) * (arrayUbound(i-1) - arrayLbound(i-1) + 1)
  end do

  memoryOffsetRemainder = memoryOffset
  do j = dims, 2, -1
     elementCount(j) = (memoryOffsetRemainder / elementsToNextIndicy(j))
     memoryOffsetRemainder = memoryOffsetRemainder - (elementCount(j) * elementsToNextIndicy(j))
  end do
  elementCount(1) = memoryOffsetRemainder


#ifdef DEBUG_MODE
  do i = 1, dims, 1
     if (elementCount(i) > (arrayUbound(i) - arrayLbound(i) + 1)) then
        print *, "[ConvertToArrayIndicies]: Invalid data passed!"
        return
     end if
  end do
#endif

  !Add the array lower bound so that the caller can pass a non-unit based array.
  elementCoord(1:dims) = elementCount(1:dims) + arrayLBound(1:dims)
  
end subroutine ut_convertToArrayIndicies
