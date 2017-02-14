!This program tests the contiguous conversion routines 
!independently from FLASH.  This is possible because the 
!conversion routines are re-entrant and have all data 
!passed by the caller.

!For testing/example purposes only.

module globalFuncAndParms

use ut_conversionInterface, ONLY : ut_convertToMemoryOffset, & 
     ut_convertToArrayIndicies

integer, parameter :: SUCCESS_VALUE = 0
integer, parameter :: FAILURE_VALUE = 1

end module globalFuncAndParms




Program UnitTest_ContiguousConversion

use globalFuncAndParms, ONLY : SUCCESS_VALUE, FAILURE_VALUE

implicit none
interface
   integer function test1d()
   end function test1d
end interface
interface
   integer function test2d()
   end function test2d
end interface
interface
   integer function test3d()
   end function test3d
end interface

integer :: iResult
character(len=*), parameter :: stringSuccess = "SUCCESS"
character(len=*), parameter :: stringFailure = "FAILURE"

!Assumes SUCCESS_VALUE=0 and FAILURE_VALUE=1.
character(len=100), dimension(SUCCESS_VALUE:FAILURE_VALUE) :: stringReturn
stringReturn(SUCCESS_VALUE) = stringSuccess
stringReturn(FAILURE_VALUE) = stringFailure


!Run each of the unit tests and report SUCCESS or FAILURE.
print *, "TESTING THE CONTIGUOUS CONVERSION ROUTINES"
iResult = test1d()
print *, " * Result of running test1d.... "//stringReturn(iResult)

iResult = test2d()
print *, " * Result of running test2d.... "//stringReturn(iResult)

iResult = test3d()
print *, " * Result of running test3d.... "//stringReturn(iResult)
print *, "DONE"


End Program UnitTest_ContiguousConversion



! ------------------------------------------------------------- !
!                        1D unit test                           !
! ------------------------------------------------------------- !
integer function test1d()

use globalFuncAndParms, ONLY : SUCCESS_VALUE, FAILURE_VALUE, &
     ut_convertToMemoryOffset, ut_convertToArrayIndicies

implicit none
integer, parameter :: dims = 1
integer, dimension(1:dims) :: elementCoord1D, newElementCoord1D
integer, dimension(3:7) :: arr1D
integer :: i, memoryOffset


!Manually fill arr1D with the expected memory offsets:
arr1D(3) = 0
arr1D(4) = 1
arr1D(5) = 2
arr1D(6) = 3
arr1D(7) = 4


!Loop through contiguous memory locations.
do i = lbound(arr1D,1), ubound(arr1D,1)
      
   elementCoord1D(1) = i


   !We want to cycle through each element of the array, and calculate 
   !each element's memory offset from the first element of the array.
   call ut_convertToMemoryOffset(dims, elementCoord1D, lbound(arr1D), ubound(arr1D), memoryOffset)
      
   if (memoryOffset /= arr1D(i)) then
      test1d = FAILURE_VALUE
      return
   end if

   !print *, "coordinates:", elementCoord1D, "memory offset=", memoryOffset

   !Now that we have the memory offset we would like to be able 
   !to restore the array coordinates.
   call ut_convertToArrayIndicies(dims, memoryOffset, lbound(arr1D), ubound(arr1D), newElementCoord1D)
            
   if ( (elementCoord1D(1) /= newElementCoord1D(1)) ) then
      test1d = FAILURE_VALUE
      return
   end if
      
end do

test1d = SUCCESS_VALUE

end function test1d



! ------------------------------------------------------------- !
!                        2D unit test                           !
! ------------------------------------------------------------- !
integer function test2d()

use globalFuncAndParms, ONLY : SUCCESS_VALUE, FAILURE_VALUE, &
     ut_convertToMemoryOffset, ut_convertToArrayIndicies

implicit none
integer, parameter :: dims = 2
integer, dimension(1:dims) :: elementCoord2D, newElementCoord2D
integer, dimension(1:2, 1:3) :: arr2D
integer :: i, j, memoryOffset


!Manually fill arr2D with the expected memory offsets:
arr2D(1, 1) = 0
arr2D(2, 1) = 1
arr2D(1, 2) = 2
arr2D(2, 2) = 3
arr2D(1, 3) = 4
arr2D(2, 3) = 5


!Loop through contiguous memory locations.
do j = lbound(arr2D,2), ubound(arr2D,2)
   do i = lbound(arr2D,1), ubound(arr2D,1)
      
      elementCoord2D(1) = i; elementCoord2D(2) = j


      !We want to cycle through each element of the array, and calculate 
      !each element's memory offset from the first element of the array.
      call ut_convertToMemoryOffset(dims, elementCoord2D, lbound(arr2D), ubound(arr2D), memoryOffset)
      
      if (memoryOffset /= arr2D(i,j)) then
         test2d = FAILURE_VALUE
         return
      end if

      !print *, "coordinates:", elementCoord2D, "memory offset=", memoryOffset

      !Now that we have the memory offset we would like to be able 
      !to restore the array coordinates.
      call ut_convertToArrayIndicies(dims, memoryOffset, lbound(arr2D), ubound(arr2D), newElementCoord2D)
            
      if ( (elementCoord2D(1) /= newElementCoord2D(1)) .or. &
           (elementCoord2D(2) /= newElementCoord2D(2))) then
         test2d = FAILURE_VALUE
         return
      end if
      

   end do
end do


!Loop through non-contiguous memory locations.
do i = lbound(arr2D,1), ubound(arr2D,1)
   do j = lbound(arr2D,2), ubound(arr2D,2)
      
      elementCoord2D(1) = i; elementCoord2D(2) = j


      !We want to cycle through each element of the array, and calculate 
      !each element's memory offset from the first element of the array.
      call ut_convertToMemoryOffset(dims, elementCoord2D, lbound(arr2D), ubound(arr2D), memoryOffset)
      
      if (memoryOffset /= arr2D(i,j)) then
         test2d = FAILURE_VALUE
         return
      end if

      !print *, "coordinates:", elementCoord2D, "memory offset=", memoryOffset

      !Now that we have the memory offset we would like to be able 
      !to restore the array coordinates.
      call ut_convertToArrayIndicies(dims, memoryOffset, lbound(arr2D), ubound(arr2D), newElementCoord2D)
            
      if ( (elementCoord2D(1) /= newElementCoord2D(1)) .or. &
           (elementCoord2D(2) /= newElementCoord2D(2))) then
         test2d = FAILURE_VALUE
         return
      end if
      

   end do
end do

test2d = SUCCESS_VALUE

end function test2d



! ------------------------------------------------------------- !
!                        3D unit test                           !
! ------------------------------------------------------------- !
integer function test3d()

use globalFuncAndParms, ONLY : SUCCESS_VALUE, FAILURE_VALUE, &
     ut_convertToMemoryOffset, ut_convertToArrayIndicies

implicit none
integer, parameter :: dims = 3
integer, dimension(1:dims) :: elementCoord3D, newElementCoord3D
integer, dimension(-2:0, 0:1, 3:4) :: arr3D
integer :: i, j, k, memoryOffset


!Manually fill arr3D with the expected memory offsets:
arr3D(-2, 0, 3) = 0
arr3D(-1, 0, 3) = 1
arr3D(0, 0, 3) = 2
arr3D(-2, 1, 3) = 3
arr3D(-1, 1, 3) = 4
arr3D(0, 1, 3) = 5
arr3D(-2, 0, 4) = 6
arr3D(-1, 0, 4) = 7
arr3D(0, 0, 4) = 8
arr3D(-2, 1, 4) = 9
arr3D(-1, 1, 4) = 10
arr3D(0, 1, 4) = 11



!Loop through contiguous memory locations.
do k = lbound(arr3D,3), ubound(arr3D,3)
   do j = lbound(arr3D,2), ubound(arr3D,2)
      do i = lbound(arr3D,1), ubound(arr3D,1)
      
         elementCoord3D(1) = i; elementCoord3D(2) = j; elementCoord3d(3) = k


         !We want to cycle through each element of the array, and calculate 
         !each element's memory offset from the first element of the array.
         call ut_convertToMemoryOffset(dims, elementCoord3D, lbound(arr3D), ubound(arr3D), memoryOffset)
      
         if (memoryOffset /= arr3D(i,j,k)) then
            test3d = FAILURE_VALUE
            return
         end if

         !print *, "coordinates:", elementCoord3D, "memory offset=", memoryOffset

         !Now that we have the memory offset we would like to be able 
         !to restore the array coordinates.
         call ut_convertToArrayIndicies(dims, memoryOffset, lbound(arr3D), ubound(arr3D), newElementCoord3D)

         if ( (elementCoord3D(1) /= newElementCoord3D(1)) .or. &
              (elementCoord3D(2) /= newElementCoord3D(2)) .or. &
              (elementCoord3D(3) /= newElementCoord3D(3))) then
            test3d = FAILURE_VALUE
            return
         end if


      end do
   end do
end do



!Loop through non-contiguous memory locations.
do j = lbound(arr3D,2), ubound(arr3D,2)
   do i = lbound(arr3D,1), ubound(arr3D,1)
      do k = lbound(arr3D,3), ubound(arr3D,3)
      
         elementCoord3D(1) = i; elementCoord3D(2) = j; elementCoord3d(3) = k


         !We want to cycle through each element of the array, and calculate 
         !each element's memory offset from the first element of the array.
         call ut_convertToMemoryOffset(dims, elementCoord3D, lbound(arr3D), ubound(arr3D), memoryOffset)
      
         if (memoryOffset /= arr3D(i,j,k)) then
            test3d = FAILURE_VALUE
            return
         end if

         !print *, "coordinates:", elementCoord3D, "memory offset=", memoryOffset

         !Now that we have the memory offset we would like to be able 
         !to restore the array coordinates.
         call ut_convertToArrayIndicies(dims, memoryOffset, lbound(arr3D), ubound(arr3D), newElementCoord3D)
         
         if ( (elementCoord3D(1) /= newElementCoord3D(1)) .or. &
              (elementCoord3D(2) /= newElementCoord3D(2)) .or. &
              (elementCoord3D(3) /= newElementCoord3D(3))) then
            test3d = FAILURE_VALUE
            return
         end if


      end do
   end do
end do

test3d = SUCCESS_VALUE

end function test3d
