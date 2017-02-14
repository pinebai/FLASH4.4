!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/UNIT_TEST
!!
!! NAME
!!  UNIT_TEST
!!
!! SYNOPSIS
!!
!!   use global_Data_and_stubs
!! 
!! DESCRIPTION
!!
!! NOTE: This is a standalone test for the subroutine:
!! gr_ptFindValidChildren which is a helper routine in gr_ptFindNegh.F90. 
!! We figure out which child blocks are next to each guard cell region 
!! of a source block.  expectedChild is an array of which child blocks 
!! should be neighbors to each guard cell region, and validChild is
!! an array of our predicted neighbors.  The tests are successful 
!! if validChild is identical to expectedChild
!!
!!***

module global_data_and_subs

  implicit none
  integer, parameter :: MDIM = 3, NONEXISTENT = -1
  integer, parameter :: LEFT_EDGE = 1, CENTER = 2, RIGHT_EDGE = 3
  integer, parameter :: LOW = 1, HIGH = 2
  integer, parameter :: TOTAL_CHILDREN = 2**MDIM
  integer, parameter, dimension(MDIM,TOTAL_CHILDREN) :: blockChild = &
       reshape (source = (/&
       LOW,  LOW,  LOW , &
       HIGH, LOW,  LOW , &
       LOW,  HIGH, LOW , &
       HIGH, HIGH, LOW , &
       LOW,  LOW,  HIGH, &
       HIGH, LOW,  HIGH, &
       LOW,  HIGH, HIGH, &
       HIGH, HIGH, HIGH  &
       /), shape = (/MDIM,TOTAL_CHILDREN/))
  integer, dimension(MDIM) :: guardCellID
  integer :: NDIM
  logical :: m_success

contains

  !All state is passed to this subroutine so that it can be
  !tested independently of FLASH.
  subroutine gr_ptFindValidChildren(maxDims, maxChildren, blockChild, &
       guardCellID, dims, actualChildren, validChild)

    implicit none
    integer, intent(IN) :: maxDims, maxChildren
    integer, dimension(maxDims,maxChildren), intent(IN) :: blockChild
    integer, dimension(maxDims), intent(IN) :: guardCellID
    integer, intent(IN) :: dims, actualChildren
    logical, dimension(actualChildren), intent(OUT) :: validChild

    integer :: eachChild, eachAxis
    integer :: numNegh, validChildren

    validChild(:) = .true.
    do eachChild = 1, actualChildren

       do eachAxis = 1, dims

          !We can only eliminate children when they are not in the center.
          if ( (guardCellID(eachAxis) == LEFT_EDGE .and. & 
               blockChild(eachAxis,eachChild) == LOW) .or. &
               (guardCellID(eachAxis) == RIGHT_EDGE .and. &
               blockChild(eachAxis,eachChild) == HIGH) ) then
             validChild(eachChild) = .false.
          end if

       end do
    end do


    !Add an extra check that the number of child blocks is correct 
    !(this code is not in the actual subroutine).
    numNegh = 2 ** (count(guardCellID(1:dims) == CENTER))
    validChildren = count(validChild(1:actualChildren) .eqv. .true.)


    if (validChildren /= numNegh) then
       print *, "Actual neighbors:", numNegh, "Predicted neighbors:", validChildren
       stop
    end if

  end subroutine gr_ptFindValidChildren


  subroutine test_gr_findWhichChildren(guardCellID, dims, validChild)

    implicit none
    integer, dimension(MDIM), intent(IN) :: guardCellID
    integer, intent(IN) :: dims
    logical, dimension(2**dims), intent(OUT) :: validChild

    integer, allocatable, dimension(:) :: whichChildren
    integer :: i, numNegh

    validChild(:) = .false.
    numNegh = 2 ** (count(guardCellID(1:dims) == CENTER))

    allocate(whichChildren(numNegh))
    call gr_findWhichChildren(numNegh, guardCellID, whichChildren)
    do i = 1, numNegh
       validChild(whichChildren(i)) = .true.
    end do
    deallocate(whichChildren)

  end subroutine test_gr_findWhichChildren


  subroutine gr_findWhichChildren(numNegh,Negh,whichChildren)

    implicit none

    integer,intent(IN) :: numNegh
    integer, dimension(MDIM),intent(IN) :: Negh
    integer, intent(OUT) :: whichChildren(numNegh)
    integer,dimension(NDIM,LOW:HIGH) :: pos
    integer :: i,j,k,n,m

    pos=0
    do i = 1, NDIM

       !! If the negh is on left edge, 
       !!that is the right half of the neighbor's parent block
       if(Negh(i)==LEFT_EDGE)pos(i,1)=1

       !! If the negh is on right edge, 
       !! that is the left half of the neighbor's parent block
       if(Negh(i)==RIGHT_EDGE)pos(i,1)=-1

       !! if center, we are interested in both halves.
       !! n here is keeping track of count of the children.
       if(Negh(i)==CENTER)then
          pos(i,LOW)=-1; pos(i,HIGH)=1
       end if
    end do

    !! Now translate those positions into the child places within the parent
    !! block
    whichChildren(:)=1
    k=1
    m=2
    do i = 1,NDIM
       if(pos(i,LOW)>=0)whichChildren(:)=whichChildren(:)+k
       if(pos(i,HIGH)>0) then
          if(m==2) then
             whichChildren(m:numNegh)=whichChildren(m:numNegh)+k
          else
             whichChildren(m:numNegh)=whichChildren(LOW:HIGH)+k
          end if
          m=m+1
       end if
       k=k*2
    end do
    return
  end subroutine gr_findWhichChildren



  subroutine CompareResults(testStr, dims, validChild, expectedChild)
    implicit none
    character (len=3), intent(IN) :: testStr 
    integer, intent(IN) :: dims
    logical, dimension(1:2**dims), intent(IN) :: validChild, expectedChild

    if (any(validChild(:) .neqv. expectedChild(:))) then 
       print *, "Failure in test:"//testStr
       m_success = .false.  !module level variable
    end if
  end subroutine CompareResults


end module global_data_and_subs


!This is the actual program.
program UNIT_TEST

  use global_data_and_subs
  implicit none
  integer :: usableChildren
  logical, allocatable, dimension(:) :: validChild, expectedChild

  m_success = .true.  !module level variable.

  !TOTAL TESTS:
  !1D (2 corner) = 2
  !2D (4 corner & 4 face) = 8
  !3D (8 corner & 12 edge & 6 face) = 26
  !NUMBER OF TESTS = 36


  !-----------------------------------------------------------------
  ! 1D test cases  (2 corner orientations)
  !-----------------------------------------------------------------
  NDIM = 1 !Module level variable.
  usableChildren = 2**NDIM
  allocate(validChild(usableChildren), expectedChild(usableChildren))


  guardCellID = (/LEFT_EDGE, NONEXISTENT, NONEXISTENT/)
  expectedChild = (/.false., .true./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('1a ', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('1b ', NDIM, validChild, expectedChild)



  guardCellID = (/RIGHT_EDGE, NONEXISTENT, NONEXISTENT/)
  expectedChild = (/.true., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('2a ', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('2b ', NDIM, validChild, expectedChild)



  deallocate(validChild, expectedChild)


  !-----------------------------------------------------------------
  ! 2D test cases  (4 possible corner guard cell orientations)
  !-----------------------------------------------------------------
  NDIM = 2 !Module level variable.
  usableChildren = 2**NDIM
  allocate(validChild(usableChildren), expectedChild(usableChildren))

  guardCellID = (/LEFT_EDGE, LEFT_EDGE, NONEXISTENT/)
  expectedChild = (/.false., .false., .false., .true./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('3a ', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('3b ', NDIM, validChild, expectedChild)




  guardCellID = (/RIGHT_EDGE, LEFT_EDGE, NONEXISTENT/)
  expectedChild = (/.false., .false., .true., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('4a ', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('4b ', NDIM, validChild, expectedChild)



  guardCellID = (/LEFT_EDGE, RIGHT_EDGE, NONEXISTENT/)
  expectedChild = (/.false., .true., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('5a ', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('5b ', NDIM, validChild, expectedChild)




  guardCellID = (/RIGHT_EDGE, RIGHT_EDGE, NONEXISTENT/)
  expectedChild = (/.true., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('6a ', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('6b ', NDIM, validChild, expectedChild)



  !-----------------------------------------------------------------
  ! 2D test cases  (4 possible edge guard cell orientations)
  !-----------------------------------------------------------------
  guardCellID = (/LEFT_EDGE, CENTER, NONEXISTENT/)
  expectedChild = (/.false., .true., .false., .true./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('7a ', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('7b ', NDIM, validChild, expectedChild)


  guardCellID = (/RIGHT_EDGE, CENTER, NONEXISTENT/)
  expectedChild = (/.true., .false., .true., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('8a ', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('8b ', NDIM, validChild, expectedChild)



  guardCellID = (/CENTER, LEFT_EDGE, NONEXISTENT/)
  expectedChild = (/.false., .false., .true., .true./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('9a ', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('9b ', NDIM, validChild, expectedChild)



  guardCellID = (/CENTER, RIGHT_EDGE, NONEXISTENT/)
  expectedChild = (/.true., .true., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('10a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('10b', NDIM, validChild, expectedChild)


  deallocate(validChild, expectedChild)


  !-----------------------------------------------------------------
  ! 3D test cases  (8 possible corner guard cell orientations)
  !-----------------------------------------------------------------
  NDIM = 3 !Module level variable.
  usableChildren = 2**NDIM
  allocate(validChild(usableChildren), expectedChild(usableChildren))

  guardCellID = (/LEFT_EDGE, LEFT_EDGE, LEFT_EDGE/)
  expectedChild = (/.false., .false., .false., .false., .false., .false., .false., .true./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('11a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('11b', NDIM, validChild, expectedChild)



  guardCellID = (/RIGHT_EDGE, LEFT_EDGE, LEFT_EDGE/)
  expectedChild = (/.false., .false., .false., .false., .false., .false., .true., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('12a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('12b', NDIM, validChild, expectedChild)



  guardCellID = (/LEFT_EDGE, RIGHT_EDGE, LEFT_EDGE/)
  expectedChild = (/.false., .false., .false., .false., .false., .true., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('13a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('13b', NDIM, validChild, expectedChild)



  guardCellID = (/RIGHT_EDGE, RIGHT_EDGE, LEFT_EDGE/)
  expectedChild = (/.false., .false., .false., .false., .true., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('14a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('14b', NDIM, validChild, expectedChild)



  guardCellID = (/LEFT_EDGE, LEFT_EDGE, RIGHT_EDGE/)
  expectedChild = (/.false., .false., .false., .true., .false., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('15a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('15b', NDIM, validChild, expectedChild)



  guardCellID = (/RIGHT_EDGE, LEFT_EDGE, RIGHT_EDGE/)
  expectedChild = (/.false., .false., .true., .false., .false., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('16a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('16b', NDIM, validChild, expectedChild)



  guardCellID = (/LEFT_EDGE, RIGHT_EDGE, RIGHT_EDGE/)
  expectedChild = (/.false., .true., .false., .false., .false., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('17a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('17b', NDIM, validChild, expectedChild)



  guardCellID = (/RIGHT_EDGE, RIGHT_EDGE, RIGHT_EDGE/)
  expectedChild = (/.true., .false., .false., .false., .false., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('18a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('18b', NDIM, validChild, expectedChild)



  !-----------------------------------------------------------------
  ! 3D test cases  (12 possible edge guard cell orientations)
  !-----------------------------------------------------------------
  guardCellID = (/CENTER, LEFT_EDGE, LEFT_EDGE/)
  expectedChild = (/.false., .false., .false., .false., .false., .false., .true., .true./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('19a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('19b', NDIM, validChild, expectedChild)



  guardCellID = (/CENTER, RIGHT_EDGE, LEFT_EDGE/)
  expectedChild = (/.false., .false., .false., .false., .true., .true., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('20a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('20b', NDIM, validChild, expectedChild)



  guardCellID = (/LEFT_EDGE, CENTER, LEFT_EDGE/)
  expectedChild = (/.false., .false., .false., .false., .false., .true., .false., .true./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('21a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('21b', NDIM, validChild, expectedChild)


  guardCellID = (/RIGHT_EDGE, CENTER, LEFT_EDGE/)
  expectedChild = (/.false., .false., .false., .false., .true., .false., .true., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('22a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('22b', NDIM, validChild, expectedChild)



  guardCellID = (/LEFT_EDGE, LEFT_EDGE, CENTER/)
  expectedChild = (/.false., .false., .false., .true., .false., .false., .false., .true./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('23a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('23b', NDIM, validChild, expectedChild)



  guardCellID = (/RIGHT_EDGE, LEFT_EDGE, CENTER/)
  expectedChild = (/.false., .false., .true., .false., .false., .false., .true., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('24a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('24b', NDIM, validChild, expectedChild)



  guardCellID = (/LEFT_EDGE, RIGHT_EDGE, CENTER/)
  expectedChild = (/.false., .true., .false., .false., .false., .true., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('25a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('25b', NDIM, validChild, expectedChild)



  guardCellID = (/RIGHT_EDGE, RIGHT_EDGE, CENTER/)
  expectedChild = (/.true., .false., .false., .false., .true., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('26a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('26b', NDIM, validChild, expectedChild)



  guardCellID = (/CENTER, LEFT_EDGE, RIGHT_EDGE/)
  expectedChild = (/.false., .false., .true., .true., .false., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('27a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('27b', NDIM, validChild, expectedChild)



  guardCellID = (/CENTER, RIGHT_EDGE, RIGHT_EDGE/)
  expectedChild = (/.true., .true., .false., .false., .false., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('28a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('28b', NDIM, validChild, expectedChild)



  guardCellID = (/LEFT_EDGE, CENTER, RIGHT_EDGE/)
  expectedChild = (/.false., .true., .false., .true., .false., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('29a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('29b', NDIM, validChild, expectedChild)



  guardCellID = (/RIGHT_EDGE, CENTER, RIGHT_EDGE/)
  expectedChild = (/.true., .false., .true., .false., .false., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('30a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('30b', NDIM, validChild, expectedChild)



  !-----------------------------------------------------------------
  ! 3D test cases  (6 possible face guard cell orientations)
  !-----------------------------------------------------------------
  guardCellID = (/CENTER, CENTER, LEFT_EDGE/)
  expectedChild = (/.false., .false., .false., .false., .true., .true., .true., .true./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('31a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('31b', NDIM, validChild, expectedChild)



  guardCellID = (/CENTER, CENTER, RIGHT_EDGE/)
  expectedChild = (/.true., .true., .true., .true., .false., .false., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('32a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('32b', NDIM, validChild, expectedChild)



  guardCellID = (/CENTER, LEFT_EDGE, CENTER/)
  expectedChild = (/.false., .false., .true., .true., .false., .false., .true., .true./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('33a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('33b', NDIM, validChild, expectedChild)



  guardCellID = (/CENTER, RIGHT_EDGE, CENTER/)
  expectedChild = (/.true., .true., .false., .false., .true., .true., .false., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('34a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('34b', NDIM, validChild, expectedChild)



  guardCellID = (/LEFT_EDGE, CENTER, CENTER/)
  expectedChild = (/.false., .true., .false., .true., .false., .true., .false., .true./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('35a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('35b', NDIM, validChild, expectedChild)



  guardCellID = (/RIGHT_EDGE, CENTER, CENTER/)
  expectedChild = (/.true., .false., .true., .false., .true., .false., .true., .false./)
  call gr_ptFindValidChildren(MDIM, TOTAL_CHILDREN, blockChild, &
       guardCellID, NDIM, usableChildren, validChild)
  call CompareResults('36a', NDIM, validChild, expectedChild)

  call test_gr_findWhichChildren(guardCellID, NDIM, validChild)
  call CompareResults('36b', NDIM, validChild, expectedChild)


  deallocate(validChild, expectedChild)

  if (m_success .eqv. .true.) then
     print *, "Test complete: All 36 tests passed!!!"
  else
     print *, "Test complete: WITH FAILURES :("
  end if

end program UNIT_TEST
