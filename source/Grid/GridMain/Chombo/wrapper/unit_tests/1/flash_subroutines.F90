#include "Flash.h"

subroutine driver_init_flash()
  implicit none
  interface
     subroutine language_interface_init &
          () bind(c,name='language_interface_init')
     end subroutine language_interface_init
  end interface

  call language_interface_init
end subroutine driver_init_flash


subroutine driver_evolve_flash()
  implicit none
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer :: blockID, i, j, k, v, loopCounter
  logical :: success
  interface
     subroutine Grid_getBlkPtr(blockID, fortranPointer)
       implicit none
       integer, intent(in) :: blockID
       real, dimension(:,:,:,:), pointer :: fortranPointer
     end subroutine Grid_getBlkPtr
  end interface

  write(6,*) "We expect to pick up the complete set of integer "//&
       "values from 0 to 23 in contiguous memory locations"

  success = .true.
  loopCounter = 0
  blockID = 1;
  call Grid_getBlkPtr(blockID, dataPtr)
  
  !Use / modify the data in any way we wish.
  do k = lbound(dataPtr,4), ubound(dataPtr,4)
     do j = lbound(dataPtr,3), ubound(dataPtr,3)
        do i = lbound(dataPtr,2), ubound(dataPtr,2)
           do v = lbound(dataPtr,1), ubound(dataPtr,1)
              if (dataPtr(v,i,j,k) /= real(loopCounter)) then
                 success = .false.
              end if
              write(6,'(a,4i3,a,i6,a,f6.0)') " Fortran array element:", v,i,j,k, &
                   ".   Expected:", loopCounter, &
                   ", actual:", dataPtr(v,i,j,k)
              loopCounter = loopCounter + 1
           end do
        end do
     end do
  end do

  if (success .eqv. .true.) then
     write(6,*) "SUCCESS!  Picked up expected C++ assigned values"
  else
     write(6,*) "FAILURE!  Picked up unexpected C++ assigned values"
  end if

end subroutine driver_evolve_flash


subroutine driver_finalize_flash()
  implicit none
  interface
     subroutine language_interface_finalize &
          () bind(c,name='language_interface_finalize')
     end subroutine language_interface_finalize
  end interface

  call language_interface_finalize
end subroutine driver_finalize_flash


subroutine Grid_getBlkPtr(blockID, fortranPointer)
  use iso_c_binding
  implicit none
  integer, intent(in) :: blockID

  !This should really have type real(c_double), but since we know we are
  !promoting FLASH reals to double precision this is OK.
  real, dimension(:,:,:,:), pointer :: fortranPointer
  interface
     subroutine language_interface_get_blk_ptr &
          (blkID,dataPtr) bind(c,name='language_interface_get_blk_ptr')
       use iso_c_binding
       integer(c_int), value, intent(in) :: blkID
       type(c_ptr), intent(out) :: dataPtr
     end subroutine language_interface_get_blk_ptr
  end interface
  type(c_ptr) :: dataMemoryAddress
  integer(c_int) :: blkID


  nullify(fortranPointer)
  if (blockID == 1) then
     blkID = blockID

     !Obtain a block's raw memory address.
     call language_interface_get_blk_ptr(blkID, dataMemoryAddress);

     !Construct a Fortran pointer object using raw memory address and shape argument.
     call c_f_pointer(dataMemoryAddress, fortranPointer, (/NVAR,NXB,NYB,NZB/))
     if (.not.associated(fortranPointer)) then
        print *, "[ERROR!] Pointer not associated"
        stop
     end if
  end if
end subroutine Grid_getBlkPtr
