!!****if* source/IO/IOMain/io_checkBlockShape
!!
!! NAME
!!
!!  io_checkBlockShape
!!
!! SYNOPSIS
!!
!!  call io_checkBlockShape(integer,intent(IN)  :: numblocks)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   numblocks : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine io_checkBlockShape(numBlocks)
  use IO_interface, ONLY: IO_getPrevScalar
  use Driver_interface, ONLY: Driver_abortFlash

  implicit none

  integer,intent(IN) :: numBlocks

  integer :: fileNxb, fileNyb, fileNzb

  if (numBlocks == 0) return

  call IO_getPrevScalar("nxb", fileNxb)
  call IO_getPrevScalar("nyb", fileNyb)
  call IO_getPrevScalar("nzb", fileNzb)

  if (fileNxb .NE. NXB) call Driver_abortFlash("NXB in checkpoint does not match this executable!")
  if (fileNyb .NE. NYB) call Driver_abortFlash("NYB in checkpoint does not match this executable!")
  if (fileNzb .NE. NZB) call Driver_abortFlash("NZB in checkpoint does not match this executable!")

end subroutine io_checkBlockShape
