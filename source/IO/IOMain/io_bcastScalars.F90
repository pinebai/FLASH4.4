!!****if* source/IO/IOMain/io_bcastScalars
!!
!! NAME
!!  io_bcastScalars
!!
!! SYNOPSIS
!!
!!  io_bcastScalars()
!!
!! DESCRIPTION
!!
!! broadcasts scalars from global me to the other processors
!! calls nameValueLL_bcast to do the implementation
!! 
!! 
!! 
!!
!! ARGUMENTS
!!
!!        
!!
!!
!!
!!***

subroutine io_bcastScalars()

  use IO_data, only : io_scalar, io_globalMe
  
implicit none

  call nameValueLL_bcast(io_scalar, io_globalMe)

  return

end subroutine io_bcastScalars


  

