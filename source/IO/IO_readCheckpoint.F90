!!****f* source/IO/IO_readCheckpoint
!!
!! NAME
!!
!!  IO_readCheckpoint
!!
!!
!! SYNOPSIS
!!
!!  IO_readCheckpoint()
!!
!!
!! DESCRIPTION
!!
!!  IO_readCheckpoint is a generic subroutine that retrieves
!!  the unklabels and then calls io_readData
!!  which is specific to pnetcdf, hdf5 and
!!  the necessary grid, UG, paramesh etc.
!!  io_readData reads a checkpoint file and reinitializes
!!  grid and scalar values to resume the run
!!  
!!
!!
!! ARGUMENTS
!!  
!!
!!
!!***

subroutine IO_readCheckpoint()


implicit none
 
  return
end subroutine IO_readCheckpoint

