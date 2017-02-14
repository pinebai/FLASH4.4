!!****if* source/physics/Hydro/localAPI/hy_uhd_biermannSource
!!
!! NAME
!!
!!  hy_uhd_biermannSource
!!
!! SYNOPSIS
!!
!!  call hy_uhd_biermannSource( integer (IN) :: blockCount,
!!                         integer (IN) :: blockList(blockCount),
!!                         real    (IN) :: dt )
!!
!! DESCRIPTION
!!
!! Implement Biermann Battery Term as a source to the magnetic field. This is a stub.
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList  -  array holding local IDs of blocks on which to advance
!!  dt         -  timestep
!!***

Subroutine hy_uhd_biermannSource ( blockCount, blockList, dt )
  
  implicit none

  ! Arguments:
  integer, intent(IN) :: blockCount
  integer, intent(IN) :: blockList(blockCount)
  real,    intent(IN) :: dt

  return

End Subroutine hy_uhd_biermannSource
