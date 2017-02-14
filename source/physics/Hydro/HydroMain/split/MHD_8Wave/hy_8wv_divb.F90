!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/hy_8wv_divb
!!
!! NAME
!!
!!  hy_8wv_divb
!!
!!
!! SYNOPSIS
!!
!!  hy_8wv_divb(
!!              integer(IN) :: blockCount,
!!              integer(IN) :: blockList(blockCount),
!!              real(IN)    :: dt)
!!
!!
!! DESCRIPTION
!!
!!  Cleans divb using the parabolic diffusion method of
!!  Marder, J. Comput. Phys., 68, 48, 1987.
!!  Stub function for diffusive and projection divb cleaning methods,
!!  whose purpose is to reduce or remove magnetic monopole contamination.
!!
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  dt         - time step
!!
!!***

subroutine hy_8wv_divb(blockCount,blockList,dt)
  implicit none
  integer, intent(IN) :: blockCount
  integer, intent(IN), dimension(blockCount) :: blockList
  real, intent(IN)    :: dt

end subroutine hy_8wv_divb
