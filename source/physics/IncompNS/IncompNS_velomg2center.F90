!!****f* source/physics/IncompNS/IncompNS_velomg2center
!!
!! NAME
!!
!!  IncompNS_velomg2center
!!
!! SYNOPSIS
!!
!!  call IncompNS_velomg2center(integer(in) :: blockList,
!!                              integer(in) :: blockCount)
!!
!! DESCRIPTION
!!
!!   Compute cell-centered (or volume-averaged) versions of some
!!   quantities, by interpolation or averaging from face-centered
!!   versions.
!!
!! ARGUMENTS
!!
!!   blockList : list of block IDs
!!
!!   blockCount : number of block IDs in the list
!!
!! NOTES
!!
!!  To write cell centered velocities and vorticity velx,vely,velz,omgx,omgy,omgz,
!!  add REQUIRES physics/IncompNS/IncompNSExtras to Config.
!!
!!***

subroutine IncompNS_velomg2center(blockList,blockCount)

  implicit none
  integer,intent(in) :: blockCount
  integer,intent(in) :: blockList(MAXBLOCKS)

end subroutine IncompNS_velomg2center
