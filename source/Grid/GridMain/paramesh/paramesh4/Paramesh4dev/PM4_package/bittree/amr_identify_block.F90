!!****if* source/Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package/bittree/amr_identify_block
!!
!! NAME
!!
!!  amr_identify_block
!!
!! SYNOPSIS
!!
!!  call amr_identify_block(procs,lev,ijk,proc,locblk)
!!
!! DESCRIPTION
!!
!!  Determines the processor a block resides on from its spatial index.
!!  If a block matching the lev,ijk coordinates exists then its proc,locblk are
!!  filled in accordingly and lev,ijk are returned untouched.  If that block
!!  doesnt exist, then the finest parent containing that space will be returned
!!  in lev,ijk,proc,locblk. Hence block existence can be determined from change
!!  in the lev argument. If that space lies outside the domain then lev=0
!!  will be returned, boundary periodicity is not considered.
!!
!! ARGUMENTS
!!
!!  procs: (in) number of processors on mesh (result of MPI_COMM_SIZE)
!!  lev: (inout) 1-based level of block being hunted
!!  ijk: (inout) 0-based block coordinate (see gr_xyzToBlockLevel)
!!  proc: (out) 0-based processor rank block resides on
!!  locblk: (out) 1-based local block index for that processor
!!
!!***
subroutine amr_identify_block(procs, lev, ijk, proc, locblk)
  use bittree, only: bittree_block_count, bittree_morton
  use paramesh_dimensions, only: ndim
  use iso_c_binding, only: c_int
  implicit none
  
  integer, intent(in) :: procs
  integer, intent(inout) :: lev
  integer, intent(inout) :: ijk(ndim)
  integer, intent(out) :: proc
  integer, intent(out) :: locblk
  
  integer(c_int) :: levc, ijkc(ndim)
  integer :: blks, mort, procs_fat, blks_in_fat
  
  call bittree_block_count(blks)
  
  levc = lev - 1
  ijkc = ijk
  call bittree_morton(levc, ijkc, mort)
  lev = levc + 1
  ijk = ijkc
  
  procs_fat = mod(blks, procs)
  blks_in_fat = procs_fat*(blks/procs + 1)
  
  if(mort < blks_in_fat) then
    proc = mort/(blks/procs + 1)
    locblk = 1 + mod(mort, blks/procs + 1)
  else
    proc = procs_fat + (mort-blks_in_fat)/(blks/procs)
    locblk = 1 + mod(mort-blks_in_fat, blks/procs)
  end if
end subroutine
