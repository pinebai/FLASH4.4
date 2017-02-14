!!****if* source/Grid/localAPI/gr_hg_amr_1blk_bcset_work
!!
!! NAME
!!  gr_hg_amr_1blk_bcset_work
!!
!! SYNOPSIS
!! 
!!  call gr_hg_amr_1blk_bcset_work(integer(IN) :: lb, 
!!                           integer(IN) :: idest, 
!!                           integer(IN) :: idiag, 
!!                           integer(IN) :: idir)
!!
!! DESCRIPTION
!!
!!  Set boundary values of the exterior guardcells of a leaf block in the
!!  work array such that the given boundary value is satisfied.  The
!!  handled cases are MG_BND_DIRICHLET, MG_BND_NEUMANN, MG_BND_GIVENVAL
!! 
!!  Periodic boundaries are handled topologically and do not need a case
!!  here.  Allegedly this routine also handles cylindrical and spherical
!!  coordinates.
!!
!!  The implementation contains the same logic as the subroutine gr_hgSetExtBoundary,
!!  which is meant to work with PARAMESH2, except that here, boundary conditions
!!  are applied to only one block. This fits somewhat better into the PARAMESH3/4
!!  way of doing things.
!!
!!  The implementation is meant for comparisons of PARAMESH2 and PARAMESH3/4 results.
!!  In general, it would not be used when building Multigrid with PARAMESH 3 or 4.
!!
!! ARGUMENTS
!!
!!  lb     - block ID of the block to operate on
!!  idest  - selects the storage space in data_1blk.fh which is to
!!           be used in this call. If the leaf node is having its
!!           guardcells filled then this is set to 1, if its parent
!!           is being filled this is set to 2.
!!  idiag  - ignored variable for diagonal guardcell case
!!  idir   - sweep all directions (0), or x, y, or z only.
!!
!! RESULTS
!!
!!  Leaves some layers of guardcells in the work array for one block filled.
!!
!! NOTES
!!
!!  "Whoever fights monsters should see to it that in the process he does
!!   not become a monster.  And when you look long into an abyss, the
!!   abyss also looks into you" - Nietzsche
!!
!! SEE ALSO
!!
!!  gr_hgSetExtBoundary   (this is derived from it)
!!
!!***

subroutine gr_hg_amr_1blk_bcset_work (lb, idest, &
     idiag, idir)
  implicit none

  integer, intent(in) :: lb
  integer, intent(in) :: idest

  integer, intent(IN) :: idiag, idir

  return
end subroutine gr_hg_amr_1blk_bcset_work

