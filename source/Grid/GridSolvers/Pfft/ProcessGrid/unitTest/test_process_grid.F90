#include "constants.h"

!This is a standalone program that prints the selected
!JAXIS and KAXIS process grid for different values of total
!processes.
!
!I'm using it to check that there are no regressions
!in gr_pfftMakePencilIn3dSpace after my optimizations:
!test_process_grid > file1   ... (before optimization)
!test_process_grid > file2   ... (after optimization)
!diff file1 file2

program test_process_grid
  use gr_pfftInterfaceTypeDecl, ONLY : PossibleGrid_t
  use gr_pfftInterface, ONLY : gr_pfftMakePencilIn3dSpace, &
       gr_pfftFnArgEasyConstraint
  implicit none
  type(PossibleGrid_t) :: pencilGrid
  integer :: pfftNumProcs

  !This is a required argument for gr_pfftMakePencilIn3dSpace but not
  !for the easy constraint function named gr_pfftFnArgEasyConstraint.
  integer, dimension(MDIM) :: pfftGlobalLen
  pfftGlobalLen = 1  

  do pfftNumProcs = 1, 10000
     pencilGrid = gr_pfftMakePencilIn3dSpace &
          (pfftGlobalLen, pfftNumProcs, gr_pfftFnArgEasyConstraint)
     print *, pfftNumProcs, pencilGrid % jprocs, pencilGrid % kprocs
  end do
end program test_process_grid
