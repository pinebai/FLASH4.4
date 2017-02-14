!!****if* source/Grid/GridSolvers/Pfft/ProcessGrid/gr_pfftMakePencilIn3dSpace
!!
!! NAME
!!
!!  gr_pfftMakePencilIn3dSpace
!!
!! SYNOPSIS
!!  
!!  gr_pfftMakePencilIn3dSpace(integer(IN) :: pencilGlobalLen(MDIM), &
!!                             integer(IN) :: totalProcs, &
!!                             logical function :: gr_pfftFnArgConstraint)
!!
!! DESCRIPTION
!!
!!  Generates our PFFT grid.  The grid is defined such that the 
!!  number of processors in the IAXIS is 1.  Additional constraints 
!!  are specified in a passed function.  gr_pfftMakePencilIn3dSpace 
!!  returns a datatype containing the proposed processor grid.
!!
!! ARGUMENTS
!!
!!  pencilGlobalLen - An array containing the number of cells in the 
!!                    global pencil grid
!!  totalProcs - the number of processors in the run
!!  gr_pfftFnArgConstraint - A function specifying the constraints 
!!                           that the generated grid must satisfy.
!!
!!***
function gr_pfftMakePencilIn3dSpace(pencilGlobalLen, totalProcs, &
     gr_pfftFnArgConstraint)
#include "constants.h"
  use gr_pfftInterfaceTypeDecl, ONLY : PossibleGrid_t
  implicit none
  integer, dimension(MDIM), intent(IN) :: pencilGlobalLen
  integer, intent(IN) :: totalProcs
  interface     
     logical function gr_pfftFnArgConstraint &
          (pencilGlobalLen, totalProcs, iProcs, jProcs, kProcs)
       implicit none
       integer, dimension(MDIM), intent(IN) :: pencilGlobalLen
       integer, intent(IN) :: totalProcs, iProcs, jProcs, kProcs
     end function gr_pfftFnArgConstraint
  end interface
  type (PossibleGrid_t) :: gr_pfftMakePencilIn3dSpace  !Return.
  type (PossibleGrid_t), dimension(2) :: possibleGrid
  integer :: n, m, jProcs, kProcs, sqrtTotalProcs

  !Ensure we return NONEXISTENT if we find no pencil grid.
  gr_pfftMakePencilIn3dSpace % jProcs = NONEXISTENT
  gr_pfftMakePencilIn3dSpace % kProcs = NONEXISTENT

  !The best (JAXIS,KAXIS) process grid is a square:
  !jprocs = kprocs = sqrt(totalProcs).
  !If we have 12 processors then the square root cast to integer is 3,
  !and so we only test (1,12), (2,6), (3,4).  We do not need to iterate
  !beyond 3 because the remaining grid combinations are just be
  !the integer values swapped around, i.e. (12,1), (6,2), (4,3).
  sqrtTotalProcs = int(sqrt(real(totalProcs)))
  do n = sqrtTotalProcs, 1, -1

     possibleGrid(1) % jProcs = n
     possibleGrid(1) % kProcs = totalProcs / n

     !The process grid must equal totalProcs
     if (possibleGrid(1) % jProcs * possibleGrid(1) % kProcs &
          == totalProcs) then

        !Some versions of the constraint function also test that
        !the global grid points are exactly divisible by the process
        !grid.  This may mean that the best grid has jprocs > kprocs, 
        !e.g. (4,3) from our earlier example.  This would not
        !be tested in the n loop because the highest value is
        !sqrtTotalProcs.  Therefore we include an m loop to test
        !the swapped around jprocs and kprocs values.
        possibleGrid(2) % jProcs = possibleGrid(1) % kProcs
        possibleGrid(2) % kProcs = possibleGrid(1) % jProcs

        do m = 1, 2
           jProcs = possibleGrid(m) % jProcs
           kProcs = possibleGrid(m) % kProcs

           if (gr_pfftFnArgConstraint(pencilGlobalLen, totalProcs, &
                1, jProcs, kProcs) .eqv. .true.) then
              !We accept the first match because the n index begins at
              !sqrtTotalProcs and so the first match is the most square.
              gr_pfftMakePencilIn3dSpace % jProcs = jProcs
              gr_pfftMakePencilIn3dSpace % kProcs = kProcs
              return
           end if
        end do
     end if
  end do
end function gr_pfftMakePencilIn3dSpace
