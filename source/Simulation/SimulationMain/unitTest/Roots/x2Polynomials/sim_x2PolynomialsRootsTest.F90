!!****if* source/Simulation/SimulationMain/unitTest/Roots/x2Polynomials/sim_x2PolynomialsRootsTest
!!
!! NAME
!!
!!  sim_x2PolynomialsRootsTest
!! 
!! SYNOPSIS
!!
!!  call sim_x2PolynomialsRootsTest (logical (out) :: perfect)
!!
!! DESCRIPTION
!!
!!  This routine tests the x2 polynomial root solver on a bunch of difficult to solve x2
!!  polynomials (widespread coefficients, near degenerate roots, etc...). The results will
!!  be printed to a special printout file. The obtained x2 polynomial roots will be compared
!!  to the exact analytical roots.
!!
!! ARGUMENTS 
!! 
!!  perfect : test status indicator (if true, test ran without errors).
!!
!!***

subroutine sim_x2PolynomialsRootsTest (perfect)

  use Roots_interface, ONLY: Roots_x2Polynomial

  use Simulation_data, ONLY: sim_maxRelativeAccuracy,   &
                             sim_numberOfx2Polynomials, &
                             sim_printUnit,             &
                             sim_rootsAnalytical,       &
                             sim_rootsRelativeAccuracy, &
                             sim_x2Polynomialx0Coeff,   &
                             sim_x2Polynomialx1Coeff,   &
                             sim_x2RelativeAccuracy

  implicit none

  logical, intent (out) :: perfect

  character (len=30) prF

  integer :: i, j, n
  integer :: nReal
  integer :: prU

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2

  real    :: A, B
  real    :: analyticRoot, obtainedRoot

  real    :: diffs (1:2,1:2)
  real    :: exact (1:2,1:2)
  real    :: roots (1:2,1:2)
!
!
!   ...Set printing controls.
!
!  
  prF = "(a,es25.16,a,es25.16,a)"      ! the printing format
  prU = sim_printUnit                  ! the printing file unit ID
!
!
!   ...Solve the set of x2 polynomials.
!
!  
  write (prU,prF) ' -----------------------------------------------------------------------'
  write (prU,prF) '            QUADRATIC POLYNOMIALS: x^2 + A x + B'
  write (prU,prF) ' -----------------------------------------------------------------------'

  do n = 1, sim_numberOfx2Polynomials

     A = sim_x2Polynomialx1Coeff (n)
     B = sim_x2Polynomialx0Coeff (n)

     write (prU,prF) ' '
     write (prU,prF) ' A = ',A
     write (prU,prF) ' B = ',B
     write (prU,prF) ' '

     call Roots_x2Polynomial (A, B, nReal, roots (1:2,1:2))
!
!
!   ...Determine all (relative) accuracies. In case an analytic root is equal to zero
!      use the absolute accuracy, i.e. the absolute value of the ontained computational
!      root.
!
!  
     do j = Re,Im
        do i = 1,2

           obtainedRoot = roots (i,j)
           analyticRoot = sim_rootsAnalytical (i,j,n)

           if (analyticRoot == 0.0) then
               sim_rootsRelativeAccuracy (i,j,n) = abs (obtainedRoot)
           else
               sim_rootsRelativeAccuracy (i,j,n) = abs ((analyticRoot - obtainedRoot) / analyticRoot)
           end if

        end do
     end do
!
!
!   ...Print out obtained roots and relative accuracies for current x3 polynomial.
!
!  
     exact (1:2,1:2) = sim_rootsAnalytical       (1:2,1:2,n)
     diffs (1:2,1:2) = sim_rootsRelativeAccuracy (1:2,1:2,n)

     write (prU,prF) ' '
     write (prU,prF) ' Exact  (analytical)  Root 1 (x + iy) = ',exact (1,Re),' + ',exact (1,Im),' i'
     write (prU,prF) ' Exact  (analytical)  Root 2 (x + iy) = ',exact (2,Re),' + ',exact (2,Im),' i'
     write (prU,prF) ' '
     write (prU,prF) ' x2 Polynomial Solver Root 1 (x + iy) = ',roots (1,Re),' + ',roots (1,Im),' i'
     write (prU,prF) ' x2 Polynomial Solver Root 2 (x + iy) = ',roots (2,Re),' + ',roots (2,Im),' i'
     write (prU,prF) ' '
     write (prU,prF) ' Relative Accuracy Root 1 (Re and Im) = ',diffs (1,Re),'   ',diffs (1,Im)
     write (prU,prF) ' Relative Accuracy Root 2 (Re and Im) = ',diffs (2,Re),'   ',diffs (2,Im)
     write (prU,prF) ' '
     write (prU,prF) ' -----------------------------------------------------------------------'

  end do
!
!
!   ...Determine overall relative accuracy (i.e. worst root) for each polynomial.
!
!  
  sim_x2RelativeAccuracy (:) = maxval (maxval (sim_rootsRelativeAccuracy, dim = 1), dim = 1)

  do n = 1, sim_numberOfx2Polynomials
     write (prU,prF) ' Relative Accuracy of x2 Polynomial = ', sim_x2RelativeAccuracy (n)
  end do
!
!
!   ...Compare the x2 polynomial relative accuracies obtained with the corresponding
!      upper bounds and decide, if the test passed.
!
!  
  perfect = all (sim_x2RelativeAccuracy (:) <= sim_maxRelativeAccuracy (:))
!
!
!   ...Ready!
!
!
  return
end subroutine sim_x2PolynomialsRootsTest
