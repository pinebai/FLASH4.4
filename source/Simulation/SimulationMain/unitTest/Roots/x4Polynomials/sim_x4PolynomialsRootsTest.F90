!!****if* source/Simulation/SimulationMain/unitTest/Roots/x4Polynomials/sim_x4PolynomialsRootsTest
!!
!! NAME
!!
!!  sim_x4PolynomialsRootsTest
!! 
!! SYNOPSIS
!!
!!  call sim_x4PolynomialsRootsTest (logical (out) :: perfect)
!!
!! DESCRIPTION
!!
!!  This routine tests the x4 polynomial root solver on a bunch of difficult to solve x4
!!  polynomials (widespread coefficients, near degenerate roots, etc...). The results will
!!  be printed to a special printout file. The obtained x4 polynomial roots will be compared
!!  to the exact analytical roots.
!!
!! ARGUMENTS 
!! 
!!  perfect : test status indicator (if true, test ran without errors).
!!
!!***

subroutine sim_x4PolynomialsRootsTest (perfect)

  use Roots_interface, ONLY: Roots_x4Polynomial

  use Simulation_data, ONLY: sim_maxRelativeAccuracy,   &
                             sim_numberOfx4Polynomials, &
                             sim_printInfo,             &
                             sim_printUnit,             &
                             sim_rootsAnalytical,       &
                             sim_rootsRelativeAccuracy, &
                             sim_x4Polynomialx0Coeff,   &
                             sim_x4Polynomialx1Coeff,   &
                             sim_x4Polynomialx2Coeff,   &
                             sim_x4Polynomialx3Coeff,   &
                             sim_x4RelativeAccuracy

  implicit none

  logical, intent (out) :: perfect

  character (len=30) prF

  integer :: i, j, n
  integer :: nReal
  integer :: prU

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2

  real    :: A, B, C, D
  real    :: analyticRoot, obtainedRoot

  real    :: diffs (1:4,1:2)
  real    :: exact (1:4,1:2)
  real    :: roots (1:4,1:2)
!
!
!   ...Set printing controls.
!
!  
  prF = "(a,es25.16,a,es25.16,a)"      ! the printing format
  prU = sim_printUnit                  ! the printing file unit ID
!
!
!   ...Solve the set of x4 polynomials.
!
!  
  write (prU,prF) ' -----------------------------------------------------------------------'
  write (prU,prF) '            QUARTIC POLYNOMIALS: x^4 + A x^3 + B x^2 + C x + D'
  write (prU,prF) ' -----------------------------------------------------------------------'

  do n = 1, sim_numberOfx4Polynomials

     A = sim_x4Polynomialx3Coeff (n)
     B = sim_x4Polynomialx2Coeff (n)
     C = sim_x4Polynomialx1Coeff (n)
     D = sim_x4Polynomialx0Coeff (n)

     write (prU,prF) ' '
     write (prU,prF) ' A = ',A
     write (prU,prF) ' B = ',B
     write (prU,prF) ' C = ',C
     write (prU,prF) ' D = ',D
     write (prU,prF) ' '

     call Roots_x4Polynomial (A, B, C, D, nReal, roots (1:4,1:2), sim_printInfo, sim_printUnit)
!
!
!   ...Determine all (relative) accuracies. In case an analytic root is equal to zero
!      use the absolute accuracy, i.e. the absolute value of the ontained computational
!      root.
!
!  
     do j = Re,Im
        do i = 1,4

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
!   ...Print out obtained roots and relative accuracies for current x4 polynomial.
!
!  
     exact (1:4,1:2) = sim_rootsAnalytical       (1:4,1:2,n)
     diffs (1:4,1:2) = sim_rootsRelativeAccuracy (1:4,1:2,n)

     write (prU,prF) ' '
     write (prU,prF) ' Exact  (analytical)  Root 1 (x + iy) = ',exact (1,Re),' + ',exact (1,Im),' i'
     write (prU,prF) ' Exact  (analytical)  Root 2 (x + iy) = ',exact (2,Re),' + ',exact (2,Im),' i'
     write (prU,prF) ' Exact  (analytical)  Root 3 (x + iy) = ',exact (3,Re),' + ',exact (3,Im),' i'
     write (prU,prF) ' Exact  (analytical)  Root 4 (x + iy) = ',exact (4,Re),' + ',exact (4,Im),' i'
     write (prU,prF) ' '
     write (prU,prF) ' x4 Polynomial Solver Root 1 (x + iy) = ',roots (1,Re),' + ',roots (1,Im),' i'
     write (prU,prF) ' x4 Polynomial Solver Root 2 (x + iy) = ',roots (2,Re),' + ',roots (2,Im),' i'
     write (prU,prF) ' x4 Polynomial Solver Root 3 (x + iy) = ',roots (3,Re),' + ',roots (3,Im),' i'
     write (prU,prF) ' x4 Polynomial Solver Root 4 (x + iy) = ',roots (4,Re),' + ',roots (4,Im),' i'
     write (prU,prF) ' '
     write (prU,prF) ' Relative Accuracy Root 1 (Re and Im) = ',diffs (1,Re),'   ',diffs (1,Im)
     write (prU,prF) ' Relative Accuracy Root 2 (Re and Im) = ',diffs (2,Re),'   ',diffs (2,Im)
     write (prU,prF) ' Relative Accuracy Root 3 (Re and Im) = ',diffs (3,Re),'   ',diffs (3,Im)
     write (prU,prF) ' Relative Accuracy Root 4 (Re and Im) = ',diffs (4,Re),'   ',diffs (4,Im)
     write (prU,prF) ' '
     write (prU,prF) ' -----------------------------------------------------------------------'

  end do
!
!
!   ...Determine overall relative accuracy (i.e. worst root) for each polynomial.
!
!  
  sim_x4RelativeAccuracy (:) = maxval (maxval (sim_rootsRelativeAccuracy, dim = 1), dim = 1)

  do n = 1, sim_numberOfx4Polynomials
     write (prU,prF) ' Relative Accuracy of x4 Polynomial = ', sim_x4RelativeAccuracy (n)
  end do
!
!
!   ...Compare the x4 polynomial relative accuracies obtained with the corresponding
!      upper bounds and decide, if the test passed.
!
!  
  perfect = all (sim_x4RelativeAccuracy (:) <= sim_maxRelativeAccuracy (:))
!
!
!   ...Ready!
!
!
  return
end subroutine sim_x4PolynomialsRootsTest
