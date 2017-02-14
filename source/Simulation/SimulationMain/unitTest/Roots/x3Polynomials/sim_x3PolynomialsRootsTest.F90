!!****if* source/Simulation/SimulationMain/unitTest/Roots/x3Polynomials/sim_x3PolynomialsRootsTest
!!
!! NAME
!!
!!  sim_x3PolynomialsRootsTest
!! 
!! SYNOPSIS
!!
!!  call sim_x3PolynomialsRootsTest (logical (out) :: perfect)
!!
!! DESCRIPTION
!!
!!  This routine tests the x3 polynomial root solver on a bunch of difficult to solve x3
!!  polynomials (widespread coefficients, near degenerate roots, etc...). The results will
!!  be printed to a special printout file. The obtained x3 polynomial roots will be compared
!!  to the exact analytical roots.
!!
!! ARGUMENTS 
!! 
!!  perfect : test status indicator (if true, test ran without errors).
!!
!!***

subroutine sim_x3PolynomialsRootsTest (perfect)

  use Roots_interface, ONLY: Roots_x3Polynomial

  use Simulation_data, ONLY: sim_maxRelativeAccuracy,   &
                             sim_numberOfx3Polynomials, &
                             sim_printInfo,             &
                             sim_printUnit,             &
                             sim_rootsAnalytical,       &
                             sim_rootsRelativeAccuracy, &
                             sim_x3Polynomialx0Coeff,   &
                             sim_x3Polynomialx1Coeff,   &
                             sim_x3Polynomialx2Coeff,   &
                             sim_x3RelativeAccuracy

  implicit none

  logical, intent (out) :: perfect

  character (len=30) prF

  integer :: i, j, n
  integer :: nReal
  integer :: prU

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2

  real    :: A, B, C
  real    :: analyticRoot, obtainedRoot

  real    :: diffs (1:3,1:2)
  real    :: exact (1:3,1:2)
  real    :: roots (1:3,1:2)
!
!
!   ...Set printing controls.
!
!  
  prF = "(a,es25.16,a,es25.16,a)"      ! the printing format
  prU = sim_printUnit                  ! the printing file unit ID
!
!
!   ...Solve the set of x3 polynomials.
!
!  
  write (prU,prF) ' -----------------------------------------------------------------------'
  write (prU,prF) '              CUBIC POLYNOMIALS: x^3 + A x^2 + B x + C'
  write (prU,prF) ' -----------------------------------------------------------------------'

  do n = 1, sim_numberOfx3Polynomials

     A = sim_x3Polynomialx2Coeff (n)
     B = sim_x3Polynomialx1Coeff (n)
     C = sim_x3Polynomialx0Coeff (n)

     write (prU,prF) ' '
     write (prU,prF) ' A = ',A
     write (prU,prF) ' B = ',B
     write (prU,prF) ' C = ',C
     write (prU,prF) ' '

     call Roots_x3Polynomial (A, B, C, nReal, roots (1:3,1:2), sim_printInfo, sim_printUnit)
!
!
!   ...Determine all (relative) accuracies. In case an analytic root is equal to zero
!      use the absolute accuracy, i.e. the absolute value of the ontained computational
!      root.
!
!  
     do j = Re,Im
        do i = 1,3

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
     exact (1:3,1:2) = sim_rootsAnalytical       (1:3,1:2,n)
     diffs (1:3,1:2) = sim_rootsRelativeAccuracy (1:3,1:2,n)

     write (prU,prF) ' '
     write (prU,prF) ' Exact  (analytical)  Root 1 (x + iy) = ',exact (1,Re),' + ',exact (1,Im),' i'
     write (prU,prF) ' Exact  (analytical)  Root 2 (x + iy) = ',exact (2,Re),' + ',exact (2,Im),' i'
     write (prU,prF) ' Exact  (analytical)  Root 3 (x + iy) = ',exact (3,Re),' + ',exact (3,Im),' i'
     write (prU,prF) ' '
     write (prU,prF) ' x3 Polynomial Solver Root 1 (x + iy) = ',roots (1,Re),' + ',roots (1,Im),' i'
     write (prU,prF) ' x3 Polynomial Solver Root 2 (x + iy) = ',roots (2,Re),' + ',roots (2,Im),' i'
     write (prU,prF) ' x3 Polynomial Solver Root 3 (x + iy) = ',roots (3,Re),' + ',roots (3,Im),' i'
     write (prU,prF) ' '
     write (prU,prF) ' Relative Accuracy Root 1 (Re and Im) = ',diffs (1,Re),'   ',diffs (1,Im)
     write (prU,prF) ' Relative Accuracy Root 2 (Re and Im) = ',diffs (2,Re),'   ',diffs (2,Im)
     write (prU,prF) ' Relative Accuracy Root 3 (Re and Im) = ',diffs (3,Re),'   ',diffs (3,Im)
     write (prU,prF) ' '
     write (prU,prF) ' -----------------------------------------------------------------------'

  end do
!
!
!   ...Determine overall relative accuracy (i.e. worst root) for each polynomial.
!
!  
  sim_x3RelativeAccuracy (:) = maxval (maxval (sim_rootsRelativeAccuracy, dim = 1), dim = 1)

  do n = 1, sim_numberOfx3Polynomials
     write (prU,prF) ' Relative Accuracy of x3 Polynomial = ', sim_x3RelativeAccuracy (n)
  end do
!
!
!   ...Compare the x3 polynomial relative accuracies obtained with the corresponding
!      upper bounds and decide, if the test passed.
!
!  
  perfect = all (sim_x3RelativeAccuracy (:) <= sim_maxRelativeAccuracy (:))
!
!
!   ...Ready!
!
!
  return
end subroutine sim_x3PolynomialsRootsTest
