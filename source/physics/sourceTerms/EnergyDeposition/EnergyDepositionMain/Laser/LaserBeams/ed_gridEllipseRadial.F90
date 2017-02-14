!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_gridEllipseRadial
!!
!! NAME
!!
!!  ed_gridEllipseRadial
!!
!! SYNOPSIS
!!
!!  call ed_gridEllipseRadial (integer (inout) :: nGridPoints,
!!                             integer (out)   :: nTicsRadial,
!!                             integer (out)   :: nTicsAngular,
!!                             real    (out)   :: deltaRadial,
!!                             real    (out)   :: deltaAngular)
!!
!! DESCRIPTION
!!
!!  Given an elliptical boundary shape, the routine places a radial grid onto it, such that the
!!  number of grid points inside the ellipse is close to the number of such grid points wanted.
!!  The radial grid is formed by placing equally spaced concentrical inner ellipses (radial tics)
!!  inside the boundary ellipse and subdividing evenly the angular part by the same angle (angular
!!  tics).
!!
!!  Features:
!!
!!   1) The routine only defines a grid in case the number of requested grid points is > 1.
!!      Nothing is done, if this number is =< 1.
!!
!!   2) The routine tries to balance the number of radial and angular tics in such a way that
!!      space inside the elliptical boundary is approximately evenly covered with grid points.
!!      This is done by applying the procedure explained below.
!!
!!  Procedure (for feature 2):
!!
!!      This procedure tries to arrange a collection of N total grid points on a radial grid
!!      such that the grid points are approximately evenly distributed in space. Consider a
!!      circular boundary shape and one angular section between consecutive radial spikes, where
!!      the grid points have been numbered:
!!
!!
!!                                          L4
!!                                   4--------------4
!!                              D ->  \     L3     /
!!                                     3----------3
!!                                D ->  \   L2   /
!!                                       2------2
!!                                  D ->  \ L1 /
!!                                         1--1
!!                                    D ->  \/
!!
!!      The distances between consecutive grid points on one spike are all equal to D and
!!      the distances between two adjacent grid points on the same ellipse are L1, L2, L3, etc.
!!      The Li (i=1,2,3,...) distances are all different, increasing in size with index i.
!!      An approximate even distribution of grid points will be achieved, if the average of the
!!      sum of all Li will be equal to D. Since each Li can be expressed in terms of D from
!!      the simple regular polygon formula:
!!
!!                                Li = 2 * i * D * sin (pi / Na)
!!
!!      where Na is the number of angular slices (the number of polygonal sides). Taking the
!!      sum over all i from 1 to Nr, where Nr is the number of radial tics, and using the fact
!!      that Nr * Na = N, we are lead to an equation in terms of Nr of the form:
!!
!!                                sin (pi * Nr / N) = 1 / (Nr + 1)
!!
!!      Rather than solving this equation exactly for each possible N, we will use the large
!!      N solution, which, when using the nearest integer value, proves to be useful for all N.
!!      Using the solution ansatz Nr = A * sqrt (N), with A to be determined, we have:
!!
!!                          sin (pi * A / sqrt (N)) = 1 / (A * sqrt (N) + 1)
!!
!!      which, for large N and using the approximation sin (x) = x for small x, gives:
!!
!!                             pi * A / sqrt (N) = 1 / (A * sqrt (N))
!!
!!      from which we obtain:
!!
!!                               A = sqrt (1 / pi) = 0.56418958....
!!
!!      and hence the large N solution is:
!!
!!                                  Nr = sqrt (N / pi)
!!
!!      of which we use the nearest integer value. Once Nr is found, we determine the number of
!!      angular slices as:
!!
!!                                  Na = nint (N / Nr)
!!
!!      and the final number of grid points will be:
!!
!!                                 N (final) = Nr * Na
!!
!!
!! ARGUMENTS
!!
!!  nGridPoints     : on input  -> the number of grid points wanted
!!                    on output -> the number of grid points that will result from the radial grid
!!  nTicsRadial     : the number of radial tics
!!  nTicsAngular    : the number of angular tics
!!  deltaRadial     : the radial tics spacing (as fraction of radius)
!!  deltaAngular    : the angular tics spacing (as fraction of 2pi)
!!
!! NOTES
!!
!!  Since the radial grid is defined in terms of concentric inner ellipses, there is no
!!  need to specify the boundary ellipses detailed structure, i.e. the length of both
!!  semiaxes. This contrasts the case of the square grid, which needs at least the ratio
!!  between both semiaxes in order to establish the square grid.
!!
!!***

subroutine ed_gridEllipseRadial (nGridPoints,               &
                                              nTicsRadial,  &
                                              nTicsAngular, &
                                              deltaRadial,  &
                                              deltaAngular  )

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"

  integer, intent (inout) :: nGridPoints
  integer, intent (out)   :: nTicsRadial
  integer, intent (out)   :: nTicsAngular
  real,    intent (out)   :: deltaRadial
  real,    intent (out)   :: deltaAngular
!
!
!     ...Set up the radial grid. 
!
!
  if (nGridPoints > 1) then

      nTicsRadial  = nint (sqrt (real (nGridPoints) / PI ))
      nTicsAngular = nint (real (nGridPoints) / real (nTicsRadial))
      nGridPoints  = nTicsRadial * nTicsAngular
      deltaRadial  = 1.0 / real (nTicsRadial)
      deltaAngular = (PI + PI) / real (nTicsAngular)

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_gridEllipseRadial
