!!****if* source/Grid/localAPI/gr_mpoleSetRadialBinData
!!
!! NAME
!!
!!  gr_mpoleSetRadialBinData
!!
!! SYNOPSIS
!!
!!  gr_mpoleSetRadialBinData ()
!!
!! DESCRIPTION
!!
!!  Sets all radial bin data that is necessary for moment and potential evaluation.
!!  It currently computes the maximum radial value for each bin and the bin
!!  damping factors. The damping factors are used to avoid possible overflow
!!  or underflow when evaluating the regular and irregular solid harmonic
!!  functions via recursion. When using these factors the recursion becomes
!!  'normalized' in the sense that the highest L-valued solid harmonic functions
!!  are close to unity.
!!
!!  For a given radial bin index Q outside of the inner zone, the maximum and
!!  minimum radius associated with that bin is calculated using the radial zone
!!  data provided. From these radii the damping factors are evaluated as follows,
!!  using the maximum angular momentum value L of the problem:
!!
!!                  regular solid harmonic = (  L-th root of L!) / rmax
!!                irregular solid harmonic = (L+1-th root of L!) / rmin
!!
!!  The x,y,z coordinates will be multiplied by the corresponding damping
!!  factors during application of the recursions. The approximations used
!!  for the roots of the factorials follow from Stirlings factorial approximation
!!  and are:
!!
!!      (  L-th root of L!) ~ (L/e) * (  2L-th root of 2*pi*L)
!!      (L+1-th root of L!) ~ (L/e) * (2L+2-th root of 2*pi*e^2/L)
!!
!!  Here e = 2.7182818...and pi = 3.14159....
!!
!!  For a given radial bin index Q inside the inner zone, there is only one
!!  radius (call it rmax) associated with it, which is sitting in the array
!!  'gr_mpoleInnerZoneDrRadii' in units of the inner zone atomic distance. The
!!  damping factors are evaluated as follows, using the maximum angular momentum
!!  value L of the problem:
!!
!!                  regular solid harmonic = (  L-th root of L!) / rmax
!!                irregular solid harmonic = (L+1-th root of L!) / rmax
!!
!!***

subroutine gr_mpoleSetRadialBinData ()

  implicit none
  
  return
end subroutine gr_mpoleSetRadialBinData
