!!****if* source/Grid/localAPI/gr_mpoleRadialSampling
!!
!! NAME
!!
!!  gr_mpoleRadialSampling
!!
!! 
!! SYNOPSIS
!!
!!  gr_mpoleRadialSampling()
!!
!!
!! DESCRIPTION
!!
!!  Determines the radial sampling for accumulating the Moments.
!!  It defines the radial bins into which each Moment falls. If the
!!  radial sampling is chosen unwisely, lots of memory for holding
!!  empty Moment bins might be wasted. This routine calls the appropriate
!!  subroutines according to the geometry specified.
!!
!!
!!***

subroutine gr_mpoleRadialSampling()

  implicit none

  return
end subroutine gr_mpoleRadialSampling
