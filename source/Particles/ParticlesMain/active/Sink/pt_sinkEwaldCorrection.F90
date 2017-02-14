!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkEwaldCorrection
!!
!! NAME
!!
!!  pt_sinkEwaldCorrection
!!
!! SYNOPSIS
!!
!!  call pt_sinkEwaldCorrection(real(in)  :: x,
!!                              real(in)  :: y,
!!                              real(in)  :: z,
!!                              real(out) :: xcorr,
!!                              real(out) :: ycorr,
!!                              real(out) :: zcorr)
!!
!! DESCRIPTION
!!
!!  For a given distance (x,y,z), this function interpolates the Ewald correction
!!  table generated in pt_sinkGenerateEwald() and returns the force correction due
!!  to fully periodic boundary conditions.
!!
!! ARGUMENTS
!!
!!   x - x distance between point of force correction and acceleration-generating mass
!!
!!   y - y distance between point of force correction and acceleration-generating mass
!!
!!   z - z distance between point of force correction and acceleration-generating mass
!!
!!   xcorr - x component of the Ewald force correction
!!
!!   ycorr - y component of the Ewald force correction
!!
!!   zcorr - z component of the Ewald force correction
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2014
!!
!!***

subroutine pt_sinkEwaldCorrection(x, y, z, xcorr, ycorr, zcorr)

  use Particles_sinkData
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  real, intent(IN)  :: x, y, z
  real, intent(OUT) :: xcorr, ycorr, zcorr

  integer :: i, j, k, ip1, jp1, kp1
  real :: ireal, jreal, kreal, p, q, r, omp, omq, omr
  real :: c1, c2, c3, c4, c5, c6, c7, c8

  ! find the nearest smaller index in the Ewald field for the given coordinate (x,y,z)

  ireal = x * sink_EwaldDxI
  jreal = y * sink_EwaldDyI
  kreal = z * sink_EwaldDzI

  i = floor(ireal)
  j = floor(jreal)
  k = floor(kreal)

  if ((i .lt. 0) .or. (i .ge. sink_EwaldNx)) then
    call Driver_abortFlash("pt_sinkEwaldCorrection: Ewald field i out of limits")
  endif
  if ((j .lt. 0) .or. (j .ge. sink_EwaldNy)) then
    call Driver_abortFlash("pt_sinkEwaldCorrection: Ewald field j out of limits")
  endif
  if ((k .lt. 0) .or. (k .ge. sink_EwaldNz)) then
    call Driver_abortFlash("pt_sinkEwaldCorrection: Ewald field k out of limits")
  endif

  ! coefficients for linear interpolation

  p = ireal - i
  q = jreal - j
  r = kreal - k

  omp = 1.0 - p
  omq = 1.0 - q
  omr = 1.0 - r

  ip1 = min(i+1, sink_EwaldNx-1) ! min() to avoid out-of-index bounds if x = Lx/2;
  jp1 = min(j+1, sink_EwaldNy-1) ! interpolation terms vanish anyway in that case
  kp1 = min(k+1, sink_EwaldNz-1)

  c1 = omp * omq * omr
  c2 = omp * omq *   r
  c3 = omp *   q * omr
  c4 = omp *   q *   r
  c5 =   p * omq * omr
  c6 =   p * omq *   r
  c7 =   p *   q * omr
  c8 =   p *   q *   r

  xcorr = ( c1 * sink_EwaldFieldX(i  , j  , k  ) + &
            c2 * sink_EwaldFieldX(i  , j  , kp1) + &
            c3 * sink_EwaldFieldX(i  , jp1, k  ) + &
            c4 * sink_EwaldFieldX(i  , jp1, kp1) + &
            c5 * sink_EwaldFieldX(ip1, j  , k  ) + &
            c6 * sink_EwaldFieldX(ip1, j  , kp1) + &
            c7 * sink_EwaldFieldX(ip1, jp1, k  ) + &
            c8 * sink_EwaldFieldX(ip1, jp1, kp1)   )

  ycorr = ( c1 * sink_EwaldFieldY(i  , j  , k  ) + &
            c2 * sink_EwaldFieldY(i  , j  , kp1) + &
            c3 * sink_EwaldFieldY(i  , jp1, k  ) + &
            c4 * sink_EwaldFieldY(i  , jp1, kp1) + &
            c5 * sink_EwaldFieldY(ip1, j  , k  ) + &
            c6 * sink_EwaldFieldY(ip1, j  , kp1) + &
            c7 * sink_EwaldFieldY(ip1, jp1, k  ) + &
            c8 * sink_EwaldFieldY(ip1, jp1, kp1)   )

  zcorr = ( c1 * sink_EwaldFieldZ(i  , j  , k  ) + &
            c2 * sink_EwaldFieldZ(i  , j  , kp1) + &
            c3 * sink_EwaldFieldZ(i  , jp1, k  ) + &
            c4 * sink_EwaldFieldZ(i  , jp1, kp1) + &
            c5 * sink_EwaldFieldZ(ip1, j  , k  ) + &
            c6 * sink_EwaldFieldZ(ip1, j  , kp1) + &
            c7 * sink_EwaldFieldZ(ip1, jp1, k  ) + &
            c8 * sink_EwaldFieldZ(ip1, jp1, kp1)   )

  return

end subroutine pt_sinkEwaldCorrection
