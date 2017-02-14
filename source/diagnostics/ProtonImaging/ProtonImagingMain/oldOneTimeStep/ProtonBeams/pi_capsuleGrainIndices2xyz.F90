!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonBeams/pi_capsuleGrainIndices2xyz
!!
!! NAME
!!
!!  pi_capsuleGrainIndices2xyz
!!
!! SYNOPSIS
!!
!!  call pi_capsuleGrainIndices2xyz (integer, intent (in)  :: i,
!!                                   integer, intent (in)  :: j,
!!                                   integer, intent (in)  :: k,
!!                                   real,    intent (in)  :: S,
!!                                   real,    intent (out) :: x,
!!                                   real,    intent (out) :: y,
!!                                   real,    intent (out) :: z)
!!
!! DESCRIPTION
!!
!!  This routine maps a capsule grain index triple (i,j,k) to a local coordinate triple
!!  (x,y,z), with local coordinate origin at the capsule center. The triple (i,j,k) denotes
!!  the labelling index of a cube inside the capsule sphere and the resulting coordinate
!!  triple (x,y,z) is the location of the cube center in local coordinates. Each index
!!  triple (i,j,k) stands for the outermost (farthest from the capsule center) coordinate
!!  location (i*S, j*S, k*S) of a cube vertex, where S is the cube side length. The center
!!  location of the cube would then be given as:
!!
!!                                    (i +/- 1/2) * S
!!                                    (j +/- 1/2) * S
!!                                    (k +/- 1/2) * S
!!
!!  where the +/- sign applies, if the signs of the i,j,k are -/+. A zero index in i,j,k
!!  leads to a zero local coordinate.
!!
!! ARGUMENTS
!!
!!  i : i-index of capsule grain index triple (i,j,k)
!!  j : j-index of capsule grain index triple (i,j,k)
!!  k : k-index of capsule grain index triple (i,j,k)
!!  S : the cube side length
!!  x : the capsule grain x-coordinate
!!  y : the capsule grain y-coordinate
!!  z : the capsule grain z-coordinate
!!
!!***

subroutine pi_capsuleGrainIndices2xyz (i,j,k,S,x,y,z)

  implicit none

  integer, intent (in)  :: i,j,k
  real,    intent (in)  :: S
  real,    intent (out) :: x,y,z
!
!
!     ...Calculate the coordinate triple (x,y,z).
!
!
  if (i > 0) then
      x = (real (i) - 0.5) * S
  else if (i < 0) then
      x = (real (i) + 0.5) * S
  else
      x = 0.0
  end if

  if (j > 0) then
      y = (real (j) - 0.5) * S
  else if (j < 0) then
      y = (real (j) + 0.5) * S
  else
      y = 0.0
  end if

  if (k > 0) then
      z = (real (k) - 0.5) * S
  else if (k < 0) then
      z = (real (k) + 0.5) * S
  else
      z = 0.0
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine pi_capsuleGrainIndices2xyz
