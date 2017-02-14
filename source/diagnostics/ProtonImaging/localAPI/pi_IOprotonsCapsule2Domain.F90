!!****if* source/diagnostics/ProtonImaging/localAPI/pi_IOprotonsCapsule2Domain
!!
!! NAME
!!
!!  pi_IOprotonsCapsule2Domain
!!
!! SYNOPSIS
!!
!!  call pi_IOprotonsCapsule2Domain ()
!!
!! DESCRIPTION
!!
!!  This routine stores beam capsule -> domain IO protons into the corresponding
!!  IO proton arrays. Two x,y,z positions need to be recorded: 1) the position
!!  of the IO proton on the capsule spherical surface and 2) the entering domain
!!  position of the IO proton. Together with the positions, the tags are also
!!  recorded.
!!
!!  Procedure:
!!
!!  When calling this routine, all protons must have been created and globally
!!  taged. They all have an entering domain position (Px,Py,Pz), an initial velocity
!!  vector (Vx,Vy,Vz), a beam identification label and a unique global tag. They
!!  do NOT have their position of origin in the beam capsule. Hence we need to work
!!  backwards in order to locate the position of origin by using the velocity vector
!!  and the beam capsule location and radius.
!!
!!                               |
!!              _                |  domain          P -> position vector of proton
!!             / \               |                  V -> velocity vector of proton
!!            /   I ************ P ---> V           C -> center vector of capsule
!!           /     \             |                  r -> radius of capsule
!!          |-r-C   |            |                  I -> intersection point on capsule
!!           \     /             |                  R -> vector (P - C)
!!            \   /                                 U -> unit velocity vector V/|V|
!!             \_/
!!
!!         beam capsule
!!
!!  After solving the two equations defining position I: 1) I = P + w * V and the
!!  sphere intersection equation 2) (Ix - Cx)^2 + (Iy - Cy)^2 + (Iz - Cz)^2 = r^2
!!  for w, we arrive at:
!!
!!                     I = P + [-U.R + sqrt ([U.R]^2 - R^2 + r^2)] U
!!
!!  If the discriminant under the square root is negative, this would mean a miss of
!!  the capsule sphere. In this case we assume the proton originated on the outermost
!!  tangential region of the capsule, which means that we set the square root equal
!!  to zero.
!!
!!***

subroutine pi_IOprotonsCapsule2Domain ()
  
  implicit none

  return
end subroutine pi_IOprotonsCapsule2Domain
