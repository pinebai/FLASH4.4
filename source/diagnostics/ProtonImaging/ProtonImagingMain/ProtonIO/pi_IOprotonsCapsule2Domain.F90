!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonIO/pi_IOprotonsCapsule2Domain
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
  
  use Driver_interface,    ONLY : Driver_abortFlash

  use ProtonImaging_data,  ONLY : pi_beams,                    &
                                  pi_IOmaxPointsPerBlock,      &
                                  pi_IOmaxProtonCount,         &
                                  pi_IOprotonCount,            &
                                  pi_IOprotonPointCount,       &
                                  pi_IOprotonPoints,           &
                                  pi_IOprotonTags,             &
                                  pi_IOprotonWriteModulo,      &
                                  pi_protonCount,              &
                                  pi_protons

  implicit none

#include "constants.h"
#include "Flash.h"
#include "ProtonImaging.h"

  logical :: IOproton

  integer :: beam
  integer :: proton
  integer :: protonTag

  real    :: Cx, Cy, Cz
  real    :: discriminant
  real    :: invMagV
  real    :: Ix, Iy, Iz
  real    :: Px, Py, Pz
  real    :: r
  real    :: Rx, Ry, Rz
  real    :: Rsqr
  real    :: UdotR
  real    :: Ux, Uy, Uz
  real    :: Vx, Vy, Vz
  real    :: w
!
!
!     ...Return immediately, if the 2 IO proton points cannot be stored.
!
!
  if (pi_IOmaxPointsPerBlock < 2) then
      print *, '[pi_IOprotonsCapsule2Domain] The capsule -> domain IO protons need 2 points, &
                 but have only ', pi_IOmaxPointsPerBlock, ' in buffer space'
      return
  end if
!
!
!     ...Loop over all protons on current processor.
!
!
  pi_IOprotonCount = 0

  do proton = 1, pi_protonCount

     protonTag  = int (pi_protons (PROTON_TAGS,proton))

     IOproton = mod (protonTag, pi_IOprotonWriteModulo) == 0
     if (IOproton) then
         IOproton = pi_IOprotonCount < pi_IOmaxProtonCount
     end if

     if (IOproton) then

         Px = pi_protons (PROTON_POSX,proton)
         Py = pi_protons (PROTON_POSY,proton)
         Pz = pi_protons (PROTON_POSZ,proton)
         Vx = pi_protons (PROTON_VELX,proton)
         Vy = pi_protons (PROTON_VELY,proton)
         Vz = pi_protons (PROTON_VELZ,proton)

         beam = int (pi_protons (PROTON_BEAM,proton))

         Cx = pi_beams (beam) % capsuleX
         Cy = pi_beams (beam) % capsuleY
         Cz = pi_beams (beam) % capsuleZ
         r  = pi_beams (beam) % capsuleRadius

         Rx = Px - Cx
         Ry = Py - Cy
         Rz = Pz - Cz

         invMagV = 1.0 / sqrt (Vx * Vx + Vy * Vy + Vz * Vz)

         Ux = Vx * invMagV
         Uy = Vy * invMagV
         Uz = Vz * invMagV

         Rsqr  = Rx * Rx + Ry * Ry + Rz * Rz
         UdotR = Ux * Rx + Uy * Ry + Uz * Rz

         discriminant = max (0.0, UdotR * UdotR - Rsqr + r * r)
         w = - UdotR + sqrt (discriminant)

         Ix = Px + w * Ux
         Iy = Py + w * Uy
         Iz = Pz + w * Uz

         pi_IOprotonCount = pi_IOprotonCount + 1

         pi_IOprotonTags       (   pi_IOprotonCount       ) = protonTag
         pi_IOprotonPointCount (   pi_IOprotonCount       ) = 2
         pi_IOprotonPoints     (1, pi_IOprotonCount, IAXIS) = Ix
         pi_IOprotonPoints     (1, pi_IOprotonCount, JAXIS) = Iy
         pi_IOprotonPoints     (1, pi_IOprotonCount, KAXIS) = Iz
         pi_IOprotonPoints     (2, pi_IOprotonCount, IAXIS) = Px
         pi_IOprotonPoints     (2, pi_IOprotonCount, JAXIS) = Py
         pi_IOprotonPoints     (2, pi_IOprotonCount, KAXIS) = Pz

     end if

  end do
!
!
!    ...Ready!
!
!
  return
end subroutine pi_IOprotonsCapsule2Domain
