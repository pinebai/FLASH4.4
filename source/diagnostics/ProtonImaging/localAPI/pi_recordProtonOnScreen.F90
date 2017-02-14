!!****if* source/diagnostics/ProtonImaging/localAPI/pi_recordProtonOnScreen
!!
!! NAME
!!
!!  pi_recordProtonOnScreen
!!
!! SYNOPSIS
!!
!!  call pi_recordProtonOnScreen (real    (in)            :: px,
!!                                real    (in)            :: py,
!!                                real    (in)            :: pz,
!!                                real    (in)            :: vx,
!!                                real    (in)            :: vy,
!!                                real    (in)            :: vz,
!!                                real    (in)            :: Jv,
!!                                real    (in)            :: Kx,
!!                                real    (in)            :: Ky,
!!                                real    (in)            :: Kz,
!!                                integer (in)            :: detector,
!!                                logical (out), optional :: onScreen
!!                                real    (out), optional :: sx,
!!                                real    (out), optional :: sy,
!!                                real    (out), optional :: sz)
!!
!! DESCRIPTION
!!
!!  Records the proton on the specified detector screen. The proton has an initial
!!  position and direction (velocity). The routine records the (x,y) location on
!!  the screen in terms of local screen coordinates as well as certain diagnostic
!!  variables (currently Jv,Kx,Ky,Kz). The (x,y) screen coordinates are rescaled to
!!  the square screen side length, such that all protons within the screen have screen
!!  coordinate paris (x,y) within the range [0,1]. Protons falling outside the screen
!!  are also saved for eventual further diagnostic studies. Protons which miss the
!!  screen plane entirely are not recorded.
!!
!!  Procedure:
!!
!!                                    |\
!!                                    | \
!!                                    |  \                   C = detector center
!!                                    |   \
!!                                    | uY|\                 n = normal unit vector
!!                                    |   | \
!!                                    |   |  \         uX,uY,n = form a right hand rule
!!                      detector      |   |   |                  orthogonal vector system:
!!                                    |   C---|---> n                  n = uY x uX
!!                       screen       \  / \  |
!!                                     \/   \ |              P = position of proton
!!                                     /\  uX\|
!!                                    /  \    |              v = velocity vector of proton
!!                                   /    \   |
!!                                  /      \  |          (x,y) = local screen coordinates
!!                                 /        \ |                  (basis = uX,uY) where proton
!!                                /          \|                   will hit the screen plane
!!                               /
!!                              /                   (sx,sy,sz) = global screen coordinates,
!!                             /                                 calculated only, if proton hits
!!                            /                                  screen plane. (0,0,0) otherwise.
!!                           /
!!                          P---> v ----------- (x,y)
!!                                         (sx,sy,sz)
!!
!!  From geometry:
!!
!!                    i) [(C - P) dot n] / [v dot n] > 0  ->  proton hits screen plane
!!
!!                   ii) x = uY dot [v cross (C - P)] / [v dot n]
!!
!!                  iii) y = uX dot [(C - P) cross v] / [v dot n]
!!
!!                   iv) if [v dot n] = 0  ->  parallel to screen plane, no hit
!!
!!
!! ARGUMENTS
!!
!!  px       : position global x-coordinate of the proton
!!  py       : position global y-coordinate of the proton
!!  pz       : position global z-coordinate of the proton
!!  vx       : velocity x-coordinate component of the proton
!!  vy       : velocity y-coordinate component of the proton
!!  vz       : velocity z-coordinate component of the proton
!!  Jv       : diagnostic variable
!!  Kx       : diagnostic variable
!!  Ky       : diagnostic variable
!!  Kz       : diagnostic variable
!!  detector : detector screen number on which to record the proton
!!  onScreen : if true, the proton hits the screen
!!  sx       : screen global x-coordinate of the proton (if screen plane is hit, 0 otherwise)
!!  sy       : screen global y-coordinate of the proton (if screen plane is hit, 0 otherwise)
!!  sz       : screen global z-coordinate of the proton (if screen plane is hit, 0 otherwise)
!!
!! NOTES
!!
!!***

subroutine pi_recordProtonOnScreen (px, py, pz,                &
                                    vx, vy, vz,                &
                                    Jv, Kx, Ky, Kz,            &
                                    detector,                  &
                                                    onScreen,  &
                                                    sx, sy, sz )

  implicit none

  real,    intent (in)            :: px, py, pz
  real,    intent (in)            :: vx, vy, vz
  real,    intent (in)            :: Jv, Kx, Ky, Kz
  integer, intent (in)            :: detector
  logical, intent (out), optional :: onScreen
  real,    intent (out), optional :: sx, sy, sz

  if (present (sx)) sx = 0.0
  if (present (sy)) sy = 0.0
  if (present (sz)) sz = 0.0
  if (present (onScreen)) onScreen = .false.

  return
end subroutine pi_recordProtonOnScreen
