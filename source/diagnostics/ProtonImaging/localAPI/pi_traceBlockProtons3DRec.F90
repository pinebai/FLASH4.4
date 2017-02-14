!!****if* source/diagnostics/ProtonImaging/localAPI/pi_traceBlockProtons3DRec
!!
!! NAME
!!
!!  pi_traceBlockProtons3DRec
!!
!! SYNOPSIS
!!
!!  call pi_traceBlockProtons3DRec (integer (in)    :: protonFirst
!!                                  integer (in)    :: protonLast,
!!                                  integer (in)    :: iminBlock,
!!                                  integer (in)    :: imaxBlock,
!!                                  integer (in)    :: jminBlock,
!!                                  integer (in)    :: jmaxBlock,
!!                                  integer (in)    :: kminBlock,
!!                                  integer (in)    :: kmaxBlock,
!!                                  real    (in)    :: xminBlock,
!!                                  real    (in)    :: xmaxBlock,
!!                                  real    (in)    :: yminBlock,
!!                                  real    (in)    :: ymaxBlock,
!!                                  real    (in)    :: zminBlock,
!!                                  real    (in)    :: zmaxBlock,
!!                                  real    (in)    :: deltaX,
!!                                  real    (in)    :: deltaY,
!!                                  real    (in)    :: deltaZ,
!!                                  real    (in)    :: deltaInvX,
!!                                  real    (in)    :: deltaInvY,
!!                                  real    (in)    :: deltaInvZ,
!!                                  logical (in)    :: blockReflectMinX,
!!                                  logical (in)    :: blockReflectMaxX,
!!                                  logical (in)    :: blockReflectMinY,
!!                                  logical (in)    :: blockReflectMaxY,
!!                                  logical (in)    :: blockReflectMinZ,
!!                                  logical (in)    :: blockReflectMaxZ)
!!
!! DESCRIPTION
!!
!!  Traces the movement of the current collection of active protons through one block
!!  for those geometries consisting formally of 3D rectangular grids (cartesian).
!!  On exit, each proton has either:
!!
!!            i)  reached a different (yet unknown) block
!!           ii)  has reached the domain boundary and exited.
!!
!! ARGUMENTS
!!
!!  protonFirst      : first proton index to be considered
!!  protonLast       : last proton index to be considered
!!  iminBlock        : minimum cell i-index limit defining the interior block
!!  imaxBlock        : maximum cell i-index limit defining the interior block
!!  jminBlock        : minimum cell j-index limit defining the interior block
!!  jmaxBlock        : maximum cell j-index limit defining the interior block
!!  kminBlock        : minimum cell k-index limit defining the interior block
!!  kmaxBlock        : maximum cell k-index limit defining the interior block
!!  xminBlock        : minimum x-coordinate limit of the block
!!  xmaxBlock        : maximum x-coordinate limit of the block
!!  yminBlock        : minimum y-coordinate limit of the block
!!  ymaxBlock        : maximum y-coordinate limit of the block
!!  zminBlock        : minimum z-coordinate limit of the block
!!  zmaxBlock        : maximum z-coordinate limit of the block
!!  deltaX           : the cell's x-dimension
!!  deltaY           : the cell's y-dimension
!!  deltaZ           : the cell's z-dimension
!!  deltaInvX        : inverse of the cell's x-dimension
!!  deltaInvY        : inverse of the cell's y-dimension
!!  deltaInvZ        : inverse of the cell's z-dimension
!!  blockReflectMinX : is the block boundary on the minimum x-side reflective ?
!!  blockReflectMaxX : is the block boundary on the maximum x-side reflective ?
!!  blockReflectMinY : is the block boundary on the minimum y-side reflective ?
!!  blockReflectMaxY : is the block boundary on the maximum y-side reflective ?
!!  blockReflectMinZ : is the block boundary on the minimum z-side reflective ?
!!  blockReflectMaxZ : is the block boundary on the maximum z-side reflective ?
!!
!! NOTES
!!        
!!  The code allows for threading to be used on the outer proton trace loop.
!!  The paths of the protons are computed using the average Lorentz force for each cell. 
!!
!!***

subroutine pi_traceBlockProtons3DRec (protonFirst,  protonLast,          &
                                      iminBlock, imaxBlock,              &
                                      jminBlock, jmaxBlock,              &
                                      kminBlock, kmaxBlock,              &
                                      xminBlock, xmaxBlock,              &
                                      yminBlock, ymaxBlock,              &
                                      zminBlock, zmaxBlock,              &
                                      deltaX, deltaY, deltaZ,            &
                                      deltaInvX, deltaInvY, deltaInvZ,   &
                                      blockReflectMinX,                  &
                                      blockReflectMaxX,                  &
                                      blockReflectMinY,                  &
                                      blockReflectMaxY,                  &
                                      blockReflectMinZ,                  &
                                      blockReflectMaxZ                   ) 

  implicit none

  integer, intent (in)    :: protonFirst, protonLast   
  integer, intent (in)    :: iminBlock, imaxBlock
  integer, intent (in)    :: jminBlock, jmaxBlock
  integer, intent (in)    :: kminBlock, kmaxBlock
  real,    intent (in)    :: xminBlock, xmaxBlock
  real,    intent (in)    :: yminBlock, ymaxBlock
  real,    intent (in)    :: zminBlock, zmaxBlock
  real,    intent (in)    :: deltaX, deltaY, deltaZ
  real,    intent (in)    :: deltaInvX, deltaInvY, deltaInvZ
  logical, intent (in)    :: blockReflectMinX
  logical, intent (in)    :: blockReflectMaxX
  logical, intent (in)    :: blockReflectMinY
  logical, intent (in)    :: blockReflectMaxY
  logical, intent (in)    :: blockReflectMinZ
  logical, intent (in)    :: blockReflectMaxZ

  return
end subroutine pi_traceBlockProtons3DRec
