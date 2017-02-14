!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_traceBlockRays1DRec
!!
!! NAME
!!
!!  ed_traceBlockRays1DRec
!!
!! SYNOPSIS
!!
!!  call ed_traceBlockRays1DRec (real    (in)    :: timeStep,
!!                               integer (in)    :: rayFirst
!!                               integer (in)    :: rayLast,
!!                               integer (in)    :: iminBlock,
!!                               integer (in)    :: imaxBlock,
!!                               real    (in)    :: xminBlock,
!!                               real    (in)    :: xmaxBlock,
!!                               real    (in)    :: deltaX,
!!                               real    (in)    :: deltaInvX,
!!                               logical (in)    :: blockReflectMinX,
!!                               logical (in)    :: blockReflectMaxX,
!!                               real    (inout) :: cellEnergyDepot (:))
!!
!! DESCRIPTION
!!
!!  Traces the movement of the current collection of active rays through one block for
!!  those geometries consisting formally of 1D rectangular grids (cartesian + spherical).
!!  On exit, each ray has either:
!!
!!            i)  reached a different (yet unknown) block
!!           ii)  has been absorbed -> is nonexistent
!!          iii)  has reached the domain boundary and exited -> is nonexistent.
!!
!! ARGUMENTS
!!
!!  timeStep         : current timestep value
!!  rayFirst         : first ray index to be considered
!!  rayLast          : last ray index to be considered
!!  iminBlock        : minimum cell i-index limit defining the interior block
!!  imaxBlock        : maximum cell i-index limit defining the interior block
!!  xminBlock        : minimum x-coordinate limit of the block
!!  xmaxBlock        : maximum x-coordinate limit of the block
!!  deltaX           : the cell's x-dimension
!!  deltaInvX        : inverse of the cell's x-dimension
!!  blockReflectMinX : is the block boundary on the minimum x-side reflective ?
!!  blockReflectMaxX : is the block boundary on the maximum x-side reflective ?
!!  cellEnergyDepot  : array collecting the ray energy deposition for each cell
!!
!! NOTES
!!
!!  The code allows for threading to be used on the outer ray trace loop.
!!  The paths of the rays are computed using the geometric optics approximation. 
!!
!!***

subroutine ed_traceBlockRays1DRec (timeStep,                          &
                                   rayFirst,  rayLast,                &
                                   iminBlock, imaxBlock,              &
                                   xminBlock, xmaxBlock,              &
                                   deltaX,                            &
                                   deltaInvX,                         &
                                   blockReflectMinX,                  &
                                   blockReflectMaxX,                  &
                                                      cellEnergyDepot ) 

  implicit none

  real,    intent (in)    :: timeStep
  integer, intent (in)    :: rayFirst,  rayLast   
  integer, intent (in)    :: iminBlock, imaxBlock
  real,    intent (in)    :: xminBlock, xmaxBlock
  real,    intent (in)    :: deltaX
  real,    intent (in)    :: deltaInvX
  logical, intent (in)    :: blockReflectMinX
  logical, intent (in)    :: blockReflectMaxX
  real,    intent (inout) :: cellEnergyDepot (iminBlock:imaxBlock)

  cellEnergyDepot (iminBlock:imaxBlock) = 0.0

  return
end subroutine ed_traceBlockRays1DRec
