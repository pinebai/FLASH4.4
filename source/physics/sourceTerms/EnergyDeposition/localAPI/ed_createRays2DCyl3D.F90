!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_createRays2DCyl3D
!!
!! NAME
!!
!!  ed_createRays2DCyl3D
!!
!! SYNOPSIS
!!
!!  call ed_createRays2DCyl3D (integer, intent (in) :: blockCount,
!!                             integer, intent (in) :: blockList (:),
!!                             real,    intent (in) :: timeStep,
!!                             real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates 3D cartesian rays and places them in their initial 3D wedge blocks for 2D cylindrical
!!  geometries. The 3D wedge blocks are considered an approximation to the true 3-dimensional
!!  cylindrical geometry.
!!
!!
!!     y axis
!!       |                                                     *  <--- upper wedge boundary
!!       |                                            *        |
!!       |                                   *        |        |
!!       |                          *        |        |        |
!!       |                 *        |        |        |        |
!!       |        *        |        |        |        |         
!!       O--------|--------|--------|--------|--------|------- X ---------------------> R or x
!!       |        *        |        |        |        |         
!!       |                 *        |        |        |        |  X = 2D cylindrical
!!       |                          *        |        |        |      domain boundary
!!       |                                   *        |        |
!!       |                                            *        |
!!       |                                                     *  <--- lower wedge boundary
!!
!!
!!
!!  The z-direction is the same for either the 2D cylindrical or the 3D cartesian picture.
!!  The ray line equations are 3D cartesian equations and we need to know at what angle they
!!  hit R-domain boundary in the (x,y) plane. All rays will be placed initially exactly
!!  along the wedge middle line (the R or x axis). The different incident angles of the rays
!!  on the 3-dimensional cylinder with radius R will be reflected in their different velocity
!!  components.
!!
!!  On exit, all rays hitting the domain boundary have been generated for the current processor.
!!  Their 2D cylindrical block ID's are not ordered as the outer loop is over all beams.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeStep       : current timestep value
!!  timeSimulation : current simulation time
!!
!! NOTES
!!
!!***

subroutine ed_createRays2DCyl3D (blockCount, blockList, timeStep, timeSimulation)

  implicit none

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation

  return
end subroutine ed_createRays2DCyl3D
