!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_traceRays2DRec
!!
!! NAME
!!
!!  ed_traceRays2DRec
!!
!! SYNOPSIS
!!
!!  call ed_traceRays2DRec (real, intent (in) :: timeStep)
!!
!! DESCRIPTION
!!
!!  Processes all rays on the collection of blocks on the current processor for those
!!  geometries consisting formally of 2D rectangular grids (cartesian + cylindrical).
!!  On exit, each ray in each block has either:
!!
!!            i)  reached a different (yet unknown) block
!!           ii)  has been absorbed -> is nonexistent
!!          iii)  has reached the domain boundary and exited -> is nonexistent.
!!
!! ARGUMENTS
!!
!!  timeStep : Current timestep value
!!
!! NOTES
!!
!!  The use of threading is possible for tracing the rays through each block. 
!!
!!***

subroutine ed_traceRays2DRec (timeStep)

  implicit none

  real, intent (in) :: timeStep

  return
end subroutine ed_traceRays2DRec
