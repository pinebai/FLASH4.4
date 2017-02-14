!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_traceRays1DRec
!!
!! NAME
!!
!!  ed_traceRays1DRec
!!
!! SYNOPSIS
!!
!!  call ed_traceRays1DRec (real, intent (in) :: timeStep)
!!
!! DESCRIPTION
!!
!!  Processes all rays on the collection of blocks on the current processor for those
!!  geometries consisting formally of 1D rectangular grids (cartesian + spherical).
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

subroutine ed_traceRays1DRec (timeStep)

  implicit none

  real, intent (in) :: timeStep

  return
end subroutine ed_traceRays1DRec
