!!****f* source/Grid/Grid_mapParticlesToMesh
!!
!! NAME
!!  Grid_mapParticlesToMesh
!!
!! SYNOPSIS
!!
!!  call Grid_mapParticlesToMesh(real(INOUT) :: particles(part_props,maxParticlesPerProc),
!!                               integer(IN) :: part_props,
!!                               integer(IN) :: numParticles,
!!                               integer(IN) :: maxParticlesPerProc,
!!                               integer(IN) :: propPart,
!!                               integer(IN) :: varGrid, 
!!                      OPTIONAL,integer(IN) :: mode,
!!                      OPTIONAL,integer(IN) :: ptInfo)
!!
!! DESCRIPTION
!!
!! Routine to map a quantity defined for particles onto the  mesh.
!!
!! ARGUMENTS
!!  particles : List of particles. It is two dimensional real array, the first dimension
!!              represents each particle's properties, and second dimension is index to
!!              particles.
!!
!!  part_props : number of properties in the particles datastructure
!!
!!  numParticles : While coming in it contains the current number of particles mapped to
!!                      this processor. After all the data structure movement, the number
!!                      of local particles might change, and the new value is put back into it
!!  maxParticlesPerProc : The size of the buffer allocated for particles at compile time. This
!!                        is the largest number of particles that a simulation can have on one
!!                        processor at any time.
!!  varGrid:   Index of gridded variable to receive interpolated
!!                           quantity
!!  propPart:   Index of particle attribute to interpolate onto the mesh
!!  mode:       (Optional) If zero (default), zero varGrid first; if
!!                              nonzero, do not zero varGrid first
!!  ptInfo:     (Optional) additional info to pass to the Particles unit
!!                              in some calls, for debugging
!!
!! PARAMETERS
!! 
!!***

subroutine Grid_mapParticlesToMesh (particles,part_props,numParticles,&
     maxParticlesPerProc,propPart, varGrid, mode, ptInfo)
  
  implicit none
 
  integer,intent(IN) :: maxParticlesPerProc,numParticles, part_props
  real,dimension(part_props,numParticles),intent(INOUT) :: particles
  integer, INTENT(in) :: propPart, varGrid
  integer, INTENT(in), optional :: mode
  integer, INTENT(in), optional :: ptInfo

  return

end subroutine Grid_mapParticlesToMesh
