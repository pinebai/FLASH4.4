!!****f* source/Particles/Particles_updateAttributes
!!
!! NAME
!!
!!  Particles_updateAttributes
!!
!! SYNOPSIS
!!
!!  Particles_updateAttributes()
!!
!!  Particle attribute advancement routines.  Use this routine to 
!!  update user-define particle attributes beyond the usual
!!  particle attributes of position, velocity, block, tag, and mass.
!!  It is usually called from  IO_output. It makes sure that guardcells
!!  are current and, if necessary, that Eos has been applied to them
!!  before particles attributes are calculated from the grid.
!!   
!!  The attributes can be specified at runtime 
!!  using runtime parameters particle_attribute_1, particle_attribute_2
!!  etc.
!!
!! ARGUMENTS
!!  
!! NOTES
!!
!! The map between particle property and a mesh variable on which the
!! property is dependent has to be specified in the Config file of 
!! the Simulation directory for this routine to work right. Please
!! see the Config file of IsentropicVortex setup for an example, and
!! also see the Setup chapter of the User's Guide.
!! 
!!
!!  
!!
!!
!!***


subroutine Particles_updateAttributes()


  implicit none


end subroutine Particles_updateAttributes

