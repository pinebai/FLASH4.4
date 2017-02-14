!!****f* source/IO/IO_writeParticles
!!
!! NAME
!!
!! IO_writeParticles
!!
!!
!! SYNOPSIS
!!
!! IO_writeParticles(logical(in) :: particlesToCheckpoint)
!! 
!!
!!
!!
!! DESCRIPTION
!!
!!     This routine writes out the particle data.  This is a general routine
!!     which will then call io_writeParticleData which is a routine specific
!!     to either hdf5 or parallel netcdf, or a plain old fortran write.
!!
!!     Particle data is written in two places.  First, particle data must be
!!     included in the checkpoint (restart) files in order to capture the 
!!     entire state of the simulation and to be able to restart.  Often, however,
!!     particle data needs to be written very frequently and is written to its 
!!     own particle file without any other mesh data. It is possible to 
!!     separate the particles from the 
!!     checkpoint files because particles are not associated with the mesh data.
!!
!!     Particles are written in double precision to a checkpoint file and
!!     single precision in particle plotfiles  
!!
!!     The functionality of writing the particle data to a checkpoint file or
!!     a particle plotfile is the same so to eliminate code duplication we 
!!     have added an argument to the IO_writeParticles interface to indicate
!!     if we are writing particles to a checkpoint file or to a particle file.
!!     (see below)
!!
!! ARGUMENTS
!!
!!
!!      particlesToCheckpoint -  logical value - .true. if particles are
!!                               written to a checkpoint. .false. if particles
!!                               are written to a particle plot file
!!    
!!
!! NOTES
!!   To control particle output there are a few different parameters
!!   Check the flash online documentation or users guide for a more detailed
!!   explanation   
!!
!!   particleFileNumber - integer.  current particle file number
!! 
!!   particleFileIntervalTime - Amount of simulation time in seconds
!!                              between particle file dumps
!!
!!   particleFileIntervalStep - Number of steps between particle file dumps
!!
!!   
!!
!!***


subroutine IO_writeParticles( particlesToCheckpoint)

  implicit none

  logical, intent(in) :: particlesToCheckpoint

  return

end subroutine IO_writeParticles
