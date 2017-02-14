!!****if* source/IO/IOParticles/IO_readParticles
!!
!! NAME
!!
!! IO_readParticles
!!
!!
!! SYNOPSIS
!!
!! IO_readParticles()
!!                  integer(in) :: numProcs)
!!
!!
!!
!! DESCRIPTION
!!
!!    This routine reads in the particle data from a checkpoint file.  
!!    This is a general routine
!!    which will then call io_readParticleData which is a routine specific
!!    to either HDF5 or parallel netcdf or a plain, nonparallel Fortran write.
!!
!!
!! ARGUMENTS
!!
!!      
!!
!!      numProcs :                total number of processors
!!
!!
!! NOTES
!!
!!***


subroutine IO_readParticles()

  use IOParticles_data, ONLY : useParticles
  use io_ptInterface, ONLY : io_ptReadParticleData
  use Particles_interface, only: Particles_sinkSyncWithParticles
  implicit none

  call io_ptReadParticleData()
  call Particles_sinkSyncWithParticles(sink_to_part=.false.)
end subroutine IO_readParticles
