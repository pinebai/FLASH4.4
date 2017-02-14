!!****f* source/IO/IO_readParticles
!!
!! NAME
!!
!! IO_readParticles
!!
!!
!! SYNOPSIS
!!
!! IO_readParticles(
!!                   integer(in) :: numProcs)
!!
!!
!!
!! DESCRIPTION
!!
!!    This routine reads in the particle data from a checkpoint file.  
!!      This is a general routine
!!      which will then call io_readParticleData which is a routine specific
!!      to either hdf5 or parallel netcdf or a plain, non parallel fortran write.
!!
!!
!! ARGUMENTS
!!
!!      
!!
!!      numProcs -                total number of processors
!!
!!
!! NOTES
!!
!!***


subroutine IO_readParticles()

  implicit none

  

  return

end subroutine IO_readParticles
