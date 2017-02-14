!!****if* source/IO/IOParticles/IOParticles_data
!!
!! NAME
!!  IOParticles_data
!!
!! SYNOPSIS
!!
!!  IOParticles_data()
!!
!! DESCRIPTION 
!!  
!!  Holds all the IO data that is needed by the IOParticles unit  
!!
!!
!!***

module IOParticles_data

 
#include "constants.h"
#include "Flash.h"


  integer, save :: io_particleFileNumber
  integer, save :: io_particleFileIntervalStep
  real, save :: io_particleFileIntervalTime
  integer, save :: io_nextParticleFileStep
  real, save :: io_nextParticleFileTime
  
  real, save :: io_particleFileIntervalZ
  real, save :: io_nextParticleFileZ

  logical, save :: useParticles
  logical, save :: io_dumpParticleFileExist
  logical, save :: io_ptSupressSinglePrecision
  logical, save :: io_writeParticleSubset
  logical, save :: io_writeParticleAll

end module IOParticles_data
