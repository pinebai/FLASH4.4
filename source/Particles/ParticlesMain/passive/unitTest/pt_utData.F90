!!****if* source/Particles/ParticlesMain/passive/unitTest/pt_utData
!!
!! NAME
!!
!!  pt_utData
!!
!! SYNOPSIS
!!  use pt_utData
!!
!! DESCRIPTION
!!
!!  Store the simulation data for unitTesting of Particles
!!   
!! ARGUMENTS
!!
!! PARAMETERS   
!!      Described in the Config file
!!
!!
!!***

module pt_utData

  implicit none
  
  real, save    :: pt_utA0, pt_utA1, pt_utInitialSimTime

  logical,save :: pt_utAnalyticParticlePositions, pt_utFakeMapMeshToParticles

end module pt_utData
