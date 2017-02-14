!!****if* source/Particles/ParticlesMain/active/DPD/pt_dpdFinalize
!!
!! NAME
!!    pt_dpdFinalize
!!
!! SYNOPSIS
!!
!!    pt_dpdFinalize()
!!
!! DESCRIPTION
!!    Local finalizeialization for particle-in-cell implementation -stub
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!***

subroutine pt_dpdFinalize()

  use pt_dpdData, ONLY : particles2

  implicit none
  
  deallocate(particles2)

End subroutine pt_dpdFinalize


