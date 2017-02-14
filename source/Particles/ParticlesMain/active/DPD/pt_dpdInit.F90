!!****if* source/Particles/ParticlesMain/active/DPD/pt_dpdInit
!!
!! NAME
!!    pt_dpdInit
!!
!! SYNOPSIS
!!
!!    pt_dpdInit()
!!
!! DESCRIPTION
!!    Local initialization for particle-in-cell implementation -stub
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!***

#include "Flash.h"

subroutine pt_dpdInit()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use pt_dpdData, ONLY : pt_dpdUpdateCycle, particles2,pt_dpdLambda
  use Particles_data, ONLY : particles, pt_maxPerProc
  implicit none
  
  call RuntimeParameters_get('pt_dpdUpdateCycle', pt_dpdUpdateCycle)
  call RuntimeParameters_get('pt_dpdLambda', pt_dpdLambda)
  allocate(particles2(NPART_PROPS,pt_maxPerProc))
  particles2=particles
End subroutine pt_dpdInit


