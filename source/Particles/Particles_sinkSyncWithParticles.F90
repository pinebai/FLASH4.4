!!****f* source/Particles/Particles_sinkSyncWithParticles
!!
!! NAME
!!
!!  Particles_sinkSyncWithParticles
!!
!! SYNOPSIS
!!
!!  call Particles_sinkSyncWithParticles(logical(in) :: sink_to_part)
!!
!! DESCRIPTION
!!
!!  Synchronizes global particle array with sink particle array.
!!
!! ARGUMENTS
!!
!!   sink_to_part -  logical flag indicating whether to sync with global
!!                   particle array or not
!!
!! NOTES
!!
!!   written by John Bachan, 2012
!!   modified to include off-domain support, Christoph Federrath, 2015
!!
!!***

subroutine Particles_sinkSyncWithParticles(sink_to_part)
  implicit none
  logical, intent(in) :: sink_to_part
end subroutine Particles_sinkSyncWithParticles
