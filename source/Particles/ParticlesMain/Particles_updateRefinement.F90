!!****if* source/Particles/ParticlesMain/Particles_updateRefinement
!!
!! NAME
!!
!!  Particles_updateRefinement
!!
!! SYNOPSIS
!!
!!  Particles_updateRefinement(real(inout) :: oldLocalNumBlocks)
!!
!! DESCRIPTION
!!   This routine provides a hook into the particle data structure
!!   for the Grid. It is called during Grid_updateRefinement processing
!!   by the Grid. The routine passes the control right back to
!!   Grid, with Particles-specific data structures in the argument
!!   list, so that Grid can operate on them.
!!
!! ARGUMENTS
!!
!!    oldLocalNumBlocks :   number of blocks on a processor before 
!!                          refinement. 
!!
!! PARAMETERS
!!  
!!
!!***
#include "constants.h"
#include "Particles.h"
#include "Flash.h"

subroutine Particles_updateRefinement(oldLocalNumBlocks)

  use Particles_data, ONLY : particles, pt_numLocal, pt_maxPerProc,useParticles, &
       pt_posInitialized, pt_logLevel, pt_meshMe,&
       pt_indexList, pt_indexCount
  use Grid_interface,ONLY : Grid_moveParticles
  use Logfile_interface,ONLY : Logfile_stamp
  use Particles_interface, only: Particles_sinkMoveParticles
  
  implicit none 
  integer,intent(INOUT) :: oldLocalNumBlocks
  logical, parameter :: regrid = .true.

  if(.not.useParticles)return
  if(.not.pt_posInitialized) then
     if (pt_logLevel > PT_LOGLEVEL_WARN_USE) then
        if (pt_meshMe==MASTER_PE) then
           print*,'WARNING: Particles_updateRefinement was called while particles positions are not yet initialized!'
        end if
        call Logfile_stamp( &
             'WARNING: Called while particles positions are not yet initialized!','Particles_updateRefinement')
     end if
     return
  end if
  call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc, pt_numLocal, &
       pt_indexList, pt_indexCount,&
       regrid)
  call Particles_sinkMoveParticles(regrid)

  return
end subroutine Particles_updateRefinement
