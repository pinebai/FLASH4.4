!!***f* source/Particles/Particles_mapToMeshOneBlk
!!
!! NAME
!!
!!  Particles_mapToMeshOneBlk
!!
!!
!! SYNOPSIS
!!
!!  call Particles_mapToMeshOneBlk(integer(IN),dimension(LOW:HIGH,MDIM)            :: blkLimitsGC,
!!                                 integer(IN),dimension(MDIM)                     :: guard,
!!                                 integer(IN)                                     :: blockID,
!!                                 real(IN)                               :: particles(NPART_PROPS,numParticles),
!!                                 integer(IN)                                     :: numParticles,
!!                                 integer(IN)                                     :: pt_attribute,
!!                                 real(INOUT),dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
!!                                                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
!!                                                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: buff,
!!                                 integer(IN),OPTIONAL                            :: particleOffset)
!! DESCRIPTION
!!
!! Routine which manages the mapping of all particles existing on a certain 
!! block to a temporary buffer (buff).  The particles in the particle buffer must be 
!! sorted prior to reaching this routine.  This is because all particles 
!! between the particle buffer limits of pStart and pEnd will be mapped to the 
!! temporary buffer.  The temporary buffer corresponds to a single block, and 
!! mapping may occur in the buffer's internal and guard cell grid points. 
!! The actual mapping takes place in the subroutine pt_mapOneParticle. 
!! 
!! ARGUMENTS
!!               blkLimitsGC:  Size of a block including guard cells.
!!               guard:  Number of guard cells.
!!               blockID:  ID of the block that will receive particle mapping.
!!               particles:  list of particles to map to this block
!!               numParticles:  length of second dim of particles argument
!!               pt_attribute:  Particle attribute that will be mapped to 
!!                              the mesh.
!!               buff:  Temporary buffer in which to accumulate particle mapping.
!!               particleOffset: offset of the dummy argument 'particles' within
!!                               the global 'particles' array in the Particles_data
!!                               module. Used only for debugging.
!!                               The caller should supply this ONLY if it really passes
!!                               an actual argument for the 'particles' dummy that is
!!                               a slice from the Particles_data 'particles' array!
!!
!! NOTES
!! 
!! In adaptive mesh simulations we may wish to identify those particles which 
!! smear across cells in blocks at different refinement levels. 
!! A motive is so that we can incorporate more elaborate particle mapping schemes
!! in future.  The routine which identifies these particles is named gr_ptParticleAtFcBdry.
!!
!!***

subroutine Particles_mapToMeshOneBlk(blkLimitsGC,guard,blockID,&
    particles,numParticles,pt_attribute,buff,particleOffset)

#include "constants.h"
#include "Flash.h"

 
  implicit none
  integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: blkLimitsGC
  integer, dimension(MDIM),intent(IN) :: guard
  integer, intent(IN) :: blockID, numParticles
  real, intent(in) :: particles(NPART_PROPS,numParticles)
  integer, intent(in) :: pt_attribute

  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                  INTENT(INOUT) :: buff
  integer, intent(in),OPTIONAL :: particleOffset

  return
end subroutine Particles_mapToMeshOneBlk

