!!***f* source/Particles/ParticlesMapping/meshWeighting/Particles_mapToMeshOneBlk
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
!!                                 real(in)   ,dimension(NPART_PROPS,numParticles) :: particles,
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

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkBoundBox
  use Driver_interface, ONLY : Driver_abortFlash

!!$#ifdef FLASH_GRID_PARAMESH
!!$  use gr_ptParameshInterface, ONLY : gr_ptParticleAtFcBdry
!!$#endif

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

  integer, dimension(MDIM) :: intPos
  real,dimension(MDIM) :: deltas,delInv,hFunc
  real,dimension(LOW:HIGH,MDIM) :: bndBox
  integer :: i,posIndStrt, posIndEnd
  logical :: fcBdry
  real :: particleAttributeValue, dvolInv

  !! Preparatory work for getting integer indices of the cell in 
  !! the particles are 

  call Grid_getDeltas(blockID,deltas)
  call Grid_getBlkBoundBox(blockID,bndBox)

  delInv = 1.0
  delInv(1:NDIM)=1.0/deltas(1:NDIM) !! take the division out of the loop
  posIndStrt = POSX_PART_PROP       !!Start and end points of particles coordinates
  posIndEnd  = posIndStrt+NDIM-1    !! in the particle data structure
  intPos=1
  hfunc = 0.0

  !! Get the inverse of the volume calculated
  dvolInv = delInv(IAXIS)*delInv(JAXIS)*delInv(KAXIS)
  !! DEVNOTE : May be if the block is on fine coarse boundary, and this is the 
  !! finer block, we may have to calculate its impact somewhat differently.  

  do i=1, numParticles

!! commented out becuase i is no longer the particle id
#ifdef DEBUG_GRIDMAPPARTICLES
     if (present(particleOffset)) &
          call pt_validateParticleState(particleOffset+i, blockID)
#endif

     !! Get the integer indices of the cell containing ith particle
     !! using hfunc as the temporary storage for offset from the block bdry
     hfunc(1:NDIM) = (particles(posIndStrt:posIndEnd,i)-&
          bndBox(LOW,1:NDIM))*delInv(1:NDIM)
     intPos(1:NDIM)= floor(hfunc(1:NDIM)) + 1 + guard(1:NDIM)


     !! Now calculate the h function
     hfunc(1:NDIM) = modulo(hfunc(1:NDIM),1.0) - 0.5

     particleAttributeValue=particles(pt_attribute,i)*dvolInv


     !! Decide whether this particle smears across a fine-coarse boundary.
     !! Commented out as we are using the traditional CIC scheme.
     !! DEV : AD -- if we ever need to use this interface, we will
     !!       have to upgrade it to a Unit API
!!$#ifdef FLASH_GRID_PARAMESH
!!$     call gr_ptParticleAtFcBdry(i, fcBdry)
!!$#endif

     !For now, just pass a .false. to calculate mapping using traditional CIC scheme.
     fcBdry = .false.
     call pt_mapOneParticle(blkLimitsGC,intPos,particleAttributeValue, &
          fcBdry,hfunc,buff)
  end do

  return

end subroutine Particles_mapToMeshOneBlk



!Problems caught here are likely to originate in particle movement.
subroutine pt_validateParticleState(particleID, blockID)

  use Particles_data, ONLY : pt_meshMe, particles
  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_outsideBoundBox
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, intent(IN) :: particleID, blockID

  real, dimension(LOW:HIGH,MDIM) :: boundBox
  integer, dimension(MDIM) :: Negh
  real, dimension(MDIM) :: particlePos
  logical :: outside
  integer, parameter :: posIndStrt = POSX_PART_PROP
  integer, parameter :: posIndEnd  = posIndStrt+NDIM-1


  !Check particleID is within the array bounds of particles array.
  !---------------------------------------------------------------------------
  if ( (particleID < lbound(particles,2)) .or. &
       (particleID > ubound(particles,2)) ) then
     print *, "[pt_validateParticleState]: Particle index is outside "//&
          "valid range.  Index ", particleID, &
          ", lower bound ", lbound(particles,2), &
          ", upper bound ", ubound(particles,2)
     call Driver_abortFlash &
          ("[pt_validateParticleState]: Particle index is outside valid range")
  end if



  !Check particle is on the correct block.
  !---------------------------------------------------------------------------
  if (int(particles(BLK_PART_PROP,particleID)) /= blockID) then
     print *, "[pt_validateParticleState]: Particle is on the wrong block"//&
          ". Process ", pt_meshMe, &
          ", particle ID ", particleID, &
          ", tag ", int(particles(TAG_PART_PROP,particleID)), &
          ", intended block ", blockID, &
          ", actual block ", int(particles(BLK_PART_PROP,particleID))
     call Driver_abortFlash &
          ("[pt_validateParticleState]: Particle is on the wrong block")
  end if



  !Check particle position is within block.
  !---------------------------------------------------------------------------
  call Grid_getBlkBoundBox(blockID, boundBox)

  !Note: Grid_outsideBoundBox expects arrays of size MDIM.
  particlePos = 0.0
  particlePos(1:NDIM) = particles(posIndStrt:posIndEnd,particleID)
  call Grid_outsideBoundBox(particlePos, boundBox, outside, Negh)

  if (outside .eqv. .true.) then
     print *, "[pt_validateParticleState]: Particle is outside of block"//&
          ". Process ", pt_meshMe, &
          ", particle ID ", particleID, &
          ", tag ", int(particles(TAG_PART_PROP,particleID)), &
          ", in direction ", negh(1:NDIM), &
          ", particle coords ", particles(posIndStrt:posIndEnd,particleID), &
          ", block coords ", boundBox(LOW:HIGH,1:NDIM)
     call Driver_abortFlash &
          ("[pt_validateParticleState]: Particle is outside of block")
  end if

end subroutine pt_validateParticleState
