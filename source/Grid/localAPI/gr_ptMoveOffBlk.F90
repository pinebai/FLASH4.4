!!****if* source/Grid/localAPI/gr_ptMoveOffBlk
!!
!! NAME
!!
!!  gr_ptMoveOffBlk
!!
!! SYNOPSIS
!!
!!  gr_ptMoveOffBlk(real(INOUT)   :: particles(part_props,maxLocalNumParticlesOnProc),
!!                 integer(in)    :: part_props,
!!                 integer(INOUT) :: localCount,
!!                 integer(in)    :: maxLocalNumParticlesOnProc,
!!                 real(INOUT)    :: destBuf(part_props,maxLocalNumParticlesOnProc),
!!                 integer(INOUT) :: numDest)
!!
!! DESCRIPTION
!!     
!!    This routine is used in moving the non stationaly data elements 
!!    associated with structures like particles and ray, when a data element
!!    moves off a block without re-gridding. Here every element currently 
!!    on the processor is examined to see if it still belongs to the same block.
!!    If it does not, it is further examimned to see if it has moved out of the physical boundary.
!!    If is out of physical boundary, it may either leave the domain, stay on the same block
!!    or be moved  to destBuf, which holds elements to be passed to the next processor, depending
!!    on the boundary conditions. If it is still in the physical domain, it may have
!!    moved to another block on the same processor, in which case only its BLK
!!    needs to change, otherwise it is moved to destBuf.
!!
!! ARGUMENTS
!!
!!     particles -           A twodimensional array of particles, containing their property
!!                           values
!!     part_props -          number of particle attributes (equal to NPART_PROPS when
!!                           dealing with regular particles of the Particles unit)
!!     localCount -          The number of particles in the particle array that
!!                           are stored on this processor.
!!     maxLocalNumParticlesOnProc - The maximum number of particles in the particle
!!                           array that can be stored on this processor. The size of
!!                           the dummy argument arrays particles and destBuf (second
!!                           dimension).
!!     destBuf   -           temporary storage for particles that need to move off processor
!!     numDest   -           number of particles in destBuf
!! 
!!
!!***
 
subroutine gr_ptMoveOffBlk(particles,part_props,localCount,maxLocalNumParticlesOnProc,destBuf,numDest)
  implicit none
  integer, intent(IN) :: part_props,maxLocalNumParticlesOnProc
  integer,intent(INOUT)::localCount,numDest
  real,dimension(part_props,maxLocalNumParticlesOnProc),intent(INOUT)::particles,destBuf
end subroutine gr_ptMoveOffBlk
