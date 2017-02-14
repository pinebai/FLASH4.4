!!****if* source/Grid/GridParticles/GridParticlesMove/UG/Directional/gr_ptMoveOffBlk
!!
!! NAME
!!
!!  gr_ptMoveOffBlk
!!
!! SYNOPSIS
!!
!!  gr_ptMoveOffBlk(real(INOUT)   :: dataBuf(propCount,bufferDim2),
!!                 integer(in)     :: propCount,
!!                 integer(INOUT) :: localCount,
!!                 integer(IN)    :: bufferDim2,
!!                 real(INOUT)    :: destBuf(propCount,bufferDim2),
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
!!     dataBuf -           A 2 dimensional array of elements and their property
!!                           values
!!     propCount - number of element attributes
!!     localCount -   The number of valid elements in the element array that
!!                           are stored on this processor.
!!     bufferDim2 -   The second dimension of the dataBuf and destBuf buffer arrays.
!!     destBuf   -           temporary storage for elements that need to move off processor
!!     numDest   -           number of elements in destBuf
!! 
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_PARTICLES
#endif
subroutine gr_ptMoveOffBlk(dataBuf,propCount,localCount,bufferDim2,destBuf,numDest)

  implicit none
  integer, intent(IN) :: propCount,bufferDim2
  integer,intent(INOUT)::localCount,numDest
  real,dimension(propCount,bufferDim2),intent(INOUT)::dataBuf,destBuf

end subroutine gr_ptMoveOffBlk
