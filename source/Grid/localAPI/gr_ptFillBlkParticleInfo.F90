!!****if* source/Grid/localAPI/gr_ptFillBlkParticleInfo
!!
!! NAME
!!  gr_ptFillBlkParticleInfo
!!
!! SYNOPSIS
!!
!!  call gr_ptFillBlkParticleInfo()
!!                    
!!  
!! DESCRIPTION 
!!
!!  This is an old interface for an obsolete GridParticles local API routine.
!!  No real implementation exists in the current FLASH code.  The routine
!!  was used as part off the interface to an old implementation of
!!  particle movement.
!!
!!  Calls to this routine should be removed; existing declarations of
!!  array gr_blkParticleInfo in various implementations of Grid_data
!!  should also be removed.
!!
!!  Obsolete documentation follows.
!!  
!!    gr_blkParticleInfo is a multidimensional array holding various info
!!    needed for updating the particle refinement
!!    gr_blkParticleInfo(1,:) holds globalID for 1st dim
!!    gr_blkParticleInfo(2,:) holds globalID for 2nd dim
!!    gr_blkParticleInfo(3,:) holds globalID for 3rd dim
!!    gr_blkParticleInfo(4,:) keeps track of if an old block has been processed or not
!!
!!  Initialize the gr_blkParticleInfo array with a globalID for each block.  
!!  The first 3 dimensions of the multi-dimensional array gr_blkParticleInfo 
!!  stores the a unique
!!  identifier for each block before refinement occurs.  This globalID
!!  is used to keep a history of the local blocks before and after
!!  grid refinement.  It is made up of a combination of the cornerID,
!!  refinement level and nodetype.
!!
!!  Dimension 4 keeps track of whether an old block has been fully processed or not.
!!  An old block is considered fully processed when
!!  its entry in array is 2**NDIM.  If block is refined
!!  then incremented by 1 for each child
!!  needs to be dimensioned as oldLocalNumBlocks, rather
!!  than oldLocalNumLeafBlocks because blockNum is used to
!!  identify the order
!!
!!
!!
!!
!! ARGUMENTS 
!!  none
!!  
!!
!! NOTES
!!   
!! SEE ALSO
!!  Grid_getBlkCornerID
!!
!!***



subroutine gr_ptFillBlkParticleInfo()
 
  implicit none

end subroutine gr_ptFillBlkParticleInfo
