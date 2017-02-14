!!****if* source/Grid/GridParticles/GridParticlesMove/paramesh/gr_ptMarkRefineDerefine
!!
!! NAME
!!  gr_ptMarkRefineDerefine
!!
!! SYNOPSIS
!!
!!  gr_ptMarkRefineDerefine()
!!
!! DESCRIPTION
!!
!! Routine to mark blocks for refinement or derefinement based upon
!! the counts of particles in them.
!!
!!
!! PARAMETERS
!! 
!!***

subroutine gr_ptMarkRefineDerefine ()

#include "constants.h"
#include "Flash.h"

  use Particles_interface, ONLY : Particles_getCountPerBlk
  use Grid_data, ONLY : gr_refineOnParticleCount,gr_refineOnPdens,&
       gr_minParticlesPerBlk, gr_maxParticlesPerBlk, gr_meshMe
  use gr_ptData, ONLY :   gr_ptMaxPerProc, &
       gr_ptMaxPerProcUpperThresh, gr_ptMaxPerProcLowerThresh, &
       gr_ptMaxPerProcBlockFactor,gr_ptMaxPerProcBlockNoFuzz, &
       gr_ptRefineOnPtMaxPerProc
  use tree, ONLY : refine,derefine,lrefine,lrefine_min,lrefine_max,stay,lnblocks
  use Grid_interface, ONLY : Grid_getListOfBlocks
  implicit none

  integer,dimension(MAXBLOCKS)::oneBlkCount,blkList
  integer :: i,blkCount,blockID, numParticlesThisProc
  integer :: maxParticlesEffective, minParticlesEffective
  integer :: maxLevelWithParticles, minLevelWithParticles
  real,dimension(:,:,:,:),pointer :: solnData
  logical :: deref


  if(.not.gr_refineOnParticleCount)return


  call Grid_getListOfBlocks(LEAF,blkList,blkCount)
  if (blkCount == 0) return

  call Particles_getCountPerBlk(oneBlkCount)
  numParticlesThisProc = sum(oneBlkCount(1:lnblocks))
  if (numParticlesThisProc == 0) return

  
  minLevelWithParticles = 999999
  maxLevelWithParticles = 0
  do i=1,blkCount
     blockID=blkList(i)
     if (oneBlkCount(blockID) > 0) then
        if (lrefine(blockID) > maxLevelWithParticles) &
             maxLevelWithParticles = lrefine(blockID)
        if (lrefine(blockID) < minLevelWithParticles) &
             minLevelWithParticles = lrefine(blockID)
     end if
  end do

  if (gr_ptRefineOnPtMaxPerProc .AND. (numParticlesThisProc > gr_ptMaxPerProcUpperThresh*gr_ptMaxPerProc)) then
     maxParticlesEffective = min( gr_maxParticlesPerBlk, &
          INT(gr_ptMaxPerProc * 1.0 * gr_ptMaxPerProcBlockFactor / (blkCount+gr_ptMaxPerProcBlockNoFuzz)) )
     if (maxParticlesEffective < gr_maxParticlesPerBlk) then
999     format ('PE',I5,' has',I7,' particles on',I3,' leaf blks of levels',I3,'..',I2, &
             ', maxParticlesEffective reduced to',I7,' from',I7,'.')
        print 999,gr_meshMe,numParticlesThisProc,blkCount,minLevelWithParticles,maxLevelWithParticles, &
             maxParticlesEffective,gr_maxParticlesPerBlk
     end if
  else
     maxParticlesEffective = gr_maxParticlesPerBlk
  endif

  ! The logic in following block can usually be ignored. In normal
  ! usage, gr_minParticlesPerBlk (the value of a runtime parameter)
  ! should be 0 or at least very small, so it will typically not need
  ! to be lowered and the else branch will apply.
  if (gr_ptRefineOnPtMaxPerProc .AND. (numParticlesThisProc >= gr_ptMaxPerProcLowerThresh*gr_ptMaxPerProc)) then
     minParticlesEffective = min( gr_minParticlesPerBlk, &
          INT(gr_ptMaxPerProc*2**(3-NDIM) * 0.125 * gr_ptMaxPerProcBlockFactor / (blkCount+gr_ptMaxPerProcBlockNoFuzz)) )
     if (minParticlesEffective < gr_minParticlesPerBlk) then
998     format ('PE',I5,' has',I7,' particles on',I3,' leaf blks of levels',I3,'..',I2, &
             ', minParticlesEffective reduced to',I7,' from',I7,'.')
        print 998,gr_meshMe,numParticlesThisProc,blkCount,minLevelWithParticles,maxLevelWithParticles, &
             minParticlesEffective,gr_minParticlesPerBlk
     end if
  else
     minParticlesEffective = gr_minParticlesPerBlk
  endif

  do i=1,blkCount
     blockID=blkList(i)
     refine(blockID)=(oneBlkCount(blockID)> maxParticlesEffective) &
          .or. refine(blockID)
     refine(blockID)=refine(blockID).and.(lrefine(blockID)<lrefine_max)
     deref=(oneBlkCount(blockID)<minParticlesEffective)
     if (deref) then            ! never changes derefine flags if minParticlesEffective == 0
        derefine(blockID) = (.not. refine(blockID)) .and. &
          (lrefine(blockID)>lrefine_min) &
          .and. (.not. stay(blockID))
     end if
  end do

  return

end subroutine gr_ptMarkRefineDerefine
