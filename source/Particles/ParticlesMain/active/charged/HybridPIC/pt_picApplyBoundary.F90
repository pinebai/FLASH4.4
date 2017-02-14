!!****if* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_picApplyBoundary
!!
!! NAME
!!
!!  pt_picApplyBoundary
!!
!! SYNOPSIS
!!
!!  call pt_picApplyBoundary(logical(in) :: tmp_part)
!!
!! DESCRIPTION
!!   Change position of particles by BCs
!!   should be done after each particle update of positions co-ordinates
!!
!!
!! ARGUMENTS
!!
!!   tmp_part : indicates whether to update intermediate values also
!!
!!
!!
!!***

subroutine pt_picApplyBoundary()

  use Particles_data, only: pt_numLocal, particles
  use pt_picData, ONLY : pt_picDomainBC, pt_picDomainBoundBox
  use Grid_interface, ONLY : Grid_getBlkBoundBox
  use Driver_interface, ONLY : Driver_abortFlash

#include "Flash.h"
#include "constants.h"  
#include "Particles.h"

  implicit none

  real, dimension(MDIM) :: domainSize
  integer :: blockList(MAXBLOCKS), blockCount, i, j
  real, dimension(LOW:HIGH,MDIM) :: boundBox
  real :: b
  
  logical :: allPeriodic, outside

  domainSize(1:MDIM) = pt_picDomainBoundBox(HIGH,1:MDIM) - pt_picDomainBoundBox(LOW,1:MDIM)

  allPeriodic=.true.
  do i = 1,NDIM
     allPeriodic = allPeriodic.and. (pt_picDomainBC(LOW,i)==PERIODIC)
     allPeriodic = allPeriodic.and. (pt_picDomainBC(HIGH,i)==PERIODIC)
  end do

  if(.not.allPeriodic)call Driver_abortFlash("HybridPic only supports periodic boundaries for now")


  ! Check that all particles are in a local block
  ! Assumes Grid_moveParticles has been called after any movement of particles

  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  if (blockCount /= 1) then
     call Driver_abortFlash("[pt_picChkBlk]: Need one local block (UG)")
  end if
  call Grid_getBlkBoundBox(blockList(1), boundBox)

  j = 0
  do i = 1, pt_numLocal
     outside = .false. 
     if (particles(POSX_PART_PROP, i) < boundBox(LOW,IAXIS) .or. &
          particles(POSX_PART_PROP, i) > boundBox(HIGH,IAXIS)) then 
        outside = .true.
     end if
     if (NDIM > 1) then 
        if (particles(POSY_PART_PROP, i) < boundBox(LOW,JAXIS) .or. &
             particles(POSY_PART_PROP, i) > boundBox(HIGH,JAXIS)) then 
           outside = .true.
        end if
     end if
     if (NDIM == 3) then 
        if (particles(POSZ_PART_PROP, i) < boundBox(LOW,KAXIS) .or. &
             particles(POSZ_PART_PROP, i) > boundBox(HIGH,KAXIS)) then 
           outside = .true.
        end if
     end if
     if (outside) then 
        j = j+1
        print *, '  pt_picChkBlk: ', &
          (particles(POSX_PART_PROP:POSX_PART_PROP+2, i)-boundBox(LOW,IAXIS:IAXIS+2)) / &
          (boundBox(HIGH,IAXIS:IAXIS+2)-boundBox(LOW,IAXIS:IAXIS+2) )   ! how far outside?  pos-low/(high-low)
     end if
  end do

  if (j > 0) then
     print *, 'pt_picChkBlk: ', j, ' particles outside block' 
     call Driver_abortFlash("[pt_picApplyBoundary]: Particles outside block")
  end if

end subroutine pt_picApplyBoundary
