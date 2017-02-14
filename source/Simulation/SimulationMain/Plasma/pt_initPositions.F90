!!****if* source/Simulation/SimulationMain/Plasma/pt_initPositions
!!
!! NAME
!!
!!  pt_initPositions
!!
!! SYNOPSIS
!!
!!  call pt_initPositions(integer(in)  :: blockID,
!!                        logical(OUT) :: success)
!!
!! DESCRIPTION
!!
!!    Initialize particle locations.  This version implements a CUSTOM method
!!    for the Plasma PIC test simulation.
!!
!! ARGUMENTS
!!
!!   blockID : ID of block in current processor
!!
!! PARAMETERS
!!
!!
!!***

subroutine pt_initPositions (blockID, success)
  ! Place the initial ions
  use Grid_interface, ONLY : Grid_getBlkType, Grid_getLocalNumBlks, &
    Grid_getBlkPhysicalSize, &
    Grid_getBlkCenterCoords, Grid_getBlkIndexLimits

  use Particles_data, only: pt_numLocal, pt_maxPerProc
  use pt_picData, ONLY : pt_picPmass_1, pt_picPcharge_1,pt_picPdensity_1,&
       pt_picPtemp_1,pt_picPvelx_1, pt_picPvely_1, &
       pt_picPvelz_1,pt_picPpc_1,&
       pt_picPmass_2, pt_picPcharge_2,pt_picPdensity_2,&
       pt_picPtemp_2,pt_picPvelx_2, pt_picPvely_2, &
       pt_picPvelz_2,pt_picPpc_2,&
       pt_picRng_seed
  use pt_picTools, only: pt_picRan_ini, pt_picRan_uni, pt_picMaxwellian, pt_picParticleAdd
!!$  use Particles_data, ONLY : pt_meshMe, pt_meshNumProcs

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(in)  :: blockID
  logical, intent(out) :: success
!!$  logical, save :: first_call = .true.
  integer :: localNumBlocks, blkType
  integer :: blockList(MAXBLOCKS), blockCount, localSize(MDIM)
  integer :: i, j, k, l
  integer :: blkLimits(2,MDIM), blkLimitsGC(2,MDIM), guard(MDIM)
  real :: blockSize(MDIM), blockCenter(MDIM), blockLo(MDIM), h(MDIM)
  real :: r(MDIM), rp(3), x, v(3), dv, pweight

  real, parameter :: amu = 1.66054e-27   ! Atomic mass unit [kg]
  real, parameter :: ec  = 1.6021773e-19 ! Elementary charge [C]

  success = .true.

!!$! Only add particles for LEAF blocks
!!$  call Grid_getBlkType(blockID, blkType) ! Will not link with this
!!$  if (blkType /= LEAF) then 
!!$     return
!!$  end if
  
  print *, 'pt_initPositions: Will put particles in block ', blockID

  call Grid_getLocalNumBlks(localNumBlocks) 

!!$  if (first_call) then
!!$     first_call = .false.
!!$     call pt_picRan_ini(pt_picRng_seed, Pt_meshMe, pt_meshNumProcs)
!!$  end if
  
  ! Get block coordinates
  call Grid_getBlkPhysicalSize(blockID, blockSize)
  call Grid_getBlkCenterCoords(blockID, blockCenter)
  blockLo = blockCenter - 0.5*blockSize ! lower block corner
  ! Get cell indices
  call Grid_getBlkIndexLimits(blockID, &
       blkLimits, blkLimitsGC, CENTER)
  localSize=blkLimits(HIGH,:)-blkLimits(LOW,:)+1   ! NXB, NYB and NZB
  guard = blkLimits(LOW,:)-blkLimitsGC(LOW,:)      ! NGUARD

  h = blockSize/localSize  ! cell size
  dv = product(h(1:NDIM))  ! cell volume [m^3]

  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           ! cell center
           r(1) = blockLo(1) + (i-blkLimits(LOW,IAXIS)+0.5)*h(1) 
           r(2) = blockLo(2) + (j-blkLimits(LOW,JAXIS)+0.5)*h(2) 
           r(3) = blockLo(3) + (k-blkLimits(LOW,KAXIS)+0.5)*h(3) 

           pweight = pt_picPdensity_1*dv/pt_picPpc_1
           do l = 1, pt_picPpc_1
              if (pt_numLocal >= pt_maxPerProc) then
                 success = .FALSE.
                 return
              end if
              ! place particle randomly in cell
              call pt_picRan_uni(rp(1));   call pt_picRan_uni(rp(2)); 
              call pt_picRan_uni(rp(3))
              rp = (rp-0.5)        ! [-0.5,0.5]
              rp = rp*h            ! [-h/2,h/2]
              rp = r + rp
              ! velocity
              call pt_picMaxwellian(v, (/ pt_picPvelx_1,pt_picPvely_1,pt_picPvelz_1 /), &
                   pt_picPtemp_1, pt_picPmass_1*amu)
              call pt_picParticleAdd(rp, v, amu*pt_picPmass_1*pweight, &
                   ec*pt_picPcharge_1*pweight, 1, blockID)
           end do

           pweight = pt_picPdensity_2*dv/pt_picPpc_2
           do l = 1, pt_picPpc_2
              if (pt_numLocal >= pt_maxPerProc) then
                 success = .FALSE.
                 return
              end if
              ! place particle randomly in cell
              call pt_picRan_uni(rp(1));   call pt_picRan_uni(rp(2)); 
              call pt_picRan_uni(rp(3))
              rp = (rp-0.5)        ! [-0.5,0.5]
              rp = rp*h            ! [-h/2,h/2]
              rp = r + rp
              ! velocity
              call pt_picMaxwellian(v, (/ pt_picPvelx_2,pt_picPvely_2,pt_picPvelz_2 /), &
                   pt_picPtemp_2, pt_picPmass_2*amu)
              call pt_picParticleAdd(rp, v, amu*pt_picPmass_2*pweight, &
                   ec*pt_picPcharge_2*pweight, 2, blockID)
           end do
        end do
     end do
  end do

end subroutine pt_initPositions
