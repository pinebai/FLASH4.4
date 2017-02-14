!!****if* source/physics/sourceTerms/Heat/HeatMain/Neutrino/Heat
!!
!! NAME
!!  
!!  Heat 
!!
!!
!! SYNOPSIS
!! 
!!  call Heat (integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount),
!!             real(IN)    :: dt,
!!             real(IN)    :: time)
!!
!!  
!!  
!! DESCRIPTION
!!
!!   Calculates local heating and cooling due to neutrinos
!!   according to the approach of Murphy & Burrows (2008, ApJ, 688, 1159).
!!
!! NOTES
!!   The FLASH implementation of this method is described in 
!!   S.M. Couch (2013, ApJ, 765, 29).  Citation of this latter 
!!   reference is appreciated if this unit is used in published work.
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***
!
!==============================================================================
!
subroutine Heat(blockCount, blockList, dt, time)

#include "Flash.h"
#include "constants.h"

  use Heat_data, ONLY : useHeat, ht_Lneut, ht_Tneut, ht_bounceTime, ht_postBounce, &
                        ht_useEntr
  use Eos_interface, ONLY : Eos_wrapped, Eos_nucOneZone
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getCellCoords
  use Deleptonize_interface, ONLY : Deleptonize_getBounce

  implicit none
  
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN)::blockList
  real,intent(IN) :: dt,time
  
  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  integer :: blockID
  integer :: i,j,k,n
  integer,dimension(MDIM)  :: dimSize
  real,allocatable, dimension(:) :: xCenter, yCenter, zCenter
  real :: xx, yy, zz
  logical :: gcell = .true.
  real :: radius, dEneut, ek, tauNu
  real, parameter ::  MeVtoK = 1.16045221d10

  real :: bounceTime, temp

  real :: xDens, xTemp, xEner, xEntr, xYe, outYe
  real :: del_ye, del_entr
  real :: abar, zbar, sumY, Ye0, Ye, dXneut, abarInv
  real :: xPres, mu_nu, xCs2
  real :: xXp, xXn, xXa, xXh,xdedt,xdpderho, tmp
  logical :: threadBlockList
  logical :: eosCall2

#ifdef ST_THREAD_BLOCK_LIST
    threadBlockList = .true.

#ifdef ST_THREAD_WITHIN_BLOCK
  call Driver_abortFlash("Cannot include both threading strategies")
#endif

#else
  threadBlockList = .false.
#endif

  if (.NOT. useHeat) return

  if (.NOT. ht_postBounce) &
       call Deleptonize_getBounce(ht_postBounce, ht_bounceTime)

  if (.NOT. ht_postBounce) return

  !$omp parallel if(threadBlockList) &
  !$omp default(none) &
  !$omp private(n,blockID,k,j,i,solnData,dimSize,zCenter,yCenter, &
  !$omp xCenter,radius,dEneut,ek,blkLimits,blkLimitsGC,tmp,xDens,xTemp, &
  !$omp xYe,xEner,xPres,xEntr,xdedt,xdpderho,mu_nu,xXp,xXn,xXa,xXh,xCs2,temp) &
  !$omp shared(blockCount,blockList,ht_Lneut,ht_Tneut,dt,ht_postBounce, &
  !$omp gcell,eosCall2)

  !$omp do schedule(static)
  do n = 1, blockCount
     blockID = blockList(n)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData)
     dimSize(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
     if (NDIM > 2)then
        allocate(zCenter(dimSize(KAXIS)))
        call Grid_getCellCoords(KAXIS,blockID,&
                                 CENTER,gcell,zCenter,dimSize(KAXIS))
     end if
     if (NDIM > 1)then
        allocate(yCenter(dimSize(JAXIS)))
        call Grid_getCellCoords(JAXIS,blockID,&
                                 CENTER,gcell,yCenter,dimSize(JAXIS))
     end if
     
     allocate(xCenter(dimSize(IAXIS)))
     call Grid_getCellCoords(IAXIS,blockID,&
                                 CENTER,gcell,xCenter,dimSize(IAXIS))

#ifdef DELE_VAR
     solnData(DELE_VAR,:,:,:) = 0.0
#endif

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              dEneut = 0.

              xDens = solnData(DENS_VAR,i,j,k)
              xTemp = solnData(TEMP_VAR,i,j,k)
              xYe   = solnData(YE_MSCALAR,i,j,k)
              xEner = solnData(EINT_VAR,i,j,k)

              radius = (xCenter(i))**2
              if (NDIM > 1) then
                 radius = radius + (yCenter(j))**2
              endif
              if (NDIM > 2) then
                 radius = radius + zCenter(k)**2
              endif
              radius = sqrt(radius)
              
              ! Calculate heating
              dEneut = 1.544e20 * (ht_Lneut/1.e52) * (1.e7 / radius)**2 * (ht_Tneut / 4.)**2 

              ! Now subtract cooling 
              dEneut = dEneut - 1.399e20 * (xTemp / (2.*MeVtoK))**6
              
              dEneut = dEneut * exp(-tauNu(xDens))              
              
              ! Now call Eos to get Yp and Yn
              call Eos_nucOneZone(xDens,xTemp,xYe,xEner,xPres,xEntr,&
                   tmp,tmp,tmp,tmp,tmp,tmp,xXp,16,MODE_DENS_TEMP)
              call Eos_nucOneZone(xDens,xTemp,xYe,xEner,xPres,xEntr,&
                   tmp,tmp,tmp,tmp,tmp,tmp,xXn,15,MODE_DENS_TEMP)

              dEneut = dEneut * (xXp + xXn)
              
              solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + dEneut*dt
              
              ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 + &
                          solnData(VELY_VAR,i,j,k)**2 + &
                          solnData(VELZ_VAR,i,j,k)**2)

              solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + ek

! Now store any auxilliary variables
#ifdef DELE_VAR
              solnData(DELE_VAR,i,j,k) = dEneut
#endif
#ifdef TAUN_VAR
              solnData(TAUN_VAR,i,j,k) = tauNu(xDens)
#endif
#ifdef YP_VAR
              solnData(YP_VAR,i,j,k) = xXp
#endif
#ifdef YN_VAR
              solnData(YN_VAR,i,j,k) = xXn
#endif
#ifdef YA_VAR
              solnData(YA_VAR,i,j,k) = xXa
#endif
#ifdef YH_VAR
              solnData(YH_VAR,i,j,k) = xXh
#endif        
           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockID,solndata)
     call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)

     deallocate(xCenter)
     if(NDIM>1)deallocate(yCenter)
     if(NDIM>2)deallocate(zCenter)

  enddo
  !$omp enddo
  !$omp end parallel

  return
end subroutine Heat


function tauNu(dens)

  implicit none

  real, intent(IN) :: dens
  real :: tauNu

  ! Other optical depth approximations may be used here.
  tauNu = dens * 1.0e-11

  return
end function tauNu
