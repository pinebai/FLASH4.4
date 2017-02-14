!!****if* source/physics/sourceTerms/Deleptonize/DeleptonizeMain/Deleptonize
!!
!! NAME
!!  
!!  Deleptonize 
!!
!!
!! SYNOPSIS
!! 
!!  call Deleptonize (integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount),
!!             real(IN)    :: dt,
!!             real(IN)    :: time)
!!
!!  
!!  
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***

subroutine Deleptonize (blockCount,blockList,dt,time)
  !
  !==============================================================================
  !
#include "Flash.h"
#include "constants.h"

  use Deleptonize_data, ONLY : useDeleptonize, KtoMeV, MeVtoK, delep_Enu, &
       delep_minDens, delep_useCool, delep_postBounce, delep_bounceTime, &
       delep_rhoTwo, delep_yTwo, delep_useEntr, delep_useRadTrans, &
       delep_threadWithinBlock
  use delep_interface, ONLY : delep_detectBounce
  use Eos_interface, ONLY : Eos_wrapped, Eos_getAbarZbar
  use Eos_nucInterface, ONLY : Eos_nucDetectBounce, Eos_nucOneZone
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits
#ifdef FLASH_MULTISPECIES
  use Multispecies_interface, ONLY : Multispecies_getSumInv, Multispecies_getSumFrac
#include "Multispecies.h"
#endif

  implicit none

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN)::blockList
  real,intent(IN) :: dt,time

  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  integer :: blockID
  integer :: i,j,k,n
  integer :: err
  real :: xDens, xTemp, xEner, xEntr, xYe, outYe
  real :: del_ye, del_entr
  real :: abar, zbar, sumY, Ye0, Ye, dXneut, abarInv
  real :: xPres, mu_nu
  real :: ek
  real :: xXp, xXn, xXa, xXh,xdedt,xdpderho
  real :: tauNu
  integer :: eosMode
  Logical :: threadBlockList
  integer,dimension(blockCount) :: doEOS
  real, parameter :: dydt_max = 20.
  integer :: flag

  call Eos_nucDetectBounce(delep_postBounce,bounceTime=delep_bounceTime)

  if (.not. useDeleptonize) return

  if (delep_postBounce) delep_useEntr = .FALSE.

  if (delep_useEntr) then
     eosMode = MODE_DENS_ENTR
  else
     eosMode = MODE_DENS_EI
  end if

#ifdef FLASH_LEAKAGE
  if (delep_postBounce .AND. delep_useRadTrans) return
#endif

  call Timers_start("delep")

  !$omp parallel if (delep_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(blockID,blkLimits,blkLimitsGC,solnData,flag,blockList,&
  !$omp blockCount,delep_minDens,delep_useEntr,delep_Enu,eosMode,doEos) &
  !$omp private(n,k,j,i,xDens,xTemp,xEner, &
  !$omp xPres,xEntr,xYe,del_ye,outYe,mu_nu,xXp,xXn,xXa,xXh,del_entr,ek,&
  !$omp xdedt,xdpderho) 

  do n = 1, blockCount
     
     !$omp single
     flag = 0
     blockID = blockList(n)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData)
     !$omp end single

     !$omp do schedule(static) reduction(+:flag)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              xDens = solnData(DENS_VAR,i,j,k)
              xTemp = solnData(TEMP_VAR,i,j,k)
              xEner = solnData(EINT_VAR,i,j,k)
              xPres = solnData(PRES_VAR,i,j,k)
              xEntr = solnData(ENTR_VAR,i,j,k)
              xYe   = solnData(YE_MSCALAR,i,j,k)
              
              del_ye = 0.
              del_entr = 0.
              if (xDens <= delep_minDens) then
                 del_ye = 0.
              else
                 call ye_ofRho(xDens,outYe)
                 del_ye = outYe - xYe
                 del_ye = min(0.0, del_ye) ! Deleptonization cannot increase Ye
              endif

              
              if (del_ye /= 0.0) flag = flag+1
              !eosMode = MODE_DENS_ENTR

              if (del_ye /= 0.0 .AND. delep_useEntr) then

                 ! Now call EOS to get chemical potential
                 call Eos_nucOneZone(xDens,xTemp,xYe,xEner,xPres,xEntr,&
                      xdedt,xdpderho,xXp,xXn,xXa,xXh,mu_nu,20,MODE_DENS_TEMP)

                 ! Use chemical potential to find new entropy
                 if (mu_nu < delep_Enu .OR. xDens >= 2.0e12) then
                    del_entr = 0.0
                 else
                    del_entr = -del_ye * (mu_nu - delep_Enu) / (xTemp*KtoMeV) 
                    !eosMode = MODE_DENS_ENTR
                 endif
              endif

              ! Now update entropy, Ye
              xYe = xYe + del_ye
              xEntr = xEntr + del_entr

              solnData(ENTR_VAR,i,j,k)   = xEntr
              solnData(YE_MSCALAR,i,j,k) = xYe 

!!$ Now put species fraction variables into Unk, if they exist
!!$#ifdef YP_VAR
!!$              solnData(YP_VAR,i,j,k) = xXp
!!$#endif
!!$#ifdef YN_VAR
!!$              solnData(YN_VAR,i,j,k) = xXn
!!$#endif
!!$#ifdef YA_VAR
!!$              solnData(YA_VAR,i,j,k) = xXa
!!$#endif
!!$#ifdef YH_VAR
!!$              solnData(YH_VAR,i,j,k) = xXh
!!$#endif        
           enddo
        enddo
     enddo
     !$omp end do

     !$omp single
     call Grid_releaseBlkPtr(blockID,solnData)
     doEos(n) = flag
     !$omp end single
  enddo
  !$omp end parallel
  
  call Timers_stop("delep")
  call Timers_start("Eos")
  do n=1,blockCount
     blockID = blockList(n)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     if (doEos(n) > 0) call Eos_wrapped(eosMode,blkLimits,blockID)
  end do
  call Timers_stop("Eos")

  return
end subroutine Deleptonize


subroutine ye_ofRho(xDens,outYe)
  
  use Deleptonize_data, ONLY : delep_rhoOne, delep_rhoTwo, &
       delep_yOne, delep_yTwo, delep_yc, &
       delep_postBounce, delep_bounceTime

  implicit none
    
  real, intent(IN) :: xDens 
  real, intent(OUT) :: outYe
  
  real :: xofrho, xofrho2

  if (delep_postBounce) then
!!$     delep_yTwo = 0.1
!!$     delep_yTwo = delep_yTwo - 2.0 * dt
!!$     delep_yTwo = max(delep_yTwo, 0.1)
!!$     delep_rhoTwo = 3.e11
  endif

  xofrho = 2.0*log10(xDens) - log10(delep_rhoTwo) - log10(delep_rhoOne)
  xofrho = xofrho / (log10(delep_rhoTwo) - log10(delep_rhoOne))

  xofrho = max(-1.0,min(1.0,xofrho))

!!$  xofrho2 = log10(xDens) - log10(delep_rhoThree)
!!$  xofrho2 = xofrho2 / (log10(delep_rhoThree) - log10(delep_rhoTwo))
!!$  xofrho2 = max(-1.,min(1.,xofrho2))

  outYe = 0.5*(delep_yTwo + delep_yOne) + 0.5*xofrho*(delep_yTwo - delep_yOne) &
          + delep_yc*(1. - abs(xofrho) &
          + 4.*abs(xofrho)*(abs(xofrho)-0.5)*(abs(xofrho) - 1.)) !&
!          + max(0.,(delep_yThree-delep_yTwo)*xofrho2)

  return
end subroutine ye_ofRho
