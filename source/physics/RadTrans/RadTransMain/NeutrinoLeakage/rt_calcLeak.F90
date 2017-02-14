!!****if* source/physics/RadTrans/RadTransMain/NeutrinoLeakage/rt_calcLeak
!!
!!  NAME 
!!
!!  rt_calcLeak
!!
!!  SYNOPSIS
!!
!!  call rt_calcLeak( integer(IN) :: nblk,
!!                      integer(IN) :: blklst(nblk),
!!
!!  DESCRIPTION 
!!
!!      Interpolates the taus and luminosities from the leakage rays
!!      back to the hydro grid and then computes local source terms.
!!
!!  ARGUMENTS
!!
!!   nblk   : The number of blocks in the list
!!   blklst : The list of blocks on which the solution must be updated
!!
!!  NOTES
!!      This unit implements ray-by-ray multispecies neutrino leakage.
!!      Parts of this unit are released under a different license than the
!!      usual FLASH license.  Specifically, some subroutines in rt_calcLeak.F90 and 
!!      rt_calcTau.F90 are released under the Creative Commons 
!!      attribution-noncommercial-share alike license.  Basically, if you use this
!!      unit in your work, the license requires that you cite the two articles 
!!      mentioned below.  More details may be found here:  stellarcollapse.org/codes.html.
!!
!!      * O'Connor, E.P., & Ott, C.D. 2010, CQGra, 27, 114103
!!      * Couch, S.M., & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***

subroutine rt_calcLeak(nblk,blklst,dt)

#include "Flash.h"
#include "constants.h"

  use Grid_interface, ONLY : Grid_releaseBlkPtr, Grid_getBlkPtr, &
       Grid_getBlkBoundBox, Grid_getBlkIndexLimits, Grid_getCellCoords, &
       Grid_getSingleCellVol
  use rt_data, ONLY : rt_leakNumRad, rt_leakNumTht, rt_leakRadii,&
       rt_leakX, rt_leakY, rt_leakArr, rt_arraySize, &
       rt_tauRuff, rt_chiRoss, rt_heatFlx, rt_heatErms, rt_heatEave, &
       rt_leakLumTot, rt_leakHeatTot, rt_leakNetTot, rt_leakEaveTot, &
       rt_dTht, temp_mev_to_kelvin, rt_dr, rt_leakTheta, rt_pnsCoord, &
       rt_threadWithinBlock, rt_dPhi, rt_leakPhi
  use Eos_interface, ONLY : Eos_wrapped
  use eosmodule, ONLY : eos_yemin, eos_yemax, e_zeroPoint
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)
  real,    intent(in) :: dt

  real, pointer :: solnData(:,:,:,:)
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer,dimension(MDIM)  :: dimSize
  real,allocatable, dimension(:) :: xCenter, yCenter, zCenter

  integer, parameter :: nfs = 9
  integer, parameter :: nfs2 = 6
  real :: fint(4,nfs),fint3D(8,nfs),fint_out(nfs)
  real :: rr(2),tt(2),pp(2)

  integer :: i,j,k, n
  real :: delX, delY, delZ
  logical :: gcell = .true.

  integer :: blockID
  real :: radius, theta, phi
  real :: qmom, qener, qye

  real  :: xchi(3),xtau(3)
  real  :: xheatflux(3), xlum(3), xeave(3), xheat(3), xnetheat(3)
  real  :: xheaterms(3),xheateave(3)
  real :: dvol
  integer :: indRad=0, indRad2, indTht, indTht2, indPhi, indPhi2

  real :: xDens,xTemp,xYe,xEner,xPres,xEntr,xCs2,xdedt, xGamc
  real :: xdpderho,xdpdrhoe,xmunu, xGame
  integer :: xMode, err

  real, parameter :: Binv = 1.e-51
  real, parameter :: precision = 1.0d-10
  logical :: threadBlockList
  real :: xZone, yZone, zZone
  integer :: flag, doEos(nblk)


  !$omp parallel if(rt_threadWithinBlock) &
  !$omp default(none) &
  !$omp private(n,i,j,k,blockID,dvol,radius,theta,qmom,qener,qye,&
  !$omp         indRad,indTht,indTht2,xlum,xeave,xheat,xnetheat,xtau,xchi,xheatflux,&
  !$omp         xheateave,xheaterms,xtemp,xener,err,xmunu,xdpdrhoe,xdpderho,xcs2,xdedt,&
  !$omp         rr,tt,fint,fint_out,xZone,yZone,zZone,phi,indrad2,indphi,indphi2, &
  !$omp         pp,fint3d) &
#ifndef LEAK_STATIC
  !$omp shared(rt_leakNumTht,rt_leakNumRad) &
#endif
  !$omp shared(nblk,blkLst,rt_leakArr,rt_leakX,rt_dphi,rt_leakPhi, &
  !$omp        rt_leakY,rt_leakRadii,rt_dTht,rt_tauRuff,rt_chiRoss,rt_heatFlx,rt_heateave,&
  !$omp        rt_heaterms,dt,eos_yemin,eos_yemax,e_zeropoint,rt_leaklumtot,rt_leakheattot,&
  !$omp        rt_leaknettot,rt_leakeavetot,gcell,rt_leakTheta,rt_pnsCoord,flag,doEos,&
  !$omp        solnData,blkLimitsGC,blkLimits,dimSize,xCenter,yCenter,zCenter,delX,delY,delZ )

  ! initialize integral quantities
  !$omp workshare
  rt_leakLumTot=0.
  rt_leakHeatTot=0. 
  rt_leakNetTot=0.
  !rt_leakEaveTot=0.
  !$omp end workshare

  do n=1,nblk

     blockID = blklst(n)
     !$omp single
     flag = 0
     call Grid_getBlkPtr(blockID,solnData)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData)
     dimSize(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
     if (NDIM > 2)then
        allocate(zCenter(dimSize(KAXIS)))
        call Grid_getCellCoords(KAXIS,blockID,&
             CENTER,gcell,zCenter,dimSize(KAXIS))
        delZ = zCenter(2) - zCenter(1)
     end if
     if (NDIM > 1)then
        allocate(yCenter(dimSize(JAXIS)))
        call Grid_getCellCoords(JAXIS,blockID,&
             CENTER,gcell,yCenter,dimSize(JAXIS))
        delY = yCenter(2) - yCenter(1)
     end if
     allocate(xCenter(dimSize(IAXIS)))
     call Grid_getCellCoords(IAXIS,blockID,&
          CENTER,gcell,xCenter,dimSize(IAXIS))
     delX = xCenter(2) - xCenter(1)
     !$omp end single
     
     !$omp do schedule(static) &
     !$omp reduction(+:rt_leaklumtot,rt_leakheattot,rt_leaknettot,rt_leakeavetot,flag)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              call Grid_getSingleCellVol(blockID, EXTERIOR, (/i,j,k/), dvol)

              xZone = xCenter(i) - rt_pnsCoord(IAXIS)
              radius = xZone**2
#if NDIM > 1
              yZone = yCenter(j) - rt_pnsCoord(JAXIS)
              radius = radius + yZone**2
              
#if NDIM > 2
              zZone = zCenter(k) - rt_pnsCoord(KAXIS)
              radius = radius + zZone**2
#endif
#endif
              radius = sqrt(radius)
              theta = 0.
              phi = 0.
#if NDIM > 1
              theta = abs(atan2(xZone,yZone))
#if NDIM == 3
              if (zZone >= 0.) then
                 phi = atan2(zZone,xZone)
              else
                 phi = 2.*PI + atan2(zZone,xZone)
              end if
#endif
#endif
              qmom = 0.
              qener = 0.
              qye = 0.

              if (radius <= rt_leakRadii(rt_leakNumRad)) then

                 call ut_hunt(rt_leakRadii, rt_leakNumRad, radius, indRad)
                 indRad = max(indRad,1)
                 indRad2 = indRad+1
                 if (rt_dTht /= 0.) then
                    indTht = int(theta/rt_dTht) + 1
                    indTht2 = indTht + 1
                 else
                    indTht = 1
                    indTht2 = 1
                 end if
                 if (rt_dPhi /= 0.) then
                    indPhi = int(phi/rt_dPhi) + 1
                    indPhi2 = indPhi + 1
                 else
                    indPhi = 1
                    indPhi2 = 1
                 end if

                 xlum(1:3) = 0.0
                 xEave(1:3) = 0.0
                 xHeat(1:3) = 0.0
                 xNetHeat(1:3) = 0.0
#if NDIM==1
                 !! 1D - interpolation 
                 call linterp(rt_leakRadii(indRad),rt_leakRadii(indRad+1),&
                      rt_tauRuff(indRad,1:3,1),rt_tauRuff(indRad+1,1:3,1),&
                      3,radius,xtau(1:3))
                 call linterp(rt_leakRadii(indRad),rt_leakRadii(indRad+1),&
                      rt_chiRoss(indRad,1:3,1),rt_chiRoss(indRad+1,1:3,1),&
                      3,radius,xchi(1:3))
                 call linterp(rt_leakRadii(indRad),rt_leakRadii(indRad+1),&
                      rt_heatFlx(indRad,1:3,1),rt_heatFlx(indRad+1,1:3,1),&
                      3,radius,xheatFlux(1:3))
                 xHeatEave(1:3) = rt_heatEave(1:3,1)
                 xHeatErms(1:3) = rt_heatErms(1:3,1)
#endif
#if NDIM==2
                 !! 2D - interpolation
                 rr(1) = rt_leakRadii(indRad)
                 rr(2) = rt_leakRadii(indRad+1)
                 tt(1) = rt_leakTheta(indTht)
                 tt(2) = rt_leakTheta(indTht+1)

                 fint(1,1:3) = rt_tauRuff(indRad,1:3,indTht)
                 fint(1,4:6) = rt_chiRoss(indRad,1:3,indTht)
                 fint(1,7:9) = rt_heatFlx(indRad,1:3,indTht)

                 fint(2,1:3) = rt_tauRuff(indRad+1,1:3,indTht)
                 fint(2,4:6) = rt_chiRoss(indRad+1,1:3,indTht)
                 fint(2,7:9) = rt_heatFlx(indRad+1,1:3,indTht)

                 fint(3,1:3) = rt_tauRuff(indRad,1:3,indTht+1)
                 fint(3,4:6) = rt_chiRoss(indRad,1:3,indTht+1)
                 fint(3,7:9) = rt_heatFlx(indRad,1:3,indTht+1)

                 fint(4,1:3) = rt_tauRuff(indRad+1,1:3,indTht+1)
                 fint(4,4:6) = rt_chiRoss(indRad+1,1:3,indTht+1)
                 fint(4,7:9) = rt_heatFlx(indRad+1,1:3,indTht+1)

                 call linterp2Dn(rr,tt,fint,9,radius,theta,fint_out)
                 xtau = fint_out(1:3)
                 xchi = fint_out(4:6)
                 xheatFlux = fint_out(7:9)

                 call linterp(tt(1),tt(2),rt_heatEave(1:3,indTht),rt_heatEave(1:3,indTht+1),&
                      3,theta,xHeatEave(1:3))
                 call linterp(tt(1),tt(2),rt_heatErms(1:3,indTht),rt_heatErms(1:3,indTht+1),&
                      3,theta,xHeatErms(1:3))
#endif
#if NDIM==3
                 ! 3D interpolation
                 rr(1) = rt_leakRadii(indRad)
                 rr(2) = rt_leakRadii(indRad+1)
                 tt(1) = rt_leakTheta(indTht)
                 tt(2) = rt_leakTheta(indTht+1)
                 pp(1) = rt_leakPhi(indPhi)
                 pp(2) = rt_leakPhi(indPhi+1)

                 fint3D(1,1:3) = rt_chiRoss(indRad,1:3,rt_leakNumTht*(indPhi-1)+indTht)
                 fint3D(1,4:6) = rt_tauRuff(indRad,1:3,rt_leakNumTht*(indPhi-1)+indTht)
                 fint3D(1,7:9) = rt_heatFlx(indRad,1:3,rt_leakNumTht*(indPhi-1)+indTht)
!
                 fint3D(2,1:3) = rt_chiRoss(indRad,1:3,rt_leakNumTht*(indPhi-1)+indTht2)
                 fint3D(2,4:6) = rt_tauRuff(indRad,1:3,rt_leakNumTht*(indPhi-1)+indTht2)
                 fint3D(2,7:9) = rt_heatFlx(indRad,1:3,rt_leakNumTht*(indPhi-1)+indTht2)
!
                 
                 fint3D(3,1:3) = rt_chiRoss(indRad,1:3,rt_leakNumTht*(indPhi2-1)+indTht)
                 fint3D(3,4:6) = rt_tauRuff(indRad,1:3,rt_leakNumTht*(indPhi2-1)+indTht)
                 fint3D(3,7:9) = rt_heatFlx(indRad,1:3,rt_leakNumTht*(indPhi2-1)+indTht)
!
                 fint3D(4,1:3) = rt_chiRoss(indRad,1:3,rt_leakNumTht*(indPhi2-1)+indTht2)
                 fint3D(4,4:6) = rt_tauRuff(indRad,1:3,rt_leakNumTht*(indPhi2-1)+indTht2)
                 fint3D(4,7:9) = rt_heatFlx(indRad,1:3,rt_leakNumTht*(indPhi2-1)+indTht2)
!
                 fint3D(5,1:3) = rt_chiRoss(indRad2,1:3,rt_leakNumTht*(indPhi-1)+indTht)
                 fint3D(5,4:6) = rt_tauRuff(indRad2,1:3,rt_leakNumTht*(indPhi-1)+indTht)
                 fint3D(5,7:9) = rt_heatFlx(indRad2,1:3,rt_leakNumTht*(indPhi-1)+indTht)
!
                 fint3D(6,1:3) = rt_chiRoss(indRad2,1:3,rt_leakNumTht*(indPhi-1)+indTht2)
                 fint3D(6,4:6) = rt_tauRuff(indRad2,1:3,rt_leakNumTht*(indPhi-1)+indTht2)
                 fint3D(6,7:9) = rt_heatFlx(indRad2,1:3,rt_leakNumTht*(indPhi-1)+indTht2)
!
                 fint3D(7,1:3) = rt_chiRoss(indRad2,1:3,rt_leakNumTht*(indPhi2-1)+indTht)
                 fint3D(7,4:6) = rt_tauRuff(indRad2,1:3,rt_leakNumTht*(indPhi2-1)+indTht)
                 fint3D(7,7:9) = rt_heatFlx(indRad2,1:3,rt_leakNumTht*(indPhi2-1)+indTht)
!
                 fint3D(8,1:3) = rt_chiRoss(indRad2,1:3,rt_leakNumTht*(indPhi2-1)+indTht2)
                 fint3D(8,4:6) = rt_tauRuff(indRad2,1:3,rt_leakNumTht*(indPhi2-1)+indTht2)
                 fint3D(8,7:9) = rt_heatFlx(indRad2,1:3,rt_leakNumTht*(indPhi2-1)+indTht2)

                 call linterp3Dn(tt,pp,rr,fint3D,nfs,theta,phi,radius,fint_out)
                 xchi = fint_out(1:3)
                 xtau = fint_out(4:6)
                 xheatflux = fint_out(7:9)

                 fint(:,:) = 0.0e0

                 fint(1,1:3) = rt_heatEave(1:3,rt_leakNumTht*(indPhi-1)+indTht)
                 fint(1,4:6) = rt_heatErms(1:3,rt_leakNumTht*(indPhi-1)+indTht)

                 fint(2,1:3) = rt_heatEave(1:3,rt_leakNumTht*(indPhi-1)+indTht2)
                 fint(2,4:6) = rt_heatErms(1:3,rt_leakNumTht*(indPhi-1)+indTht2)

                 fint(3,1:3) = rt_heatEave(1:3,rt_leakNumTht*(indPhi2-1)+indTht)
                 fint(3,4:6) = rt_heatErms(1:3,rt_leakNumTht*(indPhi2-1)+indTht)

                 fint(4,1:3) = rt_heatEave(1:3,rt_leakNumTht*(indPhi2-1)+indTht2)
                 fint(4,4:6) = rt_heatErms(1:3,rt_leakNumTht*(indPhi2-1)+indTht2)

                 call linterp2Dn(tt,pp,fint,nfs2,theta,phi,fint_out)
                 xheatEave(1:3) = fint_out(1:3)
                 xheatErms(1:3) = fint_out(4:6)
#endif

                 call calc_leak(solnData(DENS_VAR,i,j,k),solnData(TEMP_VAR,i,j,k)/temp_mev_to_kelvin,&
                      solnData(YE_MSCALAR,i,j,k),xchi,xtau,xheatFlux,xheatErms,&
                      xheatEave,qener,qye,dt,xlum,xEave,xHeat,xnetHeat,-1)


                 if (solnData(YE_MSCALAR,i,j,k)+qye*dt.le.eos_yemin*1.01d0) then
                    !need to surpress any cooling/heating if dyedt.lt.0 near boundary
                    if (qye.lt.0.0e0) then
                       qye = 0.0
                       qener = 0.0
                    endif
                 endif

                 if (solnData(YE_MSCALAR,i,j,k)+qye*dt.ge.eos_yemax*0.99d0) then
                    !need to surpress any cooling/heating if dyedt.lt.0 near boundary
                    if (qye.gt.0.0e0) then
                       qye = 0.0
                       qener = 0.0
                    endif
                 endif
                 
                 !If we are not doing delta-formulation, update solution based on source terms
                 solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + qener*dt
                 solnData(YE_MSCALAR,i,j,k) = solnData(YE_MSCALAR,i,j,k) + qye*dt

                 solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + 0.5*(solnData(VELX_VAR,i,j,k)**2 &
                      + solnData(VELY_VAR,i,j,k)**2 + solnData(VELZ_VAR,i,j,k)**2)

                 xTemp = solnData(TEMP_VAR,i,j,k)/temp_mev_to_kelvin
                 xEner = solnData(EINT_VAR,i,j,k)-e_zeroPoint

                 rt_leakLumTot(1:3) = rt_leakLumTot(1:3) + xlum(1:3)*dvol
                 rt_leakHeatTot(1:3) = rt_leakHeatTot(1:3) + xHeat(1:3)*dvol
                 rt_leakNetTot(1:3) = rt_leakNetTot(1:3) + xnetHeat(1:3)*dvol

                 if (qye /= 0. .OR. qener /= 0.) flag = flag+1
#ifdef DELE_VAR
                 solnData(DELE_VAR,i,j,k) = qener
#endif
#ifdef DYE_VAR
                 solnData(DYE_VAR,i,j,k) = qye
#endif
#ifdef LUM1_VAR
                 solnData(LUM1_VAR,i,j,k) = xlum(1)*Binv
                 solnData(LUM2_VAR,i,j,k) = xlum(2)*Binv
                 solnData(LUM3_VAR,i,j,k) = xlum(3)*Binv
#endif
#ifdef NET1_VAR
                 solnData(NET1_VAR,i,j,k) = xnetHeat(1)*Binv*dvol
                 solnData(NET2_VAR,i,j,k) = xnetHeat(2)*Binv*dvol
                 solnData(NET3_VAR,i,j,k) = xnetHeat(3)*Binv*dvol
#endif
#ifdef HET1_VAR
                 solnData(HET1_VAR,i,j,k) = xHeat(1)*Binv
                 solnData(HET2_VAR,i,j,k) = xHeat(2)*Binv
                 solnData(HET3_VAR,i,j,k) = xHeat(3)*Binv
#endif
#ifdef EAV1_VAR
                 solnData(EAV1_VAR,i,j,k) = xEave(1)
                 solnData(EAV2_VAR,i,j,k) = xEave(2)
                 solnData(EAV3_VAR,i,j,k) = xEave(3)
#endif
#ifdef TAU1_VAR
                 solnData(TAU1_VAR,i,j,k) = xtau(1)
                 solnData(TAU2_VAR,i,j,k) = xtau(2)
                 solnData(TAU3_VAR,i,j,k) = xtau(3)
#endif
              end if !radius < rt_leakRadii(rt_leakNumRad)
           end do ! i
        end do ! j
     end do ! k
     !$omp end do
     !$omp single
     doEos(n) = flag
     call Grid_releaseBlkPtr(blockID,solnData)

     deallocate(xCenter)
#if NDIM >1
     deallocate(yCenter)
#if NDIM >2
     deallocate(zCenter)
#endif
#endif
     !$omp end single
  end do ! blk
  !$omp end parallel
  call Timers_start("EOS")
  do n=1,nblk
     if (doEos(n) > 0) then
        blockID = blklst(n)     
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
     end if
  end do
  call Timers_stop("EOS")
end subroutine rt_calcLeak


  function get_fermi_integral(ifermi,eta)
    implicit none
    integer ifermi
    real*8 get_fermi_integral
    real*8 eta
    real*8 fermi_integral_analytical

    fermi_integral_analytical = 0.0e0

    ! Expressions for Fermi integrals given in Takahashi et al. 1978 
    if (eta.gt.1.D-3) then  
       select case (ifermi)
       case (0)
          fermi_integral_analytical = &
               log10(1.0e0+exp(eta))
       case (1)
          fermi_integral_analytical = &
               (eta**2/2.0E0 + 1.6449d0)/(1.0E0+EXP(-1.6855d0*eta))
       case (2)
          fermi_integral_analytical = &
               (eta**3/3.0E0 + 3.2899d0*eta)/(1.0E0-EXP(-1.8246d0*eta))
       case (3)
          fermi_integral_analytical = & 
               (eta**4/4.0E0 + 4.9348d0*eta**2+11.3644d0) / &
               (1.0E0+EXP(-1.9039d0*eta))        
       case (4)
          fermi_integral_analytical = &
               (eta**5/5.0E0 + 6.5797d0*eta**3+45.4576d0*eta) / &
               (1.0E0-EXP(-1.9484d0*eta))        
       case (5)
          fermi_integral_analytical = &
               (eta**6/6.0E0 + 8.2247d0*eta**4 + 113.6439d0*eta**2 + &
               236.5323d0)/(1.0E0+EXP(-1.9727d0*eta))
       end select

    else
       select case (ifermi)
       case (0)
          fermi_integral_analytical = &
               log10(1.0e0+exp(eta))
       case (1)
          fermi_integral_analytical = &
               EXP(eta)/(1.0E0+0.2159d0*EXP(0.8857d0*eta))
       case (2)
          fermi_integral_analytical = & 
               2.0E0*EXP(eta)/(1.0E0+0.1092d0*EXP(0.8908d0*eta))
       case (3)
          fermi_integral_analytical = & 
               6.0E0*EXP(eta)/(1.0E0+0.0559d0*EXP(0.9069d0*eta))
       case (4)
          fermi_integral_analytical = & 
               24.0E0*EXP(eta)/(1.0E0+0.0287d0*EXP(0.9257d0*eta))
       case (5)
          fermi_integral_analytical = &
               120.0E0*EXP(eta) / (1.0E0 + 0.0147d0*EXP(0.9431d0*eta))
       end select

    endif
    get_fermi_integral  = fermi_integral_analytical

    return
  end function get_fermi_integral 



  subroutine ut_bilinear(x1,x2,y1,y2,q11,q12,q21,q22,xpt,ypt,val)

    implicit none

!!$  integer, intent(in) :: nx, ny
!!$  real, intent(in) :: x(nx), y(ny), z(nx,ny)
    real, intent(in) :: xpt, ypt
    real, intent(out) :: val

!!$  integer,intent(in) :: indLo1, indLo2
    real,intent(in) :: q11, q12, q21, q22
    real,intent(in) :: x1, x2, y1, y2

!!$  indLo1 = 1
!!$  indLo2 = 1
!!$  if (nx>1) call ut_hunt(x,nx,xpt,indLo1)
!!$  if (ny>1) call ut_hunt(y,ny,ypt,indLo2)

    if (x1==x2) then
!!$     y1 = y(indLo2)
!!$     y2 = y(indLo2+1)
!!$     q11 = z(indLo1,indLo2)
!!$     q12 = z(indLo1,indLo2+1)
       val = q11 + (ypt-y1)*(q12-q11)/(y2-y1)
    else if (y1==y2) then
!!$     x1 = x(indLo1)
!!$     x2 = x(indLo1+1)
!!$     q11 = z(indLo1,indLo2)
!!$     q21 = z(indLo1+1,indLo2)
       val = q11 + (xpt-x1)*(q21-q11)/(x2-x1)
    else
!!$     x1 = x(indLo1)
!!$     x2 = x(indLo1+1)
!!$     y1 = y(indLo2)
!!$     y2 = y(indLo2+1)
!!$     q11 = z(indLo1,indLo2)
!!$     q12 = z(indLo1,indLo2+1)
!!$     q21 = z(indLo1+1,indLo2)
!!$     q22 = z(indLo1+1,indLo2+1)

       val = q11/(x2-x1)/(y2-y1)*(x2-xpt)*(y2-ypt) &
            +q21/(x2-x1)/(y2-y1)*(xpt-x1)*(y2-ypt) &
            +q12/(x2-x1)/(y2-y1)*(x2-xpt)*(ypt-y1) &
            +q22/(x2-x1)/(y2-y1)*(xpt-x1)*(ypt-y1)
    end if

    return
  end subroutine ut_bilinear



  subroutine linterp2Dn(xin,yin,fin,n,x,y,fout)

    implicit none
    real  :: xin(2),yin(2)
    integer :: n
    real  :: fin(4,n)
    real  :: x,y,fout(n)
    real  :: invdxdy

    !      y2 f3           f4
    !
    !
    !      y1 f1           f2
    !         x1           x2

    invdxdy = 1.0e0/( (xin(2)-xin(1))*(yin(2)-yin(1)) )

    fout(1:n) = invdxdy * (  fin(1,1:n) * (xin(2)-x)*(yin(2)-y) &
         + fin(2,1:n) * (x-xin(1))*(yin(2)-y) &
         + fin(3,1:n) * (xin(2)-x)*(y-yin(1)) &
         + fin(4,1:n) * (x-xin(1))*(y-yin(1)) )


  end subroutine linterp2Dn

  subroutine linterp3Dn(xin,yin,zin,fin,n,x,y,z,fout)

    implicit none
    real  :: xin(2),yin(2),zin(2)
    integer :: n
    real  :: fin(8,n)
    real  :: x,y,z,fout(n)
    real  :: invdxdydz

    !        y2 f7           f8
    !
    !  z2
    !        y1 f5           f6
    !           x1           x2

    !        y2 f3           f4
    !
    !  z1
    !        y1 f1           f2
    !           x1           x2

    invdxdydz = 1.0e0/( (xin(2)-xin(1))*(yin(2)-yin(1))*(zin(2)-zin(1)) )

    ! auto-generated via mathematica, terms slightly reordered:
    fout(1:n) = invdxdydz * (                                 &
         fin(8,1:n)*(x - xin(1))*(y - yin(1))*(z - zin(1)) +  &
         fin(7,1:n)*(xin(2) - x)*(y - yin(1))*(z - zin(1)) +  &
         fin(5,1:n)*(x - xin(2))*(y - yin(2))*(z - zin(1)) +  &
         fin(6,1:n)*(x - xin(1))*(yin(2) - y)*(z - zin(1)) +  &
         fin(4,1:n)*(x - xin(1))*(y - yin(1))*(zin(2) - z) +  &
         fin(3,1:n)*(xin(2) - x)*(y - yin(1))*(zin(2) - z) +  &
         fin(2,1:n)*(x - xin(1))*(yin(2) - y)*(zin(2) - z) +  &
         fin(1,1:n)*(xin(2) - x)*(yin(2) - y)*(zin(2) - z) ) 

  end subroutine linterp3Dn

  subroutine linterp(x1,x2,y1,y2,n,x,y)

    implicit none
    integer :: n
    real :: slope(n),x1,x2,y1(n),y2(n),x,y(n)

    if (x2.lt.x1) then
!!$     call CCTK_WARN(0,"Error in linterp!")
    endif

    slope(1:n) = (y2(1:n) - y1(1:n)) / (x2 - x1)
    y(1:n) = slope(1:n)*(x-x1) + y1(1:n)

  end subroutine  linterp

  subroutine linterp3(x1,x2,y1,y2,x,y)

    implicit none

    real :: slope(3),x1(3),x2(3),y1(3),y2(3),x,y(3)

    slope = (y2 - y1) / (x2 - x1)
    y = slope*(x-x1) + y1

  end subroutine  linterp3

subroutine calc_leak(rho,temp,ye,chi,tau,heatflux,heaterms,heateave,&
     depsdt,dyedt,ldt,lum,eave,heatout,netheatout,reflevel)
! WARNING: Be careful when changing the arguments to this function; it is also
! called from calc_tau to get the luminosity available for heating
! along the rays.


  use eosmodule, only : energy_shift, eos_yemin, eos_yemax
  use rt_data, ONLY : rt_leakDoHeat, rt_leakHeatFac
!  use rt_utils, ONLY : get_fermi_integral

  implicit none
!!  DECLARE_CCTK_PARAMETERS

  real, intent(in) :: rho ! density in g/cm^3
  real, intent(in) :: temp ! temperature in MeV
  real, intent(in) :: ye ! ye, dimensionless

  real, intent(in) :: chi(3) !chi calculated from Rosswog scheme, one for each nu, to be interpolated throughout 3D
  real, intent(in) :: tau(3) !tau calculated from Rosswog scheme, one for each nu, to be interpolated throughout 3D
  real, intent(in) :: heatflux(3) !flux used for heating, one for each nu, interpolated throughout 3D
  real, intent(in) :: heaterms(3) !rms neutrino energy at neutrinosphere, one for each nu, interpolated throughtout 3D
  real, intent(in) :: heateave(3) !average neutrino energy at neutrinosphere, one for each nu, interpolated throughtout 3D
  real, intent(out) :: lum(3) ! local dE/dt (luminosities)
  real, intent(out) :: eave(3) ! local <E> (average energy)
  real, intent(out) :: heatout(3) ! local Q+ ergs/s/cm^3 (heating)
  real, intent(out) :: netheatout(3) ! local net heating ergs/s/cm^3 (heating)
  real, intent(out) :: depsdt !change in the internal energy, ergs/cm^3/s ... Actually I think this is in ergs/g/s
  real, intent(out) :: dyedt !change in electron fraction
  integer, intent(in) :: reflevel

  !EOS & local variables
  integer :: keytemp, keyerr
  real :: precision = 1.0d-10
  real :: matter_rho,matter_temperature,matter_ye
  real :: matter_enr,matter_prs,matter_ent
  real :: matter_cs2,matter_dedt,matter_dpdrhoe
  real :: matter_dpderho,matter_xa,matter_xh
  real :: matter_xn,matter_xp,matter_abar
  real :: matter_zbar,matter_mue,matter_mun
  real :: matter_mup,matter_muhat

  real :: eta_e,eta_p,eta_n
  real :: eta_hat,eta_nue,eta_nua,eta_nux
  real :: eta_pn, eta_np


  !for heating
  real :: heat_const !preamble
  ! real :: f_heat !factor for increasing heating -- set via Cactus parameter
  real :: F(2) !Inverse flux factor
  real :: heat_eff(2) !effective heating

  !constants & parameters
  real, parameter :: Qnp = 1.293333d0 !m_n - m_p, MeV
  real, parameter :: Cv = 0.5d0 + 2.0e0*0.23d0 !vector coupling
  real, parameter :: Ca = 0.5d0 !axial coupling
  real, parameter :: alpha = 1.23d0 !gA
  real, parameter :: me_mev = 0.510998910e0 !mass of electron in MeV
  real, parameter :: sigma_0 = 1.76d-44 ! in units of cm^2
  real, parameter :: avo = 6.0221367d23 !Avogadro's number
  real, parameter :: pi = 3.141592653589793238462d0 !pi
  real, parameter :: clite = 29979245800.0e0 !speed of light
  real, parameter :: hc_mevcm = 1.23984172d-10 !hc in units of MeV*cm
  real, parameter :: gamma_0 = 5.565d-2 ! dimensionless
  real, parameter :: fsc = 1.0e0/137.036d0 ! fine structure constant, dimensionless
  real, parameter :: mev_to_erg = 1.60217733d-6 !conversion
  real, parameter :: massn_cgs = 1.674927211d-24 !neutron mass in grams
  
  !diffusion
  real :: scattering_kappa,abs_kappa
  real :: kappa_tilde_nu_scat(3,3)
  real :: kappa_tilde_nu_abs(3,3)
  real :: block_factor
  real :: zeta(3)
  real :: rate_const
  real :: R_diff(3),Q_diff(3)

  !free emission
  real :: beta
  real :: gamma_const,gamma, R_gamma
  real :: block_factor_e,block_factor_a,block_factor_x
  real :: R_pair, Q_pair, pair_const
  real :: enr_m,enr_p,enr_tilde_m,enr_tilde_p
  real :: R_loc(3),Q_loc(3)

  !leakage
  real :: R_eff(3),Q_eff(3)

  !function
  real :: get_fermi_integral

  real :: ldt

  !cactus relates stuff
  character(len=512) :: warnline

  matter_rho = rho
  matter_temperature = temp
  matter_ye = ye

  keytemp = 1
  keyerr = 0
  call nuc_eos_full(matter_rho,matter_temperature,matter_ye,matter_enr, &
       matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho, &
       matter_dpdrhoe,matter_xa,matter_xh,matter_xn,matter_xp,matter_abar, &
       matter_zbar,matter_mue,matter_mun,matter_mup,matter_muhat, &
       keytemp,keyerr,precision)
  if (keyerr.ne.0) then
!!$     !$OMP CRITICAL(leakleak1)
!!$     write(warnline,"(A15,1P10E15.6)") "rho: ", matter_rho
!!$     call CCTK_WARN(1,warnline)
!!$     write(warnline,"(A15,1P10E15.6)") "temperature: ", matter_temperature
!!$     call CCTK_WARN(1,warnline)
!!$     write(warnline,"(A15,1P10E15.6)") "ye: ", matter_ye
!!$     call CCTK_WARN(1,warnline)
!!$     write(warnline,"(A15,i10)") "eos error", keyerr
!!$     call CCTK_WARN(1,warnline)
!!$     call CCTK_WARN(0,"set_eos_variables: EOS error in leakage calc_leak")
!!$     !$OMP END CRITICAL(leakleak1)
  endif

  ! don't do anything outside the shock
  if(matter_xh.gt.0.5.and.matter_rho.lt.1.0d13) then 
     lum(1:3) = 0.0e0
     depsdt = 0.0e0
     dyedt = 0.0e0
     eave = 0.0e0
     return
  endif

  eta_e = matter_mue/temp
  eta_p = matter_mup/temp
  eta_n = matter_mun/temp

  eta_hat = eta_n - eta_p - Qnp/temp
  eta_nue = eta_e - eta_n + eta_p !fully includes effects of rest masses
  eta_nua = -eta_nue
  eta_nux = 0.0e0

  !interpolate etas like we do in GR1D
  eta_nue = eta_nue*(1.0e0-exp(-tau(1)))
  eta_nua = eta_nua*(1.0e0-exp(-tau(2)))
  
  scattering_kappa = rho*avo*0.25d0*sigma_0/me_mev**2
  kappa_tilde_nu_scat(1,1) = matter_xn*scattering_kappa
  kappa_tilde_nu_scat(1,2) = matter_xp*scattering_kappa
  kappa_tilde_nu_scat(2,1) = matter_xn*scattering_kappa
  kappa_tilde_nu_scat(2,2) = matter_xp*scattering_kappa
  kappa_tilde_nu_scat(3,1) = matter_xn*scattering_kappa
  kappa_tilde_nu_scat(3,2) = matter_xp*scattering_kappa
  
  scattering_kappa = rho*avo*0.0625d0*sigma_0/me_mev**2* &
       matter_abar*(1.0e0-matter_zbar/matter_abar)**2

  kappa_tilde_nu_scat(1,3) = matter_xh*scattering_kappa
  kappa_tilde_nu_scat(2,3) = matter_xh*scattering_kappa
  kappa_tilde_nu_scat(3,3) = matter_xh*scattering_kappa

  eta_pn = avo*rho*(matter_xn-matter_xp)/(exp(eta_hat)-1.0e0)
  eta_pn = max(eta_pn,0.0e0)
  eta_np = avo*rho*(matter_xp-matter_xn)/(exp(-eta_hat)-1.0e0)
  eta_np = max(eta_np,0.0e0)

  if (rho.lt.1.0d11) then
     !non degenerate here, use mass fractions as chemical potentials fail at low densities
     eta_pn = avo*rho*matter_xp
     eta_np = avo*rho*matter_xn
  endif

  !absorption
  abs_kappa = (1.0e0+3.0e0*alpha**2)*0.25d0*sigma_0/me_mev**2
  block_factor = 1.0e0 + exp(eta_e-get_fermi_integral(5,eta_nue)/ &
       get_fermi_integral(4,eta_nue))
  kappa_tilde_nu_abs(1,1) = eta_np*abs_kappa/block_factor
  kappa_tilde_nu_abs(2,1) = 0.0e0 !no absorption of a-type on neutrons
  kappa_tilde_nu_abs(3,1) = 0.0e0 !no absorption of x-type neutrinos
  kappa_tilde_nu_abs(1,2) = 0.0e0 !no absorption of e-type on protons
  block_factor = 1.0e0 + exp(-eta_e-get_fermi_integral(5,eta_nua)/ &
       get_fermi_integral(4,eta_nua))
  kappa_tilde_nu_abs(2,2) = eta_pn*abs_kappa/block_factor
  kappa_tilde_nu_abs(3,2) = 0.0e0 !no absorption of x-type neutrinos
  kappa_tilde_nu_abs(1,3) = 0.0e0 !no absorption on nuclei
  kappa_tilde_nu_abs(2,3) = 0.0e0 !no absorption on nuclei
  kappa_tilde_nu_abs(3,3) = 0.0e0 !no absorption on nuclei

  !sum up opacities to get zeta (again, factoring out energy dependence)
  zeta(1) = kappa_tilde_nu_scat(1,1) + kappa_tilde_nu_scat(1,2) + &
       kappa_tilde_nu_scat(1,3) + kappa_tilde_nu_abs(1,1) + &
       kappa_tilde_nu_abs(1,2) + kappa_tilde_nu_abs(1,3)

  zeta(2) = kappa_tilde_nu_scat(2,1) + kappa_tilde_nu_scat(2,2) + &
       kappa_tilde_nu_scat(2,3) + kappa_tilde_nu_abs(2,1) + &
       kappa_tilde_nu_abs(2,2) + kappa_tilde_nu_abs(2,3)

  zeta(3) = kappa_tilde_nu_scat(3,1) + kappa_tilde_nu_scat(3,2) + &
       kappa_tilde_nu_scat(3,3) + kappa_tilde_nu_abs(3,1) + &
       kappa_tilde_nu_abs(3,2) + kappa_tilde_nu_abs(3,3)

  rate_const = 4.0e0*pi*clite*zeta(1)/(hc_mevcm**3*6.0e0*chi(1)**2)
  R_diff(1) = rate_const*temp*get_fermi_integral(0,eta_nue)
  Q_diff(1) = rate_const*temp**2*get_fermi_integral(1,eta_nue)

  rate_const = 4.0e0*pi*clite*zeta(2)/(hc_mevcm**3*6.0e0*chi(2)**2)
  R_diff(2) = rate_const*temp*get_fermi_integral(0,eta_nua)
  Q_diff(2) = rate_const*temp**2*get_fermi_integral(1,eta_nua)
       
  rate_const = 16.0e0*pi*clite*zeta(3)/(hc_mevcm**3*6.0e0*chi(3)**2)
  R_diff(3) = rate_const*temp*get_fermi_integral(0,eta_nux)
  Q_diff(3) = rate_const*temp**2*get_fermi_integral(1,eta_nux)

  !now for the free emission
  beta = pi*clite*(1.0e0+3.0e0*alpha**2)*sigma_0/(hc_mevcm**3*me_mev**2)
  
  R_loc = 0.0e0
  Q_loc = 0.0e0

  R_loc(1) = beta*eta_pn*temp**5*get_fermi_integral(4,eta_e)
  Q_loc(1) = beta*eta_pn*temp**6*get_fermi_integral(5,eta_e)
  R_loc(2) = beta*eta_np*temp**5*get_fermi_integral(4,-eta_e)
  Q_loc(2) = beta*eta_np*temp**6*get_fermi_integral(5,-eta_e)

  !e-e+ pair processes from Ruffert et al.
  block_factor_e = 1.0e0+exp(eta_nue-0.5d0*( &
       get_fermi_integral(4,eta_e)/get_fermi_integral(3,eta_e) + &
       get_fermi_integral(4,-eta_e)/get_fermi_integral(3,-eta_e) &
       ))
  block_factor_a = 1.0e0+exp(eta_nua-0.5d0*( &
       get_fermi_integral(4,eta_e)/get_fermi_integral(3,eta_e) + &
       get_fermi_integral(4,-eta_e)/get_fermi_integral(3,-eta_e) &
       ))
  block_factor_x = 1.0e0+exp(eta_nux-0.5d0*( &
       get_fermi_integral(4,eta_e)/get_fermi_integral(3,eta_e) + &
       get_fermi_integral(4,-eta_e)/get_fermi_integral(3,-eta_e) &
       ))

  enr_m = 8.0e0*pi/hc_mevcm**3*temp**4*get_fermi_integral(3,eta_e)
  enr_p = 8.0e0*pi/hc_mevcm**3*temp**4*get_fermi_integral(3,-eta_e)
  
  enr_tilde_m = 8.0e0*pi/hc_mevcm**3*temp**5*get_fermi_integral(4,eta_e)
  enr_tilde_p = 8.0e0*pi/hc_mevcm**3*temp**5*get_fermi_integral(4,-eta_e)
  
  pair_const = sigma_0*clite/me_mev**2*enr_m*enr_p
  
  R_pair =  pair_const/(36.0e0*block_factor_e*block_factor_a)* &
       ((Cv-Ca)**2+(Cv+Ca)**2)
  R_loc(1) = R_loc(1) + R_pair
  Q_loc(1) = Q_loc(1) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)
  R_loc(2) = R_loc(2) + R_pair
  Q_loc(2) = Q_loc(2) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)
  
  R_pair =  pair_const/(9.0e0*block_factor_x**2)*((Cv-Ca)**2+(Cv+Ca-2.0e0)**2)
  R_loc(3) = R_loc(3) + R_pair
  Q_loc(3) = Q_loc(3) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)
  
  !plasmon decay from Ruffert et al.
  gamma = gamma_0*sqrt((pi**2+3.0e0*eta_e**2)/3.0e0)
  block_factor_e = 1.0e0 + exp(eta_nue-(1.0e0+0.5d0*gamma**2/(1.0e0+gamma)))
  block_factor_a = 1.0e0 + exp(eta_nua-(1.0e0+0.5d0*gamma**2/(1.0e0+gamma)))
  block_factor_x = 1.0e0 + exp(eta_nux-(1.0e0+0.5d0*gamma**2/(1.0e0+gamma)))
  
  gamma_const = pi**3*sigma_0*clite*temp**8/(me_mev**2*3.0e0*fsc*hc_mevcm**6)* &
       gamma**6*exp(-gamma)*(1.0e0+gamma)
  
  R_gamma = Cv**2*gamma_const/(block_factor_e*block_factor_a)
  R_loc(1) = R_loc(1) + R_gamma
  Q_loc(1) = Q_loc(1) + R_gamma*0.5d0*temp*(2.0e0+gamma**2/(1.0e0+gamma))
  R_loc(2) = R_loc(2) + R_gamma
  Q_loc(2) = Q_loc(2) + R_gamma*0.5d0*temp*(2.0e0+gamma**2/(1.0e0+gamma))

  R_gamma = (Cv-1.0e0)**2*4.0e0*gamma_const/block_factor_x**2
  R_loc(3) = R_loc(3) + R_gamma 
  Q_loc(3) = Q_loc(3) + R_gamma*0.5d0*temp*(2.0e0+gamma**2/(1.0e0+gamma))

  R_pair = 0.231d0*(2.0778d2/mev_to_erg)*0.5d0* &
       (matter_xn**2+matter_xp**2+28.0e0/3.0e0*matter_xn*matter_xp)* &
       rho**2*temp**(4.5d0)
  Q_pair = R_pair*temp/0.231d0*0.504d0
          
  R_loc(1) = R_loc(1) + R_pair
  Q_loc(1) = Q_loc(1) + Q_pair
          
  R_loc(2) = R_loc(2) + R_pair
  Q_loc(2) = Q_loc(2) + Q_pair
  R_loc(3) = R_loc(3) + 4.0e0*R_pair
  Q_loc(3) = Q_loc(3) + 4.0e0*Q_pair

  !now calculate leakage!!!!!!!
  R_eff(:) = R_loc(:)/(1.0e0+R_loc(:)/R_diff(:))
  Q_eff(:) = Q_loc(:)/(1.0e0+Q_loc(:)/Q_diff(:))

  eave(1:3) = Q_eff(1:3)/R_eff(1:3)

  if (rt_leakDoHeat) then
     
     heat_const = rt_leakHeatFac*(1.0e0+3.0e0*alpha**2) * sigma_0 / (me_mev**2 * massn_cgs) * 0.25d0 !cm^2/MeV^2/g

     F(1:2) = (4.275d0*tau(1:2)+1.15d0)*exp(-2.0e0*tau(1:2))
    
     heat_eff(1) = heat_const * rho * matter_xn * heatflux(1) * &
          heaterms(1)**2 * F(1) / mev_to_erg ! cm^2/MeV^2/g * g/cm^3 erg/s/cm^2 MeV^2 MeV/erg = MeV/cm^3/s
     heat_eff(2) = heat_const * rho * matter_xp * heatflux(2) * &
          heaterms(2)**2 * F(2) / mev_to_erg ! cm^2/MeV^2/g * g/cm^3 erg/s/cm^2 MeV^2 MeV/erg = MeV/cm^3/s

     lum(1:2) = (Q_eff(1:2)-heat_eff(1:2))*mev_to_erg !ergs/cm^3/s
! debug:
!     lum(1:2) = (Q_eff(1:2))*mev_to_erg !ergs/cm^3/s
     lum(3) = Q_eff(3)*mev_to_erg !ergs/cm^3/s

     heatout(1:2) = heat_eff(1:2)*mev_to_erg !ergs/cm^3/s
     heatout(3) = 0.0e0

     netheatout(1) = abs(min(0.0e0,lum(1)))
     netheatout(2) = abs(min(0.0e0,lum(2)))
     netheatout(3) = 0.0e0
     
     depsdt = -sum(lum(1:3))/rho
     dyedt = -(R_eff(1)-R_eff(2)+heat_eff(2)/heateave(2)-heat_eff(1)/heateave(1))*massn_cgs/rho
! debug:
!     dyedt = -(R_eff(1)-R_eff(2))*massn_cgs/rho

  else
     lum(1:3) = Q_eff(1:3)*mev_to_erg
     depsdt = -sum(Q_eff(1:3))*mev_to_erg/rho
     dyedt = -(R_eff(1)-R_eff(2))*massn_cgs/rho
  endif

  if (ye+dyedt*ldt.le.eos_yemin*1.05d0) then
     dyedt = 0.0e0
     depsdt = 0.0e0
  endif

  if (ye+dyedt*ldt.ge.eos_yemax*0.95d0) then
     dyedt = 0.0e0
     depsdt = 0.0e0
  endif


!!$  if (matter_enr+depsdt*ldt.lt.-energy_shift.and.&
!!$       (reflevel.eq.-1.or.reflevel.ge.grhydro_c2p_warn_from_reflevel)) then
!!$     !$OMP CRITICAL(leak667)
!!$     call CCTK_WARN(1,"Problem in leakage; energy change too large")
!!$     write(warnline,"(A15,i10)") "reflevel: ", reflevel
!!$     call CCTK_WARN(1,warnline)
!!$     write(warnline,"(A15,1P10E15.6)") "rho: ", rho
!!$     call CCTK_WARN(1,warnline)
!!$     write(warnline,"(A15,1P10E15.6)") "temp: ", temp
!!$     call CCTK_WARN(1,warnline)
!!$     write(warnline,"(A15,1P10E15.6)") "Y_e: ", ye
!!$     call CCTK_WARN(1,warnline)
!!$     write(warnline,"(A15,1P10E15.6)") "eps: ", matter_enr
!!$     call CCTK_WARN(1,warnline)
!!$     write(warnline,"(A15,1P10E15.6)") "depsdt*ldt:", depsdt*ldt
!!$     call CCTK_WARN(1,warnline)
!!$     write(warnline,"(A30,1P10E15.6)") "matter_enr+depsdt*ldt:", matter_enr+depsdt*ldt
!!$     call CCTK_WARN(1,warnline)
!!$     call CCTK_WARN(0,"aborting")
!!$     !$OMP END CRITICAL(leak667)
!!$  endif
  

end subroutine calc_leak

