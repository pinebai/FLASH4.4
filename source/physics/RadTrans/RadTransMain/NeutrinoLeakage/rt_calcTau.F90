!!****if* source/physics/RadTrans/RadTransMain/NeutrinoLeakage/rt_calcTau
!!
!!  NAME 
!!
!!  rt_calcTau
!!
!!  SYNOPSIS
!!
!!  call rt_calcTau
!!
!!  DESCRIPTION 
!!      This routine computes neutrino optical depths using the Rosswog
!!      and Ruffert approaches.  Kernel implementations are from GR1D.
!!      Fancy MPI jiu-jitsu is used for efficient communication.
!!
!!  ARGUMENTS
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
!!
!!***

subroutine rt_calcTau

#include "Flash.h"

  use rt_data, ONLY : rt_leakNumRad, rt_leakNumTht, rt_globalMe, &
       rt_leakArr, rt_leakRadii, rt_dr, temp_mev_to_kelvin, &
       rt_tauRuff, rt_chiRoss, rt_heatFlx, rt_heatErms, rt_heatEave,&
       rt_subMeshComm, rt_arraySize, rt_arraySmall, &
       rt_threadWithinBlock, rt_arraySizeLocal, rt_istart, rt_iend, rt_rayData, &
       rt_recvCnt, rt_dsplCnt, rt_leakEaveTot, rt_leakTheta, rt_radNu, rt_leakNumRays
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

  integer :: m,n
  real :: oldtau(rt_leakNumRad,3)
  real :: compos(rt_leakNumRad,4)
  real, dimension(rt_leakNumRad) :: xrho, xtemp, xye
  logical :: threadBlockList
  integer :: error
  real :: ns_rad(rt_leakNumRays)
  real, parameter :: twoThirds=2./3.

  include "Flash_mpi.h"


  call Timers_start("calc")
  !$omp parallel if(rt_threadWithinBlock) &
  !$omp default(none) &
  !$omp private(n,xrho,xtemp,xye,oldtau,compos) &
#ifndef LEAK_STATIC
  !$omp shared(rt_leakNumRays,rt_leakNumRad) &
#endif
  !$omp shared(rt_leakArr,rt_rayData,&
  !$omp        rt_leakRadii,rt_dr,rt_istart,rt_iend,&
  !$omp        rt_tauRuff,rt_chiRoss,rt_heatFlx,rt_heatErms,rt_heatEave) 

  !$omp do schedule(static)   
  do n=rt_istart, rt_iend
     ! input the old taus
     oldtau(:,1:3) = rt_tauRuff(:,1:3,n)
     xTemp = rt_leakArr(2,:,n)/temp_mev_to_kelvin
     ! Now calculate tau, chi, and heatFlux for the ray
     call calc_taus(rt_leakArr(1,:,n),xTemp,rt_leakArr(3,:,n),&
          oldtau,rt_tauRuff(:,1:3,n),rt_chiRoss(:,1:3,n), &
          rt_heatFlx(:,1:3,n),rt_heatErms(1:3,n),rt_heatEave(1:3,n), &
          rt_leakNumRad,rt_leakRadii,rt_dr,compos)   
  end do ! n
  !$omp enddo
  !$omp end parallel
  call Timers_stop("calc")
  
  call Timers_start("communication")     
  call MPI_Allgatherv(MPI_IN_PLACE,0, MPI_DATATYPE_NULL, &
                      rt_rayData, rt_recvCnt, rt_dsplCnt, FLASH_REAL, &
                      rt_subMeshComm, error)
  call Timers_stop("communication")
#if NDIM ==1
  rt_leakEaveTot(1:3) = rt_heatEave(1:3,1)
  do m=1,3
     rt_radNu(m) = maxval(rt_leakRadii, rt_tauRuff(:,m,1) > twoThirds) 
  end do
#else
#if NDIM==2
  do m=1,3
     do n=1,rt_leakNumRays
        ns_rad(n) = max(0.,maxval(rt_leakRadii, mask= rt_tauRuff(:,m,n) > twoThirds))
     end do
     rt_leakEaveTot(m) = sum(rt_heatEave(m,:)*sin(rt_leakTheta))/(sum(sin(rt_leakTheta)))
     rt_radNu(m) = sum(ns_rad*sin(rt_leakTheta))/(sum(sin(rt_leakTheta)))
  end do
#else
#if NDIM==3
  do m=1,3
     do n=1,rt_leakNumRays
        ns_rad(n) = max(0.,maxval(rt_leakRadii, mask= rt_tauRuff(:,m,n) > twoThirds))
     end do
     rt_leakEaveTot(m) = sum(rt_heatEave(m,:))/rt_leakNumRays
     rt_radNu(m) = sum(ns_rad)/rt_leakNumRays
  end do
#endif
#endif
#endif

end subroutine rt_calcTau

subroutine calc_taus(rho,temp,ye,oldtauruff,tauruff,chiross, &
     heatflux,heaterms,heateave,nzones,rad,ds,compos)

  use eosmodule, only: eos_tempmin
  use rt_data, ONLY : rt_leakDoHeat
!  use rt_utils, ONLY : get_fermi_integral

  implicit none
!  DECLARE_CCTK_PARAMETERS
  integer, intent(in) :: nzones ! number of radial zones in ray
  real, intent(inout) :: rho(nzones) ! density in g/cm^3
  real, intent(inout) :: temp(nzones) ! temperature in MeV
  real, intent(inout) :: ye(nzones) ! ye, dimensionless

  real, intent(in) :: rad(nzones) !radial points, cm
  real, intent(in) :: ds(nzones) !line element, sqrt(g_{rr}) * dr

  real, intent(inout) :: oldtauruff(nzones,3) !tau used for leakage from last iteration, one for each nu
  real, intent(out) :: tauruff(nzones,3) !new tau used for leakage from current iteration, one for each nu

  real, intent(out) :: compos(nzones,4)
  real, intent(out) :: chiross(nzones,3) !chi calculated from Rosswog scheme, one for each nu, to be interpolated throughout 3D

  real, intent(out) :: heatflux(nzones,3) !flux of neutrinos used in heating
  real, intent(out) :: heaterms(3) !rms energy of neutrinos at neutrinosphere, one for each nu, to be interpolated throughout 3D
  real, intent(out) :: heateave(3) !average energy of neutrinos at neutrinosphere, one for each nu, to be interpolated throughout 3D

  logical :: have_old_tau !whether we have rufftau or need to calculate from stratch

  !EOS & local variables
  integer keytemp, keyerr
  real :: precision = 1.0d-10
  real :: matter_rho,matter_temperature,matter_ye
  real :: matter_enr,matter_prs,matter_ent
  real :: matter_cs2,matter_dedt,matter_dpdrhoe
  real :: matter_dpderho,matter_xa,matter_xh
  real :: matter_xn,matter_xp,matter_abar
  real :: matter_zbar,matter_mue,matter_mun
  real :: matter_mup,matter_muhat

  real :: eta_e(nzones),eta_p(nzones),eta_n(nzones)
  real :: eta_nue(nzones),eta_nua(nzones),eta_nux(nzones),eta_hat(nzones)
  real :: eta_pn(nzones),eta_np(nzones)
  real :: massfracs_xa(nzones),massfracs_xh(nzones)
  real :: massfracs_xp(nzones),massfracs_xn(nzones)
  real :: massfracs_abar(nzones),massfracs_zbar(nzones)

  integer :: i !counter
  integer :: rl

  !heating variables
  real :: radial_luminosity(2),lumrad(nzones,3)
  integer :: ns_location(3)
  real :: lum(3)
  real :: leak_dummy1,leak_dummy2,leak_dummy3,leak_dummy4(3),leak_dummy5(3)
  real :: leak_dummy6(3)
  real :: lepton_blocking(nzones,2)

  !constants & parameters
  real, parameter :: Qnp = 1.293333d0 !m_n - m_p, MeV
  real, parameter :: Cv = 0.5d0 + 2.0e0*0.23d0 !vector coupling
  real, parameter :: Ca = 0.5d0 !axial coupling
  real, parameter :: alpha = 1.23d0 !gA
  real, parameter :: me_mev = 0.510998910e0 !mass of electron in MeV
  real, parameter :: sigma_0 = 1.76d-44 ! in units of cm^2
  real, parameter :: avo = 6.0221367d23 !Avogadro's number
  real, parameter :: pi = 3.1415926535897932384d0
  real, parameter :: twothirds = 2.0e0/3.0e0

  !Ruffert tau stuff
  real :: kappa_tot(nzones,3) !total kappa, 1/cm
  real :: kappa_tot_p(nzones,3) !previous kappa, for comparing, 1/cm
  real :: kappa_scat_n(nzones,3) !total scattering kappa on neutrons, 1/cm
  real :: kappa_scat_p(nzones,3) !total scattering kappa on protons, 1/cm
  real :: kappa_abs_n(nzones) !total absorption kappa on neutrons, 1/cm
  real :: kappa_abs_p(nzones) !total absorption kappa on protons, 1/cm

  real :: local_eta_nue(nzones) !interpolated eta's, nue
  real :: local_eta_nua(nzones) !interpolated eta's, nua
  real :: local_eta_nux(nzones) !interpolated eta's, nux

  real :: csn_0, csp_0 !crosssection preambles
  real :: kappa_const_scat_p(nzones)
  real :: kappa_const_scat_n(nzones)
  real :: kappa_const_abs(nzones)

  real :: xerr
  real :: xerr_out = 1.0d-10
  real :: xlye,xyn,xyp,xynp,xypn,t1,t2
  
  integer :: icount
  integer, parameter :: icount_max = 200

  !function declarations
  real :: get_fermi_integral

  !Rosswog chi variables
  real :: scattering_kappa
  real :: kappa_tilde_nu_scat(nzones,3,3)
  real :: kappa_tilde_nu_abs(nzones,3,3)
  real :: block_factor
  real :: abs_kappa

  real :: zeta(nzones,3) !zeta
  real :: chi(nzones,3) !chi

  ! cactus related stuff
  character(len=512) :: warnline

!#############################

  if(oldtauruff(1,1).gt.0.0e0) then
     have_old_tau = .true.
  else
     have_old_tau = .false.
  endif

  !first get all the EOS variables
  keytemp = 1
  keyerr = 0
  do i=1,nzones

     ! fix potential undershoots near shock
     if(rho(i).lt.1.0d3.and.i.lt.nzones-1.and.i.gt.1) then
        rho(i) = 0.5d0*(rho(i-1)+rho(i+1))
     endif

     ! fix potential undershoots near shock
     ! need to make sure that the temperature in particular does
     ! nothing crazy near the shock
     ! 1) make sure it is not within a factor of 2 of the table
     ! bound
     if(temp(i).lt.2.0e0*eos_tempmin.and.i.lt.nzones-1.and.i.gt.1) then
        temp(i) = 0.5d0*(temp(i-1)+temp(i+1))
     endif
     ! 2) make sure it does not drop by more than 50% from one zone to the
     !    next
     if(i.lt.nzones .and. i.gt.1) then
        if( temp(i) .lt. 0.5d0*temp(i-1) ) then
           temp(i) = 0.5d0*(temp(i-1)+temp(i+1))
        endif
     endif

     ! fix potential undershoots near shock
     if(ye(i).lt.0.0e0.or.ye(i).gt.0.53d0.and.(i.lt.nzones-1.and.i.gt.1)) then
        ye(i) = 0.5d0*(ye(i-1)+ye(i+1))
     endif
     
     matter_rho = rho(i)
     matter_temperature = temp(i)
     matter_ye = ye(i)

     call nuc_eos_full(matter_rho,matter_temperature,matter_ye,matter_enr, &
          matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho, &
          matter_dpdrhoe,matter_xa,matter_xh,matter_xn,matter_xp,matter_abar, &
          matter_zbar,matter_mue,matter_mun,matter_mup,matter_muhat, &
          keytemp,keyerr,precision)
#if 0
     if (keyerr.ne.0) then
        !$OMP CRITICAL(leaktau1)
        write(warnline,"(A15,1P10E15.6)") "rho: ", matter_rho
        call CCTK_WARN(1,warnline)
        write(warnline,"(A15,1P10E15.6)") "temperature: ", matter_temperature
        call CCTK_WARN(1,warnline)
        write(warnline,"(A15,1P10E15.6)") "ye: ", matter_ye
        call CCTK_WARN(1,warnline)
        write(warnline,"(A15,i10)") "eos error", keyerr
        call CCTK_WARN(1,warnline)
        write(warnline,"(A15,i10,1P10E15.6)") "location: ", i, rad(i)
        call CCTK_WARN(1,warnline)
        call CCTK_WARN(0,"set_eos_variables: EOS error in leakage calc_tau")
        !$OMP END CRITICAL(leaktau1)
     endif
#endif
     compos(i,1) = matter_xh
     compos(i,2) = matter_xn
     compos(i,3) = matter_xp
     compos(i,4) = matter_xa

     ! in our EOS the rest mass difference is in the chemical potentials of the neucleons
     eta_e(i) = matter_mue/matter_temperature
     eta_p(i) = matter_mup/matter_temperature
     eta_n(i) = matter_mun/matter_temperature
          
     massfracs_xa(i) = matter_xa
     massfracs_xh(i) = matter_xh
     massfracs_xp(i) = matter_xp
     massfracs_xn(i) = matter_xn
     massfracs_abar(i) = matter_abar
     massfracs_zbar(i) = matter_zbar
          
     eta_hat(i) = eta_n(i)-eta_p(i) - Qnp/matter_temperature
     eta_nue(i) = eta_e(i) - eta_n(i) + eta_p(i) !fully includes effects of rest masses
     eta_nua(i) = -eta_nue(i)
     eta_nux(i) = 0.0e0     
     
  enddo

!#############################
  
  !now calculate Ruffert tau

  !use previous tauruff
  if (have_old_tau) tauruff = oldtauruff

  !initialize
  kappa_tot(:,:)   = 1.0e0 ! 1/cm
  kappa_tot_p(:,:) = 1.0e0 ! 1/cm
  kappa_scat_n(:,:) = 1.0d-5 ! 1/cm
  kappa_scat_p(:,:) = 1.0d-5 ! 1/cm
  kappa_abs_n(:)  = 1.0d-5 ! 1/cm
  kappa_abs_p(:)  = 1.0d-5 ! 1/cm

  local_eta_nux(:) = 0.0e0
  local_eta_nue(:) = 0.0e0
  local_eta_nua(:) = 0.0e0

  !cross section coeffs
  csn_0 = (1.0e0 + 5.0e0*alpha**2) / 24.0e0
  csp_0 = (4.0e0*(Cv-1.0e0)**2 + 5.0e0*alpha**2) / 24.0e0

  do i=1,nzones
     ! constant parts of kappa (A6)
     t1 = sigma_0 * avo * rho(i) * (temp(i)/me_mev)**2
     kappa_const_scat_n(i) = csn_0 * t1
     kappa_const_scat_p(i) = csp_0 * t1
     ! (A11) constant part
     kappa_const_abs(i) = (1.0e0+3.0e0*alpha**2)/4.0e0 * t1
  enddo
  
  ! Loop to get converged result for tau.
  ! This is discussed in the text between equations (A5) and
  ! (A6). Note that for the initial iteration tau is set to
  ! 1 and the etas are set to 1.0d-5
  
  icount = 1
  xerr = 1.0e0

  do while(xerr.gt.xerr_out .and. icount.lt.icount_max)
     ! copy over current into previous kappa
     kappa_tot_p = kappa_tot

     ! set up new kappa based on individual contributions
     ! nu_e; (A17)
     kappa_tot(:,1) = &
          kappa_scat_p(:,1) &
          + kappa_scat_n(:,1) &
          + kappa_abs_n(:)
     ! antis; (A18)
     kappa_tot(:,2) = &
              kappa_scat_p(:,2) &
            + kappa_scat_n(:,2) &
            + kappa_abs_p(:)

     ! nu_xs (A19)
     kappa_tot(:,3) = &
          + kappa_scat_p(:,3) &
          + kappa_scat_n(:,3) !&
       
     ! Integrate optical depths: Equation (A20)
     ! Note that this is not done for particle and energy transport
     if(icount.gt.2.or..not.have_old_tau) then
        tauruff(:,:) = 0.0e0
        do i=nzones-1,1,-1
           tauruff(i,1:3) = tauruff(i+1,1:3) + kappa_tot(i,1:3)*ds(i)
           compos(i,4) = kappa_tot(i,1)
        enddo
     endif

     do i=1,nzones
        ! (A5) equilibrium eta, we have rest masses in our chemical potentials
        ! no need to include mass difference in eta
        local_eta_nux(i) = 0.0e0   ! (A2)
        local_eta_nue(i) = eta_nue(i) * (1.0e0-exp(-tauruff(i,1))) ! (A3); note that the ^0 etas are set to 0.0e0
        local_eta_nua(i) = eta_nua(i) * (1.0e0-exp(-tauruff(i,2))) ! (A4)
        
        
        !assuming completely dissociated
        xlye = ye(i)
        xyn = (1.0e0-xlye) / (1.0e0 + 2.0e0/3.0e0 * max(eta_n(i),0.0e0)) ! (A8)
        xyp = xlye / (1.0e0 + 2.0e0/3.0e0*max(eta_p(i),0.0e0))

        if(massfracs_xh(i).lt.0.5d0) then
           t1 = exp(-eta_hat(i))
           xynp = max((2.0e0*xlye-1.0e0)/ (t1-1.0e0),0.0e0) ! (A13)
           xypn = max(xynp * t1,0.0e0) ! (A14)
        else
           xypn = massfracs_xp(i)
           xynp = massfracs_xn(i)
        endif


        ! electron neutrinos
        t1 = get_fermi_integral(5,local_eta_nue(i)) / & 
             get_fermi_integral(3,local_eta_nue(i))
        t2 = 1.0e0 + exp(eta_e(i)-get_fermi_integral(5,local_eta_nue(i)) / &
             get_fermi_integral(4,local_eta_nue(i))) ! (A15)
        
        kappa_scat_n(i,1) = kappa_const_scat_n(i) * xyn  * t1
        kappa_scat_p(i,1) = kappa_const_scat_p(i) * xyp  * t1
        kappa_abs_n(i) = kappa_const_abs(i) * xynp * t1 / t2 ! (A11)
        
        ! anti-electron neutrinos
        t1 = get_fermi_integral(5,local_eta_nua(i)) / & 
             get_fermi_integral(3,local_eta_nua(i))
        t2 = 1.0e0 + exp(-eta_e(i)-get_fermi_integral(5,local_eta_nua(i)) / &
             get_fermi_integral(4,local_eta_nua(i))) ! (A16)
        
        kappa_scat_n(i,2) = kappa_const_scat_n(i) * xyn  * t1 ! (A6)
        kappa_scat_p(i,2) = kappa_const_scat_p(i) * xyp  * t1
        kappa_abs_p(i) = kappa_const_abs(i) * xypn * t1 / t2 ! (A12)
        
        ! nux neutrinos
        t1 = get_fermi_integral(5,local_eta_nux(i)) / & 
             get_fermi_integral(3,local_eta_nux(i))
        kappa_scat_n(i,3) = kappa_const_scat_n(i) * xyn * t1 ! (A6)
        kappa_scat_p(i,3) = kappa_const_scat_p(i) * xyp * t1
        
     enddo
     
     ! compute relative change xerr
     xerr = 0.0e0
     do i=1,nzones
        xerr = max(xerr,abs(kappa_tot(i,1)/kappa_tot_p(i,1)-1.0e0))
        xerr = max(xerr,abs(kappa_tot(i,2)/kappa_tot_p(i,2)-1.0e0))
        xerr = max(xerr,abs(kappa_tot(i,3)/kappa_tot_p(i,3)-1.0e0))
     enddo
     
     icount = icount + 1
     
  enddo

  if(icount.ge.icount_max) then
     write(6,"(i5,1P10E15.6)") icount,xerr,xerr_out
     stop "icount > icount_max in leakage; leak.F90"
  endif

  ! Recompute tau based on the most recent kappa_tot
  tauruff = 0.0e0

  do i=nzones-1,1,-1
     tauruff(i,1:3) = tauruff(i+1,1:3) + kappa_tot(i,1:3)*ds(i)
  enddo

  have_old_tau = .true.

  !set degeneracy factors to interpolated values
  eta_nue(:) = local_eta_nue(:)
  eta_nua(:) = local_eta_nua(:)
  eta_nux(:) = local_eta_nux(:)


!#######################################
  
!Rosswog chi
  
  zeta(:,:) = 0.0e0
  chi(:,:) = 0.0e0
      
  do i=1,nzones

     !scattering
     scattering_kappa = rho(i)*avo*0.25d0*sigma_0/me_mev**2
     kappa_tilde_nu_scat(i,1,1) = massfracs_xn(i)*scattering_kappa
     kappa_tilde_nu_scat(i,1,2) = massfracs_xp(i)*scattering_kappa
     kappa_tilde_nu_scat(i,2,1) = massfracs_xn(i)*scattering_kappa
     kappa_tilde_nu_scat(i,2,2) = massfracs_xp(i)*scattering_kappa
     kappa_tilde_nu_scat(i,3,1) = massfracs_xn(i)*scattering_kappa
     kappa_tilde_nu_scat(i,3,2) = massfracs_xp(i)*scattering_kappa
       

     scattering_kappa = rho(i)*avo*0.0625d0*sigma_0/me_mev**2* &
          massfracs_abar(i)*(1.0e0-massfracs_zbar(i)/massfracs_abar(i))**2 ! only have 1 factor of A because kappa multiples the number fraction, not mass fractions
     kappa_tilde_nu_scat(i,1,3) = massfracs_xh(i)*scattering_kappa
     kappa_tilde_nu_scat(i,2,3) = massfracs_xh(i)*scattering_kappa
     kappa_tilde_nu_scat(i,3,3) = massfracs_xh(i)*scattering_kappa

     eta_pn(i) = avo*rho(i)*(massfracs_xn(i)-massfracs_xp(i))/(exp(eta_hat(i))-1.0e0)
     eta_pn(i) = max(eta_pn(i),0.0e0)
     eta_np(i) = avo*rho(i)*(massfracs_xp(i)-massfracs_xn(i))/(exp(-eta_hat(i))-1.0e0)
     eta_np(i) = max(eta_np(i),0.0e0)

     if (rho(i).lt.1.0d11) then
        !non degenerate here, use mass fractions as chemical potentials fail at low densities
        eta_pn(i) = avo*rho(i)*massfracs_xp(i)
        eta_np(i) = avo*rho(i)*massfracs_xn(i)
     endif

     !absorption
     abs_kappa = (1.0e0+3.0e0*alpha**2)*0.25d0*sigma_0/me_mev**2
     block_factor = 1.0e0 + exp(eta_e(i)-get_fermi_integral(5,eta_nue(i))/ &
          get_fermi_integral(4,eta_nue(i)))
     kappa_tilde_nu_abs(i,1,1) = eta_np(i)*abs_kappa/block_factor
     kappa_tilde_nu_abs(i,2,1) = 0.0e0 !no absorption of a-type on neutrons
     kappa_tilde_nu_abs(i,3,1) = 0.0e0 !no absorption of x-type neutrinos
     kappa_tilde_nu_abs(i,1,2) = 0.0e0 !no absorption of e-type on protons
     block_factor = 1.0e0 + exp(-eta_e(i)-get_fermi_integral(5,eta_nua(i))/ &
          get_fermi_integral(4,eta_nua(i)))
     kappa_tilde_nu_abs(i,2,2) = eta_pn(i)*abs_kappa/block_factor
     kappa_tilde_nu_abs(i,3,2) = 0.0e0 !no absorption of x-type neutrinos
     kappa_tilde_nu_abs(i,1,3) = 0.0e0 !no absorption on nuclei
     kappa_tilde_nu_abs(i,2,3) = 0.0e0 !no absorption on nuclei
     kappa_tilde_nu_abs(i,3,3) = 0.0e0 !no absorption on nuclei

     !sum up opacities to get zeta (again, factoring out energy dependence)
     zeta(i,1) = kappa_tilde_nu_scat(i,1,1) + kappa_tilde_nu_scat(i,1,2) + &
          kappa_tilde_nu_scat(i,1,3) + kappa_tilde_nu_abs(i,1,1) + &
          kappa_tilde_nu_abs(i,1,2) + kappa_tilde_nu_abs(i,1,3)

     zeta(i,2) = kappa_tilde_nu_scat(i,2,1) + kappa_tilde_nu_scat(i,2,2) + &
          kappa_tilde_nu_scat(i,2,3) + kappa_tilde_nu_abs(i,2,1) + &
          kappa_tilde_nu_abs(i,2,2) + kappa_tilde_nu_abs(i,2,3)

     zeta(i,3) = kappa_tilde_nu_scat(i,3,1) + kappa_tilde_nu_scat(i,3,2) + &
          kappa_tilde_nu_scat(i,3,3) + kappa_tilde_nu_abs(i,3,1) + &
          kappa_tilde_nu_abs(i,3,2) + kappa_tilde_nu_abs(i,3,3)
     
  enddo
  
  do i=nzones-1,1,-1
     !integrate zeta to get chi, tau with energy dependence factored out
     chi(i,1) = chi(i+1,1) + zeta(i,1)*ds(i)
     chi(i,2) = chi(i+1,2) + zeta(i,2)*ds(i)
     chi(i,3) = chi(i+1,3) + zeta(i,3)*ds(i)
  enddo

  chi(nzones,:) = chi(nzones-1,:)
  
  chiross = chi

  !neutrinosphere located at tau = 2/3
  ns_location(:) = 1
  do i=1,nzones
     if (tauruff(i,1).gt.twothirds) then
        ns_location(1) = i
     endif
     if (tauruff(i,2).gt.twothirds) then
        ns_location(2) = i
     endif
     if (tauruff(i,3).gt.twothirds) then
        ns_location(3) = i
     endif
  enddo

  !rms energy at neutrino sphere
  heaterms(1) = temp(ns_location(1))* &
       sqrt(get_fermi_integral(5,eta_nue(ns_location(1)))/ &
       get_fermi_integral(3,eta_nue(ns_location(1))))
  heaterms(2) = temp(ns_location(2))* &
       sqrt(get_fermi_integral(5,eta_nua(ns_location(2)))/ &
       get_fermi_integral(3,eta_nua(ns_location(2))))
  heaterms(3) = temp(ns_location(3))* &
       sqrt(get_fermi_integral(5,eta_nux(ns_location(3)))/ &
       get_fermi_integral(3,eta_nux(ns_location(3))))

  !mean energy at neutrino sphere
  heateave(1) = temp(ns_location(1))* &
       get_fermi_integral(5,eta_nue(ns_location(1)))/ &
       get_fermi_integral(4,eta_nue(ns_location(1)))
  heateave(2) = temp(ns_location(2))* &
       get_fermi_integral(5,eta_nua(ns_location(2)))/ &
       get_fermi_integral(4,eta_nua(ns_location(2)))
  heateave(3) = temp(ns_location(3))* &
       get_fermi_integral(5,eta_nux(ns_location(3)))/ &
       get_fermi_integral(4,eta_nux(ns_location(3)))

  !now we leak along this lines to determine luminosity, will use this for heating!

  radial_luminosity(:) = 0.0e0
  lumrad(:,:) = 0.0e0
  heatflux = 0.0e0

  if(rt_leakDoHeat) then
     do i =1,nzones
        lepton_blocking(i,1) = 1.0e0/(1.0e0 + exp(eta_e(i) - &
             get_fermi_integral(5,eta_nue(ns_location(1)))/ &
             get_fermi_integral(4,eta_nue(ns_location(1)))))
        lepton_blocking(i,2) = 1.0e0/(1.0e0 + exp(-eta_e(i) - & 
             get_fermi_integral(5,eta_nua(ns_location(2)))/ &
             get_fermi_integral(4,eta_nua(ns_location(2)))))
     enddo
     
     rl = -1 
     do i=1,nzones-1
     
        leak_dummy1 = 0.0e0
        leak_dummy2 = 0.0e0
        leak_dummy3 = 0.0e0 !this is ldt, the actual leakage will take care of overflow of dyedt
        leak_dummy4 = 0.0e0 
        leak_dummy5 = 0.0e0 
        leak_dummy6 = 0.0e0 

        !call leak with heating turned on to get the change in
        !luminosity, then we'll interpolate that and use for heating
        call calc_leak(rho(i),temp(i),ye(i),chi(i,:),tauruff(i,:),heatflux(i,:), &
             heaterms,heateave,leak_dummy1,leak_dummy2,leak_dummy3,&
             lum,leak_dummy4,leak_dummy5,leak_dummy6,rl)

        !lum coming from calc_leak is not actually luminosity, rather, 
        !must multiply by volume (of spherical shell)
        !then add to integrated luminosity
        radial_luminosity(1) = radial_luminosity(1) + lum(1)*4.0e0*pi*((rad(i)+rad(i+1))/2.0e0)**2*ds(i)
        radial_luminosity(2) = radial_luminosity(2) + lum(2)*4.0e0*pi*((rad(i)+rad(i+1))/2.0e0)**2*ds(i)
        lumrad(i,1) = leak_dummy5(1)
        lumrad(i,2) = leak_dummy5(2)
     
        !this is what gets interpolated and put into the leak routine.
        heatflux(i+1,1) = radial_luminosity(1)/(4.0e0*pi*((rad(i)+rad(i+1))/2.0e0)**2)&
             * lepton_blocking(i,1) !luminosity/4pir^2
        heatflux(i+1,2) = radial_luminosity(2)/(4.0e0*pi*((rad(i)+rad(i+1))/2.0e0)**2)&
             * lepton_blocking(i,2) !luminosity/4pir^2

     enddo

#if 0
     !$OMP CRITICAL(leak333)
     open(666,file="lum.dat")
     do i=1,nzones-1
        write(666,"(1P10E15.6)") rad(i),heatflux(i,1),heatflux(i,2),lumrad(i,1:2)*ds(i)
     enddo
     close(666)
     call CCTK_WARN(0,"end debug")
     ! fix call to this routine
     !$OMP END CRITICAL(leak333)
#endif

  else
     heaterms = 0.0e0
     heateave = 0.0e0
  endif


end subroutine calc_taus
