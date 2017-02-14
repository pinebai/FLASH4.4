!!****if* source/physics/sourceTerms/Flame/FlameSpeed/BuoyancyCompensation/CONe/fl_fsAtwoodInitTable
!!
!! NAME
!!
!!  fl_fsAtwoodInitTable
!!
!! SYNOPSIS
!!
!!  call fl_fsAtwoodInitTable()
!!
!! DESCRIPTION
!!
!! Aaron Jackson, Dean M. Townsley, Alan Calder 2008
!!
!! Tabulate an estimate of the atwood number for the flame
!! as a function of density.  Only the exansion for the first
!! stage of burning, C12->Mg24, is considered.
!!
!! For simplicity this version makes several assumptions about the
!! fuel:
!! pure CO, with constant carbon fraction taken from the "c_frac" parameter declared
!!          in the parametric burn unit
!! fuel temperature = 10^8 K -- this be degenerate for the density range addressed
!! 
!! the density range calculated is drawn from the quenching density
!! used for turning off the buoyancy compensating flamespeed enhancement
!! and 4x10^9 g/cc, which is about as dens as we would want a WD to be
!! for a standard Type Ia simulation.  Also at higher densities the atwood number
!! is quite small.
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***


subroutine fl_fsAtwoodInitTable()

  use fl_fsData, ONLY : fl_fsQuenchDens0
  use fl_fsAtwoodData, ONLY : fl_fsAtwoodTabMinLdens, fl_fsAtwoodTabMaxLdens, &
                              fl_fsAtwoodTabDLdens, fl_fsAtwoodTabLdensNum, &
                              fl_fsAtwoodTabMinXc12, fl_fsAtwoodTabMaxXc12, &
                              fl_fsAtwoodTabDXc12, fl_fsAtwoodTabXc12Num, &
                              fl_fsAtwoodTabA
  use Flame_interface, ONLY : Flame_rhJump
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use bn_paraInterface, ONLY : bn_paraFuelAshProperties

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer :: istat, i, j
  real, parameter :: MeVperAMUtoErgperGram = 9.6485e+17
  real :: X_C12, q, q_u, q_b, s
  real, save :: X_Ne22
  real :: temp_u, dens_u, ye_u, yi_u, sumy_u, pres_u, eint_u
!! now..  eosData_u(EOS_TEMP), eosData_u(EOS_DENS),&
!!        eosData_u(EOS_ZBAR)/eosData_u(EOS_ABAR), ....,&
!!        eosData_u(EOS_PRES), eosData_u(EOS_EINT) 
  real :: temp_b, dens_b, ye_b, yi_b, sumy_b, pres_b, eint_b
!! now..  eosData_b(EOS_TEMP), eosData_b(EOS_DENS),&
!!        eosData_b(EOS_ZBAR)/eosData_u(EOS_ABAR), ....,&
!!        eosData_b(EOS_PRES), eosData_b(EOS_EINT)
  real, dimension(EOS_NUM) :: eosData_u, eosData_b
  real, dimension(NSPECIES) :: mfrac_u, mfrac_b

  s = 0.0
! get representative ne_frac
  call RuntimeParameters_get("rep_ne_frac", X_Ne22)
  call RuntimeParameters_get("max_dens", fl_fsAtwoodTabMaxLdens)
  fl_fsAtwoodTabMaxLdens = log(fl_fsAtwoodTabMaxLdens)
  fl_fsAtwoodTabMinLdens = log(fl_fsQuenchDens0)
  call RuntimeParameters_get("max_c_frac", fl_fsAtwoodTabMaxXc12)
  call RuntimeParameters_get("min_c_frac", fl_fsAtwoodTabMinXc12)
  call RuntimeParameters_get("num_ldens",fl_fsAtwoodTabLdensNum)
  call RuntimeParameters_get("num_c_frac",fl_fsAtwoodTabXc12Num)

  ! determine table gridding from fl_fsAtwoodTabNum
  fl_fsAtwoodTabDLdens = (fl_fsAtwoodTabMaxLdens-fl_fsAtwoodTabMinLdens)/fl_fsAtwoodTabLdensNum
  fl_fsAtwoodTabDXc12 = (fl_fsAtwoodTabMaxXc12-fl_fsAtwoodTabMinXc12)/fl_fsAtwoodTabXc12Num


  allocate(fl_fsAtwoodTabA(fl_fsAtwoodTabLdensNum+1,fl_fsAtwoodTabXc12Num+1),stat=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate AtwoodTabA in fl_fsAtwoodInitTable")


  do j = 1, fl_fsAtwoodTabXc12Num+1
     X_C12 = fl_fsAtwoodTabMinXc12+(j-1)*fl_fsAtwoodTabDXc12
     ! only first stage energy release (C12->Mg24)
!     call Flame_heatRelease(q,1)
     call bn_paraFuelAshProperties(X_C12,X_Ne22,ye_u,ye_b,yi_u,yi_b,q_u,q_b)
     q = (q_b - q_u)*MeVperAMUtoErgperGram
!     temp_u = 1.0e8
     eosData_u(EOS_TEMP) = 1.0e8
!     ye_u = 0.5
     ! assume CO mixture
     eosData_u(EOS_ABAR) = 1.0 / yi_u
     eosData_u(EOS_ZBAR) = ye_u * eosData_u(EOS_ABAR)
!     ye_b = 0.5
     ! carbon goes to mg24
     eosData_b(EOS_ABAR) = 1.0 / yi_b
     eosData_b(EOS_ZBAR) = ye_b * eosData_b(EOS_ABAR)

     ! fill in table of densities
     do i = 1, fl_fsAtwoodTabLdensNum+1
        eosData_u(EOS_DENS) = exp(fl_fsAtwoodTabMinLdens+(i-1)*fl_fsAtwoodTabDLdens)
        !print *, dens_u
        call Flame_rhJump(eosData_u, eosData_b, q, s, MODE_DENS_TEMP)
        fl_fsAtwoodTabA(i,j) = (eosData_u(EOS_DENS) - eosData_b(EOS_DENS)) / &
                               (eosData_u(EOS_DENS) + eosData_b(EOS_DENS))
        !print *, fl_fsAtwoodTabA(i,j)
     enddo
  enddo

end subroutine
