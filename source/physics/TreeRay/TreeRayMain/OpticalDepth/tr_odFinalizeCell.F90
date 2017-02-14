!!****if* source/physics/TreeRay/TreeRayMain/OpticalDepth/tr_odFinalizeCell
!!
!! NAME
!!
!!  tr_odFinalizeCell
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***
#include "Flash.h"

subroutine tr_odFinalizeCell(solnPoint, dl_poc, cdMaps)
  use TreeRay_data, ONLY : tr_nCd, tr_nPix
  use tr_odData, ONLY : tr_odICDTO, tr_odICDH2, tr_odICDCO
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
  use Chemistry_data, ONLY : ch_muC, ch_mH, ch_abar, ch_mf_scale
#endif

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  real, DIMENSION(:), POINTER, intent(INOUT) :: solnPoint
  real, DIMENSION(MDIM), intent(IN) :: dl_poc
  real, DIMENSION(tr_nCd, 0:tr_nPix-1), intent(IN) :: cdMaps

  integer :: ipix
  real :: minCD, yn, temp, abh2, abco, rho, dl
  real :: AV_mean, chi_mean, fshield_H2, fshield_CO
  real, DIMENSION(1:tr_nPix) :: ynCdMap, nH2CdMap, nCOCdMap

  ! No chemistry, just calculate minimum column density of gas
#ifndef CHEMISTRYNETWORK
  minCD = 1.d99
  do ipix = 1, tr_nPix-1
    if (cdMaps(tr_odICDTO,ipix) < minCD) then
      minCD = cdMaps(tr_odICDTO,ipix)
    endif
  enddo
  !print *, minCD, cdMaps
  solnPoint(CDTO_VAR) = minCD
#endif

  ! Chemistry networks of SILCC
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
  rho = solnPoint(DENS_VAR)
  yn  = rho/ch_mH/ch_abar
  temp = solnPoint(TEMP_VAR)
  abh2 = 0.5 * solnPoint(IH2_SPEC)

  ynCdMap(1:tr_nPix)  = cdMaps(tr_odICDTO, 0:tr_nPix-1)/ch_mH/ch_abar
  nH2CdMap(1:tr_nPix) = cdMaps(tr_odICDH2, 0:tr_nPix-1) /ch_mH/2.e0/ch_abar

#if CHEMISTRYNETWORK == 5
  abh2 = solnPoint(IH2_SPEC) * ch_mf_scale/ 2.e0 
  abco = solnPoint(ICO_SPEC) * ch_mf_scale/ch_muC

  nH2CdMap(1:tr_nPix)= cdMaps(tr_odICDH2, 0:tr_nPix-1) /ch_mH/2.e0*ch_mf_scale/ch_abar
  nCOCdMap(1:tr_nPix)= cdMaps(tr_odICDCO, 0:tr_nPix-1) /ch_mH/ch_muC*ch_mf_scale/ch_abar
#else
  abco = 0.e0
  nH2CdMap(1:tr_nPix) = 0.0
#endif

  dl = 0.5e0*(dl_poc(IAXIS)+dl_poc(JAXIS)+dl_poc(KAXIS))/3.e0
  call calc_shielding(yn, dl, temp, abh2, abco, tr_nPix, ynCdMap, &
  & nH2CdMap, nCOCdMap, fshield_H2, fshield_CO, AV_mean, chi_mean)

  solnPoint(CDTO_VAR) = AV_mean
  solnPoint(CHID_VAR) = chi_mean
  solnPoint(CDH2_VAR) = fshield_H2
#if CHEMISTRYNETWORK == 5
  solnPoint(CDCO_VAR) = fshield_CO
#endif

#endif
  

  return
end subroutine tr_odFinalizeCell

