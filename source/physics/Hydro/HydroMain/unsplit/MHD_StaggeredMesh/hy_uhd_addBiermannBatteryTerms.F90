!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_addBiermannBatteryTerms
!!
!! NAME
!!
!!  hy_uhd_addBiermannBatteryTerms
!!
!!
!! SYNOPSIS
!!
!!  call hy_uhd_addBiermannBatteryTerms(integer(IN) :: blockID,
!!                            integer(IN) :: blkLimitsGC(LOW:HIGH,MDIM),
!!                            integer(IN) :: ix,
!!                            integer(IN) :: iy,
!!                            integer(IN) :: iz,
!!                            real(IN)    :: Flux,
!!                            integer(IN) :: sweepDir)
!!
!!
!! DESCRIPTION
!!
!!  Adds Biermann battery terms to total MHD fluxes. This is a stub. We are not supporting this yet. 
!!
!!  d  / Bx \    d /  0  \    d /  Ez \     d / -Ey \
!! ----| By | + ---| -Ez | + ---|  0  |  + ---|  Ex | = 0, 
!!  dt \ Bz /   dx \  Ey /   dy \ -Ex /    dz \  0  /
!!
!! dEner/dt + grad( E x B ) (mind the sign!)
!!
!! where Ohm's law for Biermann battery terms give electric fields by
!! E = -gradPe/(qele * n_e) = -1/(qele * n_e)*[i*dPe/dx + j*dPe/dy + k*dPe/dz] (method 1)
!! E = kboltz * ln(Pe)  gradTe / (qele) = kboltz * LOG(Pe)/(qele) * [i*dTe/dx + j*dTe/dy + k*dTe/dz] (method 2)
!!
!! For the energy update we either compute the flux as ExB from method 1 or use the flux 
!! F = k/q (- lnPe B x gradTe + Te lnPe curlB)
!! i,j,k here are unit vectors in x,y,z directions.
!!
!! In cylindrical coords (2D) we can only have dtBphi + dzEr - drEz
!!
!!
!! ARGUMENTS
!!
!!  blockID     - a local blockID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section 
!!                of block with the guard cells 
!!  ix,iy,iz    - indices of the line along which the sweep is made
!!  Flux        - array containing MHD fluxes
!!  sweepDir    - direction of sweep
!!
!!***


Subroutine hy_uhd_addBiermannBatteryTerms(blockID,blkLimitsGC,ix,iy,iz,Flux,sweepDir)


  implicit none

#include "constants.h"
#include "UHD.h"

  !! Argument List ----------------------------------------------------------
  integer, INTENT(IN) :: blockID,ix,iy,iz
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
  real, dimension(HY_VARINUM), intent(INOUT) :: Flux
  integer, INTENT(IN) :: sweepDir
  !! ----------------------------------------------------------------------

  return
  
End Subroutine hy_uhd_addBiermannBatteryTerms


Subroutine get_upwind(vel,pe_L,pe_R,pe_Up)

  use hy_uhd_slopeLimiters, ONLY : signum

  implicit none
  real, intent(IN)  :: vel,pe_L,pe_R
  real, intent(OUT) :: pe_Up

  real :: velP,velN

  velP = 0.5*(1.+signum(vel)) !*abs(velo)
  velN = 0.5*(1.-signum(vel)) !*abs(velo)

  pe_Up = velP*pe_L + velN+pe_R

end Subroutine get_upwind
