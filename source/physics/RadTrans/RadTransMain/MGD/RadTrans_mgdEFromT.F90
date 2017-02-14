!!****if* source/physics/RadTrans/RadTransMain/MGD/RadTrans_mgdEFromT
!!
!!  NAME 
!!
!!  RadTrans_mgdEFromT
!!
!!  SYNOPSIS
!!
!!  call RadTrans_mgdEFromT( integer(IN) :: blockId,
!!                           integer(IN) :: axis(MDIM),
!!                           real(IN)    :: trad,
!!                           real(OUT), optional :: tradActual )
!!
!!  DESCRIPTION 
!!
!!  This routine uses a radiation temperature to set the group
!!  specific radiation energies for all of the groups owned by this
!!  process in a given cell and block. It also sets the total specific
!!  radiation energy.
!!
!!  This routine was created to make it easy to initialize the group
!!  energies from a temperature in Simulation_initBlock.
!!
!!  The final, optional argument is tradActual. This is the radiation
!!  temperature computed from the specific energy. This is necessary
!!  to handle the case where the highest radiation energy group
!!  boundary is too small. In this case part of the black body
!!  spectrum will be cutoff and ERAD_VAR will not exactly equal
!!  a*trad**4. Thus, this routine will optionally return the correct
!!  value of trad which accounts for this possibility. Users can use
!!  this value to set TRAD_VAR in Simulation_initBlock.
!!
!!  ARGUMENTS
!!
!!    blockId    : The blockId of the cell
!!    axis       : An array storing the i,j,k coordinate of the cell
!!    trad       : The radiation temperature (K)
!!    tradActual : The actual radiation temperature computed from erad
!!***
subroutine RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
  use RadTrans_data, ONLY: rt_boltz, rt_speedlt, rt_radconst, rt_acrossMe, &
       rt_meshCopyCount, rt_useRadTrans
  use rt_data, ONLY: rt_mgdNumGroups
  
  use RadTrans_interface, ONLY: RadTrans_mgdGetBound, RadTrans_planckInt
  use Grid_interface, ONLY: Grid_putPointData, Grid_getPointData
  implicit none

#include "constants.h"
#include "Flash.h"

  ! Arguments:
  integer, intent(in) :: blockId
  integer, intent(in) :: axis(MDIM)
  real,    intent(in) :: trad
  real,    intent(out), optional :: tradActual

  ! Local Variables:
  integer :: n        ! Loop counter for number of groups
  integer :: g        ! The group number
  integer :: ig       ! The group unk index
  real    :: rho      ! cell mass density (g/cc)
  real    :: erad     ! the specific radiation energy for this cell (ergs/g)
  real    :: eg       ! energy for group boundary g
  real    :: egp1     ! energy for group boundary g+1
  real    :: xg       ! h nu/k Te for group boundary g
  real    :: xgp1     ! h nu/k Te for group boundary g+1
  real    :: pxg      ! Planck integral for group boundary g
  real    :: pxgp1    ! Planck integral for group boundary g+1
  real    :: erad_tot ! Total specific radiation energy for this cell (ergs/g)

  ! Physical constants:
  real :: C  ! Speed of light [cm/s]
  real :: A  ! Radiation constant [ergs/cm^3/K^4]
  real :: KB ! Boltzmann constant [ergs/K]  

  if(.not. rt_useRadTrans) then 
     if(present(tradActual)) tradActual = trad
     return 
  end if

  KB = rt_boltz
  C  = rt_speedlt
  A  = rt_radconst

  erad_tot = 0.0
  do n = 1, NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)
     ! Get the group number for this group:
     g = NONREP_LOC2GLOB(n, rt_acrossMe, rt_meshCopyCount)
     ig = MGDR_NONREP_LOC2UNK(n)

     call RadTrans_mgdGetBound(g, eg)
     call RadTrans_mgdGetBound(g+1, egp1)
     xg   = eg / (KB*trad)
     xgp1 = egp1 / (KB*trad)
     call RadTrans_planckInt(xg, pxg)
     call RadTrans_planckInt(xgp1, pxgp1)
     pxg = pxg * 15.0/PI**4
     pxgp1 = pxgp1 * 15.0/PI**4

     call Grid_getPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
     erad = A * trad**4 * (pxgp1 - pxg) / rho

     call Grid_putPointData(blockId, CENTER, ig, EXTERIOR, axis, erad)
     erad_tot = erad_tot + erad
  end do

  call Grid_putPointData(blockId, CENTER, ERAD_VAR, EXTERIOR, axis, erad_tot)

  if(present(tradActual)) then
     ! Compute the actual radiation temperature. This accounts for the
     ! fact that some of the energy may have been cutoff by the fact
     ! that the upper radiation group boundary was too low.
     call RadTrans_mgdGetBound(1, eg)
     call RadTrans_mgdGetBound(rt_mgdNumGroups+1, egp1)
     xg   = eg / (KB*trad)
     xgp1 = egp1 / (KB*trad)
     call RadTrans_planckInt(xg, pxg)
     call RadTrans_planckInt(xgp1, pxgp1)
     tradActual = (trad**4 * (pxgp1-pxg)*15.0/PI**4)**0.25
  end if
  return

end subroutine RadTrans_mgdEFromT
