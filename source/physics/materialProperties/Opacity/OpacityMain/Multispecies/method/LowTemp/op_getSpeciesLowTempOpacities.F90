!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_getSpeciesLowTempOpacities
!!
!! NAME
!!
!!  op_getSpeciesLowTempOpacities
!!
!! SYNOPSIS
!!
!!  call op_getSpeciesLowTempOpacities (integer (in)  :: species,
!!                                      real    (in)  :: speciesTemperature,
!!                                      integer (in)  :: speciesEnergyGroup)
!!
!! DESCRIPTION
!!
!!  Computes absorption, emission and transport opacities from tabulated low temperature
!!  opacity values for a particular (temperature, energy group, species) triple.
!!  Extracts via interpolation low temperature Planck and Rosseland opacities for a
!!  specific temperature, energy group and species and copies them to the absorption,
!!  emission and transport opacities. The interpolation method is linear.
!!
!!  The mapping between the low temperature Planck and Rosseland opacities and the
!!  absorption, emission and transport opacities is as follows:
!!
!!              absorption opacity = low temperature Planck
!!              emission   opacity = low temperature Planck
!!              transport  opacity = low temperature Rosseland
!!
!! ARGUMENTS
!!
!!  species            : The species index
!!  speciesTemperature : The species temperature (in K)
!!  speciesEnergyGroup : The species energy group index
!!
!!***
subroutine op_getSpeciesLowTempOpacities (species,            &
                                          speciesTemperature, &
                                          speciesEnergyGroup  )

  use Driver_interface,  ONLY : Driver_abortFlash

  use Opacity_data,      ONLY : op_totalSpecies,    &
                                op_nEnergyGroups,   &
                                op_speciesOpacities

  use op_lowTempData,    ONLY : op_minLowTempTemperature,  &
                                op_maxLowTempTemperature,  &
                                op_maxNstepsLowTemp,       &
                                op_tableLowTemp,           &
                                op_PlanckLowTempTables,    &
                                op_RosselandLowTempTables
  implicit none

#include "Opacity.h"

  integer, intent (in)  :: species
  integer, intent (in)  :: speciesEnergyGroup
  real,    intent (in)  :: speciesTemperature

  integer :: i,j,t

  real :: T1,T2
  real :: o1P,o2P,o1R,o2R
  real :: opacityPlanck, opacityRosseland
!
!
!   ...Check if species index  and energy group index are in range.
!
!
  if (species < 1 .or. species > op_totalSpecies) then
      call Driver_abortFlash ('[op_getSpeciesLowTempOpacity] ERROR: bad species index')
  end if

  if (speciesEnergyGroup < 1 .or. speciesEnergyGroup > op_nEnergyGroups) then
      call Driver_abortFlash ('[op_getSpeciesLowTempOpacity] ERROR: bad species energy group index')
  end if
!
!
!   ...Determine the following:
!
!        1) the temperature pair (T1,T2) boundary containing the species's temperature
!        2) the low temp table index pair (i,j) of the temperature boundary
!        3) the two tabulated low temp opacities (o1,o2) corresponding to the T1 and T2
!           temperatures
!
!      If the supplied species temperature lays outside the tabulated boundaries, the
!      routine stops the calculation. The criteria as to when the species's temperature
!      belong within the boundary are:
!
!                         T1 =< speciesTemperature =< T2
!
!
  if (     (speciesTemperature < op_minLowTempTemperature) &
      .or. (speciesTemperature > op_maxLowTempTemperature) ) then
       call Driver_abortFlash ('[op_getSpeciesLowTempOpacity] ERROR: bad species temperature')
  end if

  do t = 2,op_maxNstepsLowTemp
     if (op_tableLowTemp (t) >= speciesTemperature) then
         i  = t - 1
         j  = t
         T1 = op_tableLowTemp (i)
         T2 = op_tableLowTemp (j)
         exit
     end if
  end do

  o1P = op_PlanckLowTempTables    (i,species,speciesEnergyGroup)
  o2P = op_PlanckLowTempTables    (j,species,speciesEnergyGroup)
  o1R = op_RosselandLowTempTables (i,species,speciesEnergyGroup)
  o2R = op_RosselandLowTempTables (j,species,speciesEnergyGroup)
!
!
!   ...Do the linear interpolation.
!
!
  opacityPlanck    = ((o1P - o2P) / (T2 - T1)) * (speciesTemperature - T1) + o1P
  opacityRosseland = ((o1R - o2R) / (T2 - T1)) * (speciesTemperature - T1) + o1R
!
!
!   ...Set the absorption, emission and transport opacities.
!
!
  op_speciesOpacities (ABSORPTION,species) = opacityPlanck
  op_speciesOpacities (  EMISSION,species) = opacityPlanck
  op_speciesOpacities ( TRANSPORT,species) = opacityRosseland
!
!
!   ...Ready! 
!
!
  return
end subroutine op_getSpeciesLowTempOpacities
