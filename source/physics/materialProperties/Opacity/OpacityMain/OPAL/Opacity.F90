!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/Opacity
!!
!! NAME
!!
!!  Opacity
!!
!! SYNOPSIS
!!
!!  call Opacity (real    (in)  :: soln(:),
!!                integer (in)  :: ngrp,
!!                real    (out) :: opacityAbsorption,
!!                real    (out) :: opacityEmission,
!!                real    (out) :: opacityTransport )
!!
!! DESCRIPTION
!!
!!  Computes absorption, emission and transport opacities for a particular cell.
!!
!! ARGUMENTS
!!
!!   soln              : The solution vector for the cell
!!   ngrp              : The energy group number
!!   opacityAbsorption : the absorption opacity (in 1/cm)
!!   opacityEmission   : the emission opacity (in 1/cm)
!!   opacityTransport  : the transport opacity (in 1/cm)
!!
!!***
subroutine Opacity (soln, ngrp, opacityAbsorption, opacityEmission, opacityTransport)

  use Driver_interface,  ONLY : Driver_abortFlash

  use Opacity_data,      ONLY : op_useOpacity,               &
                                op_absorptionKind,           &
                                op_emissionKind,             &
                                op_transportKind,            &
                                op_ionNumberDensity,         &
                                op_opalNumHydrogenAbundances,             &
                                op_nEnergyGroups,            &
                                op_cellMassDensity,          &
                                op_cellTemperature,          &
                                op_hydrogenMassFracVar,      &
                                op_hydrogenMassFrac,         &
                                op_speciesLowTempCutoff

  use Opacity_data,       ONLY: op_emitConst,  &
                                op_absorbConst
  use op_opalData,  ONLY : OP_OPAL_LOWT,OP_OPAL_HIGHT,       &
                                op_opalMaxLowT,              &
                                op_opalTableAbundMax

  use op_interface,      ONLY : op_getOpalTableOpacity,          &
                                op_computeMultispeciesOpacities, &
                                op_computeIonNumberDensities
  implicit none
  
#include "Flash.h"  
#include "constants.h"
#include "Opacity.h"
  
  integer, intent (in)  :: ngrp
  real,    intent (out) :: opacityAbsorption
  real,    intent (out) :: opacityEmission
  real,    intent (out) :: opacityTransport

  real,    intent (in), dimension (:) :: soln

  integer :: species
  integer :: whichTable
  integer :: mfH, m

  real    :: speciesLowTempCutoff
  real    :: speciesMinTempLowTable
  real    :: speciesMinTempHighTable
  real    :: speciesMinTempROTable
  real    :: temperature
  real    :: density
  real    :: hydrogenMassFrac
!
!
!   ...Check, if opacities need to be evaluated at all.
!
!
  if (.not.op_useOpacity) then
       return
  end if
!
!
!   ...Check the energy group and abort, if out of range.
!
!
  if ( (ngrp < 1) .or. (ngrp > op_nEnergyGroups) ) then
      call Driver_abortFlash ('[Opacity] ERROR: no corresponding energy group')
  end if
!
!
!   ...Store all cell specific properties into internal data.
!
!
  op_cellMassDensity = soln (DENS_VAR)
  op_cellTemperature = soln (TELE_VAR)

  if (op_cellTemperature > op_opalMaxLowT) then
     whichTable = OP_OPAL_HIGHT
  else
     whichTable = OP_OPAL_LOWT
  end if


  if (op_hydrogenMassFracVar > 0) then
     hydrogenMassFrac = soln (op_hydrogenMassFracVar)
  else
     hydrogenMassFrac = op_hydrogenMassFrac
  end if

!
!
!   ...Generate the ion number densities for all the species in current cell.
!
!
  call op_computeIonNumberDensities ()

!!$  print*,'[Opacity] op_opalNumHydrogenAbundances,op_opalTableAbundMax=',op_opalNumHydrogenAbundances,op_opalTableAbundMax
  mfH = op_opalNumHydrogenAbundances !default
  do m = 1,op_opalNumHydrogenAbundances
     if (op_opalTableAbundMax(m) .GE. hydrogenMassFrac) then
        mfH = m
        EXIT
     end if
  end do
!!$  print*,'[Opacity] hydrogenMassFrac,m,mfH=',hydrogenMassFrac,m,mfH

  temperature    = op_cellTemperature

  density = op_cellMassDensity

  call op_getOpalTableOpacity (mfH, whichTable,            &
                                        temperature, &
                                        density,     &
                                        ngrp,               &
                                        opacityTransport         )
!
!
!   ...Calculate the total cell opacities.
!
!
  opacityAbsorption = op_cellMassDensity * op_absorbConst
  opacityEmission   = op_cellMassDensity * op_emitConst
  opacityTransport  = op_cellMassDensity * opacityTransport
!
!
!   ...Ready!
!
!
  return
end subroutine Opacity

