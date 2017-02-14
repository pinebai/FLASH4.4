!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_setLowTemperatureBoundaries
!!
!! NAME
!!
!!  op_setLowTemperatureBoundaries
!!
!! SYNOPSIS
!!
!!  call op_setLowTemperatureBoundaries ()
!!
!! DESCRIPTION
!!
!!  This routine sets the temperature boundaries (min and max) for the low temperature opacity
!!  unit. It considers the lowest temperatures for each tabulated opacities for each species
!!  and the low temperature cutoff value for each species provided by the user. The maximum
!!  boundary is set such that it ensures no undefined temperature sections for any species.
!!
!! ARGUMENTS
!!
!!***
subroutine op_setLowTemperatureBoundaries ()

  use Opacity_data,         ONLY : op_totalSpecies,              &
                                   op_absoluteLowestTemperature, &
                                   op_speciesLowTempCutoff

  use op_lowTempData,       ONLY : op_minLowTempTemperature, &
                                   op_maxLowTempTemperature

  use op_tabulatedData,     ONLY : op_initializedTabulated,  &
                                   op_speciesMinTempPATable, &
                                   op_speciesMinTempPETable, &
                                   op_speciesMinTempROTable

  use Driver_interface,     ONLY : Driver_abortFlash

  implicit none

  real    :: maxLowestTabulatedTemp
  real    :: maxLowestTabulatedTempPA
  real    :: maxLowestTabulatedTempPE
  real    :: maxLowestTabulatedTempRO
  real    :: maxLowTempCutoff
!
!
!    ...Check first, if the tabulated opacity section has been initialized. Otherwise
!       we have no info about the tabulated temperature grid and the calculation must be
!       halted. It is a programming error in the sense that this routine was called
!       in the main opacity unit initialization routine before inititalization of the
!       tabulated opacity section, and it should be an easy fix.
!
!
  if (.not.op_initializedTabulated) then
       call Driver_abortFlash ('[op_setLowTemperatureBoundaries] PROGRAMMING ERROR: Tabulated not initialized')
  end if
!
!
!    ...Set the lower low temperature bound as the lowest representable temperature of the
!       opacity unit.
!
!
  op_minLowTempTemperature = op_absoluteLowestTemperature
!
!
!    ...Form the maximum of the all the lowest tabulated grid temperatures.
!
!
  maxLowestTabulatedTempPA = maxval (op_speciesMinTempPATable)
  maxLowestTabulatedTempPE = maxval (op_speciesMinTempPETable)
  maxLowestTabulatedTempRO = maxval (op_speciesMinTempROTable)

  maxLowestTabulatedTemp   = max (maxLowestTabulatedTempPA, &
                                  maxLowestTabulatedTempPE, &
                                  maxLowestTabulatedTempRO)
!
!
!    ...Form the maximum of the all the low temperature cutoffs.
!
!
  maxLowTempCutoff = maxval (op_speciesLowTempCutoff)
!
!
!    ...Set the upper low temperature bound as the maximum between both previous maxima.
!
!
  op_maxLowTempTemperature = max (maxLowestTabulatedTemp , maxLowTempCutoff)
!
!
!    ...Ready!
!
!
  return
end subroutine op_setLowTemperatureBoundaries
