!!****if* source/physics/materialProperties/Opacity/localAPI/op_setLowTemperatureBoundaries
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

  implicit none

  return
end subroutine op_setLowTemperatureBoundaries
