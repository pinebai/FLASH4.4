!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/Constant/MagneticResistivity_init
!!
!! NAME
!!  MagneticResistivity_init
!!
!! SYNOPSIS
!!  MagneticResistivity_init()
!!
!! DESCRIPTION
!!  Initialize constant magnetic resistivity.
!!
!! ARGUMENTS
!!
!!***

subroutine MagneticResistivity_init()

  use MagneticResistivity_data,    ONLY : mResistivity, mUnit, &
        mag_useMagneticResistivity
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  

  ! Everybody should know these
  

  call RuntimeParameters_get("useMagneticResistivity", mag_useMagneticResistivity)
  call RuntimeParameters_get("resistivity", mResistivity)
  call RuntimeParameters_get("UnitSystem",  mUnit)

  return
end subroutine MagneticResistivity_init
