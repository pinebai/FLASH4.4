!!****f* source/physics/materialProperties/MagneticResistivity/MagneticResistivity
!!
!! NAME
!!  MagneticResistivity
!!
!! SYNOPSIS
!!  MagneticResistivity(real,    intent(IN)  :: temp,
!!                      real,    intent(IN)  :: dens,
!!                      real,    intent(IN)  :: xn(NSPECIES),
!!                      real,    intent(OUT) :: magResist)
!!
!! DESCRIPTION
!!  A stub implementation for MagneticResistivity.
!!  Returns magResist = 0.0
!!
!! ARGUMENTS
!!  temp      - Plasma temperature
!!  dens      - Plasma density
!!  xn        - Species
!!  magResist - Magnetic resistivity
!!
!!***

subroutine MagneticResistivity(temp,dens,xn,magResist)

#include "Flash.h"

  implicit none

  real, intent(IN) :: temp, dens
  real, intent(IN), dimension(NSPECIES) :: xn
  real, intent(OUT):: magResist

  magResist = 0.0

  return
end subroutine MagneticResistivity
