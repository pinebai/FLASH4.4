!!****if* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_picSetCurrents
!!
!! NAME
!!
!!  pt_picSetCurrents
!!
!! SYNOPSIS
!!
!!  call pt_picSetCurrents()
!!
!! DESCRIPTION
!!
!! Change fields, e.g., charge density   
!!
!! ARGUMENTS
!!
!!
!!
!!***

subroutine pt_picSetCurrents()
  ! Compute current for later deposit on grid
  use Particles_data, only: pt_numLocal, particles

#include "Flash.h"

  implicit none !! Added by fix script
  integer         :: i

  do i = 1, pt_numLocal
     particles(JX_PART_PROP, i) &
          = particles(CHARGE_PART_PROP, i)*particles(VELX_PART_PROP, i)
     particles(JY_PART_PROP, i) &
          = particles(CHARGE_PART_PROP, i)*particles(VELY_PART_PROP, i)
     particles(JZ_PART_PROP, i) &
          = particles(CHARGE_PART_PROP, i)*particles(VELZ_PART_PROP, i)
  end do
end subroutine pt_picSetCurrents
