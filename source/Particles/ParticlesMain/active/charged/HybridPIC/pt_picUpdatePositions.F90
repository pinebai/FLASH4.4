!!****if* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_picUpdatePositions
!!
!! NAME
!!
!!  pt_picUpdatePositions
!!
!! SYNOPSIS
!!
!!  call pt_picUpdatePositions(real(in) :: dt)
!!
!! DESCRIPTION
!!
!!   Move particles by dt in the physical domain
!!
!! ARGUMENTS
!!
!!   dt : time step
!!
!!
!!
!!***

subroutine pt_picUpdatePositions(dt, halfdt)
  ! Move particles by dt
  use Particles_data, only: pt_numLocal, particles

#include "Flash.h"

  implicit none
  real, intent(in)    :: dt
  logical, intent(in) :: halfdt
  integer :: i
  

  if(halfdt) then
  ! Update position to time level n+1/2
     do i = 1, pt_numLocal
        particles(POSX_PART_PROP, i) = particles(TMPX_PART_PROP, i) 
        particles(POSY_PART_PROP, i) = particles(TMPY_PART_PROP, i) 
        particles(POSZ_PART_PROP, i) = particles(TMPZ_PART_PROP, i) 
     end do
  else
     ! Update position to time level n (and save n+1/2)
     do i = 1, pt_numLocal
        particles(TMPX_PART_PROP, i) = particles(POSX_PART_PROP, i) &
             + dt*particles(VELX_PART_PROP, i)
        particles(POSX_PART_PROP, i) = 0.5*(particles(TMPX_PART_PROP, i) &
             + particles(POSX_PART_PROP, i))
        particles(TMPY_PART_PROP, i) = particles(POSY_PART_PROP, i) &
             + dt*particles(VELY_PART_PROP, i)
        particles(POSY_PART_PROP, i) = 0.5*(particles(TMPY_PART_PROP, i) &
             + particles(POSY_PART_PROP, i))
        particles(TMPZ_PART_PROP, i) = particles(POSZ_PART_PROP, i) &
             + dt*particles(VELZ_PART_PROP, i)
        particles(POSZ_PART_PROP, i) = 0.5*(particles(TMPZ_PART_PROP, i) &
             + particles(POSZ_PART_PROP, i))
     end do
  end if
end subroutine pt_picUpdatePositions
