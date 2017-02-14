!!****if* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_picAdvanceVel
!!
!! NAME
!!
!!  pt_picAdvanceVel
!!
!! SYNOPSIS
!!
!!  call pt_picAdvanceVel(real(in) :: dt,
!!                     logical(in) :: halfDt)
!!
!! DESCRIPTION
!!  
!!   Advance particle velocities in time by dt
!!
!! ARGUMENTS
!!
!!   dt : delta t
!!   halfDt : indicates whether to advance by half timestep
!!
!!
!!
!!***

subroutine pt_picAdvanceVel(dt, halfDt)
  
  use Particles_data, only: particles, pt_numLocal, pt_posAttrib 
  use Grid_interface, only: Grid_mapMeshToParticles

#include "Flash.h"
#include "constants.h"  
#include "Particles.h"

  implicit none
  real, intent(in)    :: dt
  logical, intent(in) :: halfDt
  integer :: numAttrib, i
  integer, dimension(2,99) :: attrib
  real :: qm, vxb(3), b(3)

  ! Map E and B to particles
  numAttrib = 9
  attrib(GRID_DS_IND, 1) = GEPX_VAR
  attrib(PART_DS_IND, 1) = EX_PART_PROP
  attrib(GRID_DS_IND, 2) = GEPY_VAR
  attrib(PART_DS_IND, 2) = EY_PART_PROP
  attrib(GRID_DS_IND, 3) = GEPZ_VAR
  attrib(PART_DS_IND, 3) = EZ_PART_PROP
  attrib(GRID_DS_IND, 4) = GRBX_VAR
  attrib(PART_DS_IND, 4) = BX_PART_PROP
  attrib(GRID_DS_IND, 5) = GRBY_VAR
  attrib(PART_DS_IND, 5) = BY_PART_PROP
  attrib(GRID_DS_IND, 6) = GRBZ_VAR
  attrib(PART_DS_IND, 6) = BZ_PART_PROP
  attrib(GRID_DS_IND, 7) = GBX1_VAR
  attrib(PART_DS_IND, 7) = BX1_PART_PROP
  attrib(GRID_DS_IND, 8) = GBY1_VAR
  attrib(PART_DS_IND, 8) = BY1_PART_PROP
  attrib(GRID_DS_IND, 9) = GBZ1_VAR
  attrib(PART_DS_IND, 9) = BZ1_PART_PROP

  call Grid_mapMeshToParticles(particles, NPART_PROPS, BLK_PART_PROP, &
       pt_numLocal, pt_posAttrib, numAttrib, attrib, WEIGHTED)

  if(halfDt) then
  ! Update velocity to time level n+1/2 (and save n)
     do i = 1, pt_numLocal
        b(1) = particles(BX_PART_PROP, i)+particles(BX1_PART_PROP, i)
        b(2) = particles(BY_PART_PROP, i)+particles(BY1_PART_PROP, i)
        b(3) = particles(BZ_PART_PROP, i)+particles(BZ1_PART_PROP, i)
        
        ! Save v^n
        particles(TMPX_PART_PROP, i) = particles(VELX_PART_PROP, i)
        particles(TMPY_PART_PROP, i) = particles(VELY_PART_PROP, i)
        particles(TMPZ_PART_PROP, i) = particles(VELZ_PART_PROP, i)
        ! compute v x B
        vxb(1) = particles(VELY_PART_PROP, i)*b(3) &
             -  particles(VELZ_PART_PROP, i)*b(2)
        vxb(2) = particles(VELZ_PART_PROP, i)*b(1) &
             -  particles(VELX_PART_PROP, i)*b(3)
        vxb(3) = particles(VELX_PART_PROP, i)*b(2) &
             -  particles(VELY_PART_PROP, i)*b(1)
        qm = particles(CHARGE_PART_PROP, i)/particles(MASS_PART_PROP, i)
        ! Advance v
        particles(VELX_PART_PROP, i) = particles(VELX_PART_PROP, i) &
             + 0.5*dt*qm*(particles(EX_PART_PROP, i)+vxb(1))
        particles(VELY_PART_PROP, i) = particles(VELY_PART_PROP, i) &
             + 0.5*dt*qm*(particles(EY_PART_PROP, i)+vxb(2))
        particles(VELZ_PART_PROP, i) = particles(VELZ_PART_PROP, i) &
             + 0.5*dt*qm*(particles(EZ_PART_PROP, i)+vxb(3))
     end do

  else
  ! Update velocity to time level n+1
     do i = 1, pt_numLocal
        b(1) = particles(BX_PART_PROP, i)+particles(BX1_PART_PROP, i)
        b(2) = particles(BY_PART_PROP, i)+particles(BY1_PART_PROP, i)
        b(3) = particles(BZ_PART_PROP, i)+particles(BZ1_PART_PROP, i)
        
        ! compute v x B
        vxb(1) = particles(VELY_PART_PROP, i)*b(3) &
             -  particles(VELZ_PART_PROP, i)*b(2)
        vxb(2) = particles(VELZ_PART_PROP, i)*b(1) &
             -  particles(VELX_PART_PROP, i)*b(3)
        vxb(3) = particles(VELX_PART_PROP, i)*b(2) &
             -  particles(VELY_PART_PROP, i)*b(1)
        qm = particles(CHARGE_PART_PROP, i)/particles(MASS_PART_PROP, i)
        ! Advance v
        particles(VELX_PART_PROP, i) = particles(TMPX_PART_PROP, i) &
             + dt*qm*(particles(EX_PART_PROP, i)+vxb(1))
        particles(VELY_PART_PROP, i) = particles(TMPY_PART_PROP, i) &
             + dt*qm*(particles(EY_PART_PROP, i)+vxb(2))
        particles(VELZ_PART_PROP, i) = particles(TMPZ_PART_PROP, i) &
             + dt*qm*(particles(EZ_PART_PROP, i)+vxb(3))
     end do
  end if


end subroutine pt_picAdvanceVel
