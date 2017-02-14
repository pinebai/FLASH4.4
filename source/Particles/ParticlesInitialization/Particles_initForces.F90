!!****if* source/Particles/ParticlesInitialization/Particles_initForces
!!
!! NAME
!!
!!  Particles_initForces
!!
!! SYNOPSIS
!!
!!  call Particles_initForces()
!!
!! DESCRIPTION
!!
!!  When Active particles are being used for a simulation, this 
!!  routine will calculate the long range and short range forces
!!  caused by the gravity at initialization.
!!
!!***

subroutine Particles_initForces()
  
  use Grid_interface, ONLY: Grid_sortParticles

  use Particles_interface, ONLY : Particles_shortRangeForce, &
       Particles_longRangeForce

  use pt_interface, ONLY: pt_updateTypeDS

  use Particles_data, ONLY : particles, useParticles, &
       pt_typeInfo,pt_numLocal, pt_maxPerProc

#include "Flash.h"
#include "Particles.h"

  implicit none
  integer :: p_begin,p_end,p_count,mapType,i
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk
  logical :: needForces

  if (useParticles) then
#ifdef TYPE_PART_PROP
     call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
          pt_maxPerProc,particlesPerBlk,BLK_PART_PROP, TYPE_PART_PROP)
#else
     call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
          pt_maxPerProc,particlesPerBlk,BLK_PART_PROP)
#endif
  
     call pt_updateTypeDS(particlesPerBlk)

     do i = 1, NPART_TYPES
        needForces = (pt_typeInfo(PART_ADVMETHOD,i)==LEAPFROG)
        needForces = (pt_typeInfo(PART_ADVMETHOD,i)==LEAPFROG_COSMO).or.needForces
        needForces = (pt_typeInfo(PART_ADVMETHOD,i)==EULER_MAS).or.needForces
        if(needForces) then
           p_begin=pt_typeInfo(PART_TYPE_BEGIN,i)
           p_count=sum(pt_typeInfo(PART_LOCAL,i:NPART_TYPES))
           p_end=p_begin+p_count-1
           mapType=WEIGHTED

           call Particles_longRangeForce(particles(:,p_begin:p_end),p_count,mapType)
        
           call Particles_shortRangeForce
        end if
     end do
  end if
end subroutine Particles_initForces
