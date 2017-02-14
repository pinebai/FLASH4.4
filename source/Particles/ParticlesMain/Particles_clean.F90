!!****if* source/Particles/ParticlesMain/Particles_clean
!!
!! NAME
!!
!!  Particles_clean
!!
!! SYNOPSIS
!!
!!  Particles_clean(real(inout) :: particles(NPART_PROP,part_num) 
!!                  real(inout) :: part_num)
!!
!! DESCRIPTION
!!
!!  Particles cleaning routine for the particle module.
!!  cleans the virtual particles after the forcing has been done 
!!  
!! ARGUMENTS
!!
!!   particles(NPART_PROP,part_num) 
!!   part_num`
!!  
!!
!!***

!===============================================================================

subroutine Particles_clean()

  use Particles_data, ONLY : particles, pt_numLocal
  implicit none
  
#include "constants.h"  
#include "Flash.h"
  
  ! local variables
  integer :: i,j,local_num


  local_num=pt_numLocal
 
  j=1
  do i=1,local_num
     if(particles(TAG_PART_PROP,j)<0.) then
        particles(:,j)=particles(:,pt_numLocal)
        pt_numLocal=pt_numLocal-1
     else
        j=j+1
     end if
  end do
#ifdef DEBUG_VPARTICLES
  write(*,*) 'Element=',i,'TAG=',int(particles(TAG_PART_PROP,i)),'Blk=',int(particles(BLK_PART_PROP,i)),particles(POSX_PART_PROP,i),particles(POSY_PART_PROP,i)
  print *, 'This particles is a virtual particle and will be deleted'
  write(*,*)'Initial local no.=',local_num,i,'Current local no.=',part_num
#endif
  
end subroutine Particles_clean
