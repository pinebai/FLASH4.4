!!****if* source/Particles/ParticlesMain/Particles_getCountPerBlk
!!
!! NAME
!!  Particles_getCountPerBlk
!!
!! SYNOPSIS
!!
!!  Particles_getCountPerBlk(integer(OUT)  :: perBlkCount(:))
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Finds the number of particles in each block.  
!! 
!!
!! ARGUMENTS 
!!
!! perBlkCount : integer array containing number of particles on each blk.
!!
!!***

#include "Flash.h"
subroutine Particles_getCountPerBlk(perBlkCount)

  use Particles_data, ONLY : particles, pt_numLocal
#ifdef DEBUG_PARTICLES
  use Logfile_interface, ONLY: Logfile_open, Logfile_close
#endif

  implicit none
#include "constants.h"

  integer,dimension(MAXBLOCKS),intent(OUT) :: perBlkCount
  integer :: i,j
#ifdef DEBUG_PARTICLES
  integer :: logunit
#endif

  perBlkCount=0
#ifdef DEBUG_PARTICLES
  call Logfile_open(logUnit, .TRUE.)
  write(logUnit,*) 'the number of particles',pt_numLocal
  call Logfile_close(.TRUE.)
#endif
  do i = 1,pt_numLocal
     j=int(particles(BLK_PART_PROP,i))
     if(j/=NONEXISTENT) then
        perBlkCount(j)=perBlkCount(j)+1
     end if
  end do
  
  return
end subroutine Particles_getCountPerBlk
