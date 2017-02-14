!!****f* source/Particles/Particles_getCountPerBlk
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
subroutine Particles_getCountPerBlk(perBlkCount)

  implicit none
#include "Flash.h"

  integer,dimension(MAXBLOCKS),intent(OUT) :: perBlkCount

  perBlkCount=0
  
  return
end subroutine Particles_getCountPerBlk
