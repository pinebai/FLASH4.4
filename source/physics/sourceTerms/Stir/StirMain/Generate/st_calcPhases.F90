!!****if* source/physics/sourceTerms/Stir/StirMain/Generate/st_calcPhases
!!
!! NAME
!!
!!  st_calcPhases
!!
!! SYNOPSIS
!!
!!  st_calcPhases()
!!
!! DESCRIPTION
!!
!!     This routine updates the stirring phases from the OU phases.
!!     It copies them over, then subtracts out the divergence part.
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***



subroutine st_calcPhases()

  use Stir_data, ONLY : st_nmodes, st_mode, st_aka, st_akb, st_OUphases

#include "Flash.h"

  implicit none

  real bjiR, bjiI, kb, kk
  integer i,j

  do i = 1, st_nmodes 
     kb = 0.
     kk = 0. 
     do j=1,NDIM
        kk = kk + st_mode(j,i)*st_mode(j,i) 
        kb = kb + st_mode(j,i)*st_OUphases(6*(i-1)+2*(j-1)+0+1)
     enddo
     do j=1,NDIM
        bjiR = st_OUphases(6*(i-1)+2*(j-1) + 0 + 1)
        bjiI = st_OUphases(6*(i-1)+2*(j-1) + 1 + 1)
        st_aka(j,i) = bjiR - st_mode(j,i)*kb/kk
        st_akb(j,i) = bjiI

     enddo
  enddo

  return
end subroutine st_calcPhases
