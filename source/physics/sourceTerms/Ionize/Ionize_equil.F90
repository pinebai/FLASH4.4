!!****f* source/physics/sourceTerms/Ionize/Ionize_equil
!!
!!  NAME
!!
!!   Ionize_equil
!!
!!  SYNOPSIS
!!
!!   Ionize_equil(real(INOUT)  :: tx, 
!!                integer(IN)  :: nel, 
!!                integer(OUT) :: nion, 
!!                real(OUT)    :: delem(ION_NIMAX))
!!
!!  DESCRIPTION
!!
!!   Computes the ions density under equilibrium ionization.
!!
!!  ARGUMENTS
!!
!!           tx -    REAL         Extrapolated temperature
!!           nel -    INTEGER      Element I.D.
!!                               1,  2, 3, 4, 5,  6,  7,  8, 9,  10, 11, 12
!!                               he, c, n, o, ne, mg, si, s, ar, ca, fe, ni
!!
!!           nion -    INTEGER     Number of ionization states
!!           delem -   REAL        Population fraction
!!***

subroutine Ionize_equil(tx,nel,nion,delem)
  implicit none
#include "Ionize.h"
  real, intent(INOUT) :: tx
  integer, intent(IN) :: nel
  integer, intent(OUT) :: nion
  real, dimension(ION_NIMAX), intent(OUT) :: delem
  nion = 0
  delem(1:ION_NIMAX) = 0.0
  return
end subroutine Ionize_equil
