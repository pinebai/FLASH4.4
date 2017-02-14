!!****if* source/physics/sourceTerms/Ionize/IonizeMain/Ionize_equil
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
!!   INPUTS
!!           tx   -  REAL         Extrapolated temperature
!!           nel  -  INTEGER      Element I.D.
!!                               1,  2, 3, 4, 5,  6,  7,  8, 9,  10, 11, 12
!!                               he, c, n, o, ne, mg, si, s, ar, ca, fe, ni
!!
!!   OUTPUTS
!!           nion  -   INTEGER     Number of ionization states
!!           delem -   REAL        Population fraction
!!    [tx may be modified on output to force it into the allowed temperature range]
!!***

subroutine Ionize_equil(tx,nel,nion,delem)
  use ion_interface, ONLY : ion_intCoeff

  !
  !..declare
  implicit none
#include "Ionize.h"
  real, intent(INOUT) :: tx
  integer, intent(IN) :: nel
  integer, intent(OUT) :: nion
  real, dimension(ION_NIMAX), intent(OUT) :: delem
  
  real, dimension(0:ION_NIMAX) :: cz
  real, dimension(ION_NIMAX) :: al
  integer :: i,j,nw
  real :: cc

  !
  !---------------------------------------------------------------------
  !
  if (tx.ge.1.e4) then
     call ion_intCoeff(nel,tx,al,cz,nion)
     nw = nion-1
     cc = 1.
     do i = nw,1,-1
        cc = cz(i)/al(i+1)*cc+1.
     end do
     delem(1) = 1./cc
     do i = 1,nw
        delem(i+1) = cz(i)/al(i+1)*delem(i)
     end do
  else
     delem(1) = 1.
     do i = 2,ION_NIMAX
        delem(i) = 0.
     end do
  end if
  !
  !
  return
end subroutine Ionize_equil
