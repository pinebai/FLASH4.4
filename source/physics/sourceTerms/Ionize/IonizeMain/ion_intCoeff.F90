!!****if* source/physics/sourceTerms/Ionize/IonizeMain/ion_intCoeff
!!
!! NAME
!!  
!!  ion_intCoeff
!!
!!
!! SYNOPSIS
!! 
!!  call ion_intCoeff( integer(IN)   :: nel,
!!                         real(INOUT)  :: tx,
!!                        real(OUT)  :: al,
!!                        real(OUT)  :: cz
!!                     integer(OUT)  :: nion)
!!
!!  
!! DESCRIPTION
!!
!!  This entry interpolates the ionization and recombination
!!           coefficients to the hydro temperatures. The tables of the
!!           coefficients are those of Summers stored in the current 
!!           directory for the Ionize unit.
!!
!!
!! ARGUMENTS
!!
!!  INPUTS
!!           nel   --            Element I.D.
!!                               1,  2, 3, 4, 5,  6,  7,  8, 9,  10, 11, 12
!!                               he, c, n, o, ne, mg, si, s, ar, ca, fe, ni
!!           tx    --            Temperature
!!
!!  OUTPUTS
!!           al    --            Recombination Coefficients
!!           cz    --            Ionization coefficients
!!           nion  --            Number of ionic states
!!
!!  NOTES: 
!!
!!           WRITTEN BY: S. ORLANDO AUGUST 2001
!!
!!***

subroutine ion_intCoeff(nel,tx,al,cz,nion)
#include "Ionize.h"
  use ion_interface, only : ion_wint

  !The following module level variables are not modified
  !in this subroutine:
  use Ionize_data, only : ion_tp,ion_cfinz,ion_cfric,ion_nion12
  !
  !..declare
  implicit none
  integer, intent(IN) :: nel
  real, intent(INOUT) :: tx
  real, dimension(ION_NIMAX), intent(OUT) :: al
  real, dimension(0:ION_NIMAX), intent(OUT) :: cz
  integer, intent(OUT) :: nion
  integer :: j,i
  real :: ang_cf
  integer, parameter :: ntemp = ION_NTEMP

  !..interpolate at temperature tx
  tx = max(tx,1.e4)
  tx = min(tx,3.e7)

  call ion_wint(ion_tp,ntemp,tx,j)

  ang_cf = (tx-ion_tp(j))/(ion_tp(j+1)-ion_tp(j))

  nion = ion_nion12(nel)
  do i = 1,nion+1
     al(i)   = ang_cf*(ion_cfric(nel,j+1,i)-ion_cfric(nel,j,i))+             &
          &            ion_cfric(nel,j,i)
     cz(i-1) = ang_cf*(ion_cfinz(nel,j+1,i-1)-ion_cfinz(nel,j,i-1))+         &
          &            ion_cfinz(nel,j,i-1)
  end do
  
  return
end subroutine ion_intCoeff
