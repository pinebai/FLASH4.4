#include "Ionize.h"
!..Solve the ionization equilibrium problem in flash 2.0
!
!..routine neimn drives the ioniz evolution for the set of 12 elements
!..routine initioniz initializes the ioniz problem
!
!
!
subroutine neimn(tstep,tt,dd,xin,xout,sdot,idx,selnel)
!
! *****************************************************************
!
!  LAST AUTHOR UPDATE: CHICAGO AUGUST 28, 2001
!  LAST UPDATE:       CHICAGO OCTOBER 15, 2010
!
!  NAME:    NEIMN
!
!  PURPOSE: given time step TSTEP, temperature TT, density DD
!           this routine returns the ion species XOUT in the 
!           approximation of ionization equilibrium 
!           [energy generation rate SDOT - TBD].
!
!  INPUTS:
!           TSTEP    REAL      time step
!           TT       REAL      temperature at t0 = t0+tstep
!           DD       REAL      density at t0 = t0+tstep
!           XIN      REAL      ion species at t0 [dummy variable]
!           IDX      INTEGER   Elements ID vector
!           SELNEL   INTEGER   N. Elements considered
!
!  OUTPUTS:
!           XOUT     REAL      ion species at t0 = t0+tstep
!          [SDOT     REAL      energy generation rate - TBD]
!   [TT may be modified on output to force it into the allowed temperature range]
!
!  SUBROUTINES:
!           EQUIL - returns the ion species for a single element in
!                   the approximation of ionization equilibrium
!
!  REFERENCES:
!           none
!
!  MODIFICATION HISTORY:
!           WRITTEN BY: S. ORLANDO August 2001
!           MODIFIED BY: N. Flocke, K. Weide October 2010 (added intent)
!                        K. Weide            January 2012 (fixed TT intent)
!
! *****************************************************************
!
!
  use Ionize_interface, ONLY : Ionize_equil
  use Ionize_data,  ONLY : nion12 => ion_nion12
  use Driver_interface, ONLY : Driver_abortFlash
!
!
!..declare the pass
  implicit none
  real             epscheck
  parameter        (epscheck = 1.0e-4)
!
  integer, intent(in) ::  idx(ION_NELEM)

  real,    intent(in) ::  tstep
  real,    intent(inout) ::  tt
  real,    intent(in) ::  dd
  real,    intent(in) ::  xin(*)
  real,    intent(out) :: xout(*)
  real,    intent(out) :: sdot
  real,    intent(in) ::  selnel

  integer          nel,i,jj1,jj2,nion
  real             elsum,xout2(ION_NIMAX),sdot2
!
!
!
!..do the integration for each element
  jj1 = 0
  jj2 = 0
  sdot = 0
!
  do nel = 1,ION_NELEM
     if (idx(nel).eq.1) then
        nion   = nion12(nel)
        !
        call Ionize_equil(tt,nel,nion,xout2)
        sdot2 = 0.0E0
        !
        elsum = 0.0
        do i = 1,nion
           jj2 = jj2+1
           xout(jj2) = xout2(i)/selnel
           elsum = elsum+xout2(i)
        enddo
        if (abs(elsum-1.0).gt.epscheck) then
           write(*,*) ' '
           write(*,*) 'Element I.D.: ',nel,' elsum =',elsum
           write(*,*) 'error in routine neimn: '//                     &
    &                 'wrong sum of population fraction'
           call Driver_abortFlash('error in routine neimn: &
                              wrong sum of pop. fraction')
        endif
!
        sdot = sdot+sdot2
     endif
  enddo
!
!
  return
end subroutine neimn
