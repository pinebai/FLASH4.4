!!****if* source/physics/sourceTerms/Ionize/IonizeMain/Nei/nei
!!
!!  NAME
!!
!!   nei.F90
!!
!!
!!  DESCRIPTION
!!
!!..Collection of routines to solve the non-equilibrium ionization problem 
!!
!!..routine neimn drives the NEI evolution for the set of 12 elements
!!..routine neimn2 drives the NEI evolution for a single element
!!..routine dernei sets up the odes
!!..routine sjacnei sets up the sparse dernei jacobian
!!..routine bjacnei builds the nonzero locations for sjacnei
!!..routine initioniz initializes the NEI problem
!!
!!    
!! 
!!***

#include "Ionize.h"
!
!
subroutine neimn(tstep, tt, dd, xin, xout, sdot, idx, xfrac)
!
! *****************************************************************
!
!  LAST UPDATE: CHICAGO AUGUST 28, 2001
!
!  NAME:    NEIMN2
!
!  PURPOSE: given time step TSTEP, temperature TT, density DD, and
!           ion species XIN, this routine returns the evolved
!           ion species XOUT [and the energy generation rate SDOT - TBD].
!
!  INPUTS:
!           TSTEP    REAL    time step
!           TT       REAL    temperature (in K) at t0 = t0+tstep
!           DD       REAL    density (in particles/cm^3) at t0 = t0+tstep
!           XIN      REAL    ion species at t0
!           IDX      INTEGER   Elements ID vector
!
!  OUTPUTS:
!           XOUT     REAL    evolved ion species at t0 = t0+tstep
!          [SDOT     REAL    energy generation rate - TBD]
!   [TT may be modified on output to force it into the allowed temperature range]
!
!  SUBROUTINES:
!           NEIMN2 - returns the evolved ion species for a single element
!
!  NOTE:    The sum of ALL the fluids abundances (the population
!           fractions) has to be 1 in the Flash code 
!           [see the real NELSEL]
!
!  REFERENCES:
!           none
!
!  MODIFICATION HISTORY:
!           WRITTEN BY: S. ORLANDO August 2001
!           MODIFIED BY: K. Weide  January 2012 (fixed TT intent)
!
! *****************************************************************
!
!
!
  use Ionize_data, ONLY : nion12 => ion_nion12, &
       ifirst => ion_ifirst, ion_ELEC, ion_HYD
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  
  integer, intent(in) :: idx(ION_NELEM)
  real, intent(in) :: tstep,dd,xin(*),xfrac(*)
  real, intent(inout) :: tt
  real, intent(out) :: xout(*),sdot
  
  real, parameter :: epscheck = 1.0e-4
  
  integer :: nel,i,jj1,jj2,nion
  real :: elsum,xin2(ION_NIMAX),xout2(ION_NIMAX),sdot2
  
  
  !
  !..return the hydrogen and electron populations without changes
  xout(ion_HYD) = xin(ion_HYD)
  xout(ion_ELEC) = xin(ion_ELEC)
  !
  !..do the integration for each element
  jj1  = 2
  jj2  = 2
  sdot = 0.0
  !
  do nel = 1,ION_NELEM

     !idx array selects which ions we use.  

     if (idx(nel).eq.1) then 
        !
        !..initialize ifirst (used in stifbs4) to zero
        ifirst = 0
        !
        nion = nion12(nel)
        do i = 1,nion
           jj1 = jj1+1
           xin2(i) = xin(jj1)/xfrac(jj1)
        enddo
        !
        call neimn2(tstep,tt,dd,nel,xin2,xout2,sdot2)
        !
        elsum = 0.0
        do i = 1,nion
           jj2 = jj2+1
           xout(jj2) = xout2(i)*xfrac(jj2)
           elsum = elsum+xout2(i)
        enddo
        !
        if (abs(elsum-1.0).gt.epscheck) then
           write(*,*) ' '
           write(*,*) 'Element I.D.: ',nel,' elsum =',elsum
           write(*,*) 'error in routine neimn: '//                     &
                &                 'wrong sum of population fraction'
           call Driver_abortFlash('error in routine neimn: wrong sum of pop. fraction')
        endif
        !
        sdot = sdot+sdot2
     endif
  enddo
  !
  !
  return
end subroutine neimn
!
!
!
!
!
!
!
!
subroutine neimn2(tstep,tt,dd,nel,xin,xout,sdot)
  !
! *****************************************************************
!
!  LAST UPDATE: CHICAGO AUGUST 27, 2001
!
!  NAME:    NEIMN2
!
!  PURPOSE: given time step TSTEP, temperature TT, density DD, and
!           ion species XIN, this routine returns the evolved
!           ion species XOUT [and the energy generation rate SDOT - TBD].
!
!  INPUTS:
!           TSTEP    REAL      time step
!           TT       REAL      temperature at t0 = t0+tstep
!           DD       REAL      density at t0 = t0+tstep
!           NEL      INTEGER   Element I.D.
!                              1,  2, 3, 4, 5,  6,  7,  8, 9,  10, 11, 12
!                              he, c, n, o, ne, mg, si, s, ar, ca, fe, ni
!           XIN      REAL      ion species at t0
!
!  OUTPUTS:
!           XOUT     REAL      evolved ion species at t0 = t0+tstep
!          [SDOT     REAL      energy generation rate - TBD]
!
!  SUBROUTINES:
!           NEIINT - ODE integrator for stiff ODEs 
!                    with an analytic and sparse jacobian
!
!  REFERENCES:
!           none
!
!  MODIFICATION HISTORY:
!           WRITTEN BY: S. ORLANDO August 2001
!
! *****************************************************************
!
!
  use ion_interface, ONLY : ion_intCoeff
  use Ionize_data, ONLY : btemp => ion_btemp, den => ion_den, &
       nion => ion_nion, cz => ion_cz ,al => ion_al
!
!..declare the pass
  implicit none
  integer, intent(in) :: nel
  real, intent(in) :: tstep,dd,xin(*)
  real, intent(inout) :: tt
  real, intent(out) :: xout(*),sdot
!
!..local variables
  external         dernei,sjacnei,bjacnei,stifbs4
!
  integer          nostore,tdim,nok,nbad,kount,iprint,i,k,nin
  parameter        (nostore=0, tdim=10, iprint=0)
  !
  real             beg,stptry,stpmin,ys2(ION_NIMAX),tol,ttime(tdim),    &
       &                 elem(ION_NIMAX,tdim),odescal,avo,ev2erg,conv
  parameter        (beg     = 0.0e0,      tol  = 1.0e-5,            &
       &                  odescal = 1.0e-6,     avo  = 6.0221367e23,      &
       &                  ev2erg  = 1.602e-12,  conv = ev2erg*1.0e6*avo)
  integer, parameter :: nimax = ION_NIMAX

!
!..get the ionization and recombination coefficients
  call ion_intCoeff(nel,tt,al,cz,nin)
!
!..set the ion and plasma variables
  nion  = nin
  btemp = tt
  den   = dd
  do i = 1,nion
     ys2(i) = xin(i)
  enddo
!
!..set the time step variables for a single point 
  stptry = tstep
  stpmin = stptry * 1.0e-20
!
!..bader-deuflhard integration
!..with ma28 algebra
!
!..Inputs:
!..   beg     = beginning integration point
!..   stptry  = suggested first step size
!..   stpmin  = minimum allowable step size
!..   tstep   = ending integration point
!..   ys2     = initial conditions
!..   tol     = desired fraction error during the integration
!..   beg[2]  = incremental value of indep var. at which to store the sol.
!..             if zero -> sol. is stored at every intermediate point
!..   nostore = maximum number of sol. points to store
!..   odescal = error scaling factor
!..   iprint  = determines if the sol. is printed as it evolves
!..   dernei  = name of the routine that contains the odes
!..   sjacnei  = name of the rout. that contains the jacobian of the odes
!..   bjacnei = name of the rout. that sets the pointers of the sparse jacobian
!..   stifbs4 = name of the rout. that will take a single step
!
!..Outputs:
!..   ttime   = indep.var.sol.array of phys.dim.(TDIM), log.dim.(TDIM)
!..   elem    = dep.var.sol.ar. of phys.dim.(NIMAX,TDIM), log.dim.(NION,TDIM)
!..   nok     = number of succesful steps taken
!..   nbad    = number of bad steps taken, bad but retried and then succesful
!..   kount   = total number of steps stored in arrays ttime and elem
!..Note: At the end of the computation, ys2 is updated
!
  call neiint(beg,stptry,stpmin,tstep,ys2,                          &
     &            tol,beg,nostore,                                      &
     &            ttime,elem,tdim,nimax,tdim,nion,                      &
     &            nok,nbad,kount,odescal,iprint,                        &
     &            dernei,sjacnei,bjacnei,stifbs4)
!
!..the average energy generated over the time step (TBD)
  sdot = 0.0e0
  !      do k=1,nion
  !       sdot = sdot + (ys2(k) - ymass(k)) * bion(k)
  !      enddo
  !      sdot = sdot * conv/tstep
  !
  !..update the composition
  do i = 1,nion
     xout(i) = ys2(i)
  enddo
  !
!
  return
end subroutine neimn2
!
!
!
!
!
!
subroutine dernei(tt,y,dydt)
!
! *****************************************************************
!
!  LAST UPDATE: CHICAGO AUGUST 24, 2001 
!
!  NAME:    DERNEI
!
!  PURPOSE: given time step TT, Population fraction Y, this routine 
!           sets up the system of ODE's for the NEI problem.
!  
!  NOTE:    the source coefficients are adequate to optically thin 
!           plasma in the coronal approximation
!
!  INPUTS:
!           TT      REAL      time step
!           Y       REAL      Population fraction
!
!  OUTPUTS:
!           DYDT    REAL      right hand side of the ODE's
!
!  SUBROUTINES:
!           none
!
!  REFERENCES:
!           none
!
!  MODIFICATION HISTORY:
!           WRITTEN BY: S. ORLANDO August 2001
!
! *****************************************************************
!
!

  use Ionize_data, ONLY : den => ion_den, &
       nion => ion_nion, cz => ion_cz, al => ion_al

!
!..declare
  implicit none
  real, intent(in) :: tt,y(*)
  real, intent(out) :: dydt(*)
  
  real :: pop,depop
  integer :: i
  !
  !..set up the system of odes:
  do i = 1,nion
     pop = 0.0
     if (i.gt.1) then 
        pop = y(i-1)*cz(i-1)
     endif
     if (i.lt.nion) then 
        pop = pop+y(i+1)*al(i+1)
     endif
     depop   = y(i)*(cz(i)+al(i))
     dydt(i) = den*(pop-depop)
  end do
  !
  !
  return
end subroutine dernei
!
!
!
!
! *****************************************************************
!
subroutine sjacnei(tt,y,dfdy,nzo)
!
! *****************************************************************
!
!  LAST UPDATE: CHICAGO AUGUST 24, 2001 
!
!  NAME:    SJACNEI
!
!  PURPOSE: this routine sets up the sparse DERNEI jacobian.
!
!  INPUTS:
!           TT      REAL      time step (irrelevant here)
!           Y       REAL      Population fraction
!           NZO     INTEGER   matrix element location  [DFDY(NZO)]
!
!  OUTPUTS:
!           DFDY    REAL      jacobian
!
!
!  MODIFICATION HISTORY:
!           WRITTEN BY: S. ORLANDO August 2001
!
! *****************************************************************
!
!
  use Ionize_data, ONLY : den => ion_den, &
       nion => ion_nion, cz => ion_cz ,al => ion_al
  use Driver_interface, ONLY : Driver_abortFlash

!
!..declare
  implicit none
  real, intent(in) :: tt,y(*)
  real, intent(out) :: dfdy(*)
  integer, intent(in) :: nzo
  
  integer :: i,k,nt,iat,ik1,ik2
  real :: aa(3)
  !
  !..communicate with the jacobian builder BJACNEI
  integer          neloc
  parameter        (neloc = ION_NIMAX*3)   
  integer          eloc(neloc),nterms
  common /locnei/  eloc,nterms
  
  !
  !..initialize
  nt = 0
  do i = 1,nzo
     dfdy(i) = 0.0e0
  enddo
  !
  !..jacobian elements
  do i = 1,nion
     aa(1) =  den * cz(i-1)
     aa(2) = -den * (al(i)+cz(i))
     aa(3) =  den * al(i+1)
     if (i.eq.1) then
        ik1 = 2
     else 
        ik1 = 1
     endif
     if (i.eq.nion) then
        ik2 = 2
     else 
        ik2 = 3
     endif
     do k = ik1,ik2
        nt  = nt + 1
        iat = eloc(nt)
        dfdy(iat) = aa(k)
     enddo
  enddo
!
!..bullet check the counting
  if (nt .ne. nterms) then
     write(6,*) 'nt =',nt,'  nterms =',nterms
     write(6,*) 'error in routine sjacnei: nt .ne. nterms'
     call Driver_abortFlash('error in routine sjacnei: nt .ne. nterms')
  end if
!
!
  return
end subroutine sjacnei
!
!
!
!
! *****************************************************************
!
subroutine bjacnei(iloc,jloc,nzo,np)
!
! *****************************************************************
!
!  LAST UPDATE: CHICAGO AUGUST 24, 2001 
!
!  NAME:    BJACNEI
!
!  PURPOSE: this routine builds the nonzero matrix locations for SJACNEI
!
!  INPUTS:
!           ILOC    INTEGER   matrix element location [I]
!           JLOC    INTEGER   matrix element location [J]
!           NP      INTEGER   dimension of ILOC and JLOC
!
!  OUTPUTS:
!           ILOC    INTEGER   matrix element location [I]
!           JLOC    INTEGER   matrix element location [J]
!           NZO     INTEGER   number of matrix element location
!
!  SUBROUTINES:
!           MCORD - marks the locations of a sparse matrix, 
!                   and the contributions to the matrix element.
!
!  MODIFICATION HISTORY:
!           WRITTEN BY: S. ORLANDO August 2001
!
! *****************************************************************
!
!
  use ion_interface, ONLY : ion_mcord
  use Ionize_data, ONLY : nion => ion_nion

!
!..declare
  implicit none
  integer, intent(in) :: np
  integer, intent(inout) :: iloc(*),jloc(*)
  integer, intent(out) :: nzo
  
  integer :: i,j
!
!..communicate with SJACNEI
  integer          neloc
  parameter        (neloc = ION_NIMAX*3)
  integer          eloc(neloc),nterms
  common /locnei/  eloc,nterms

!
!..initialize
  nterms = 0
  nzo    = 0
  do i = 1,neloc
     eloc(i) = 0
  enddo
!
!..tag the nonzero locations

  do i = 1,nion
     if (i.gt.1)                                                     &
          &   call ion_mcord(i,i-1,iloc,jloc,nzo,np,eloc,nterms,neloc)

     call ion_mcord(i,i,iloc,jloc,nzo,np,eloc,nterms,neloc)

     if (i.lt.nion)                                                  &
          &   call ion_mcord(i,i+1,iloc,jloc,nzo,np,eloc,nterms,neloc)
  enddo

!
!..
!..write a diagnostic
!      write(6,*) ' '
!      write(6,*) nzo,' matrix elements'
!      write(6,*) nterms,' jacobian contributions'
!      write(6,*) ' '
!
!
  return
end subroutine bjacnei
