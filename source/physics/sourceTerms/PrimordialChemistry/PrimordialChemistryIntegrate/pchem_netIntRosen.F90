!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryIntegrate/pchem_netIntRosen
!!
!! NAME
!!   pchem_rosenMa28  (pchem_rosenGift is below)
!!
!! SYNOPSIS
!!  subroutine pchem_rosenMa28(real(IN)       ::y(:),
!!                          real(IN)       ::dydx(:),
!!                          integer(IN)    ::n,
!!                          real(IN)       ::x,
!!                          real(IN)       ::htry,
!!                          real(IN)       ::eps,
!!                          real(IN)       ::yscal(:),
!!                          real(OUT)      ::hdid,
!!                          real(OUT)      ::hnext, 
!!                          procedure(IN)  :: derivs,
!!                          procedure(IN)  :: jakob,
!!                          procedure(IN)  :: bjakob) 
!!
!! DESCRIPTION
!!
!!  fourth order rosenbrock step for integrating stiff ode's with monitoring 
!!  of local truncation error to adjust the stepsize. input are the dependent  
!!  variable y(1:n) and its derivative dydx(1:n), at the starting value of the  
!!  independent variable x. also input are the stepsize to be attempted htry,  
!!  the desired fractional accuracy eps, and the vector yscal(1:n) against  
!!  which the error is scaled. derivs is a routine that computes the righthand  
!!  side of the system of ode's. on output the y and x are replaced by their  
!!  new values, hdid is the stepsize that was actually accomplished and hnext  
!!  is the estimate for the next stepsize. 
!!   
!!  nmax is the maximum number of ode's.  
!!  grow and shrnk are the extremes of the time step growth  
!!  errcon = (grow/safety)**(1/pgrow) handles the case when errmax = 0.0 
!!  
!!  routine modified to assume dfdx = 0.0 for reaction networks
!!  
!!  sparse ma28 pchem_algebra version
!!
!! ARGUMENTS
!!
!!   y       - dependent variable, array of size y(1:n)
!!   dydx    - derivative of dependent variable, array of size dydx(1:n)
!!   n       - number of dependent variables
!!   x       - independent variable
!!   htry    - attempted stepsize 
!!   eps     - desired fractional accuracy
!!   yscal   - vector of size yscal(1:n) for scaling error
!!   hdid    - stepsize used
!!   hnext   - estimate for next stepsize
!!   derivs  - procedure(IN) name of the routine that contains the odes
!!   jakob   - procedure(IN) name of the routine that contains the jacobian of the odes
!!   bjakob  - procedure(IN) name of the routine that sets the pointers of the sparse jacobian
!!
!! NOTES
!!  this file contains 2 routines to drive the integration of 
!!  nuclear reaction networks with rosenbrock time stepper
!!  and 2 choices for the linear pchem_algebra.
!!  routine pchem_rosenMa28 drives a rosenbrock step with ma28 pchem_algebra
!!  routine pchem_rosenGift drives a rosenbrock step with gift pchem_algebra
!!
!!***
!!---------------------------------------------------------------------------------



!!---------------------------------------------------------------------------------

!!
!! NAME
!!   pchem_rosenGift
!!
!! SYNOPSIS
!!  subroutine pchem_rosenGift(real(INOUT)    ::y(:),
!!                          real(INOUT)    ::dydx(:),
!!                          integer(IN)    ::n,
!!                          real(INOUT)    ::x,
!!                          real(IN)       ::htry,
!!                          real(IN)       ::eps,
!!                          real(IN)       ::yscal(:),
!!                          real(OUT)      ::hdid,
!!                          real(OUT)      ::hnext, 
!!                          procedure(IN)  :: derivs,
!!                          procedure(IN)  :: jakob,
!!                          procedure(IN)  :: bjakob) 
!!
!! DESCRIPTION
!!
!!  fourth order rosenbrock step for integrating stiff ode's with monitoring 
!!  of local truncation error to adjust the stepsize. input are the dependent  
!!  variable y(1:n) and its derivative dydx(1:n), at the starting value of the  
!!  independent variable x. also input are the stepsize to be attempted htry,  
!!  the desired fractional accuracy eps, and the vector yscal(1:n) against  
!!  which the error is scaled. derivs is a routine that computes the righthand  
!!  side of the system of ode's. on output the y and x are replaced by their  
!!  new values, hdid is the stepsize that was actually accomplished and hnext  
!!  is the estimate for the next stepsize. 
!!   
!!  nmax is the maximum number of ode's.  
!!  grow and shrnk are the extremes of the time step growth  
!!  errcon = (grow/safety)**(1/pgrow) handles the case when errmax = 0.0 
!!  
!!  routine modified to assume dfdx = 0.0 for reaction networks
!!  
!!  sparse gift pchem_algebra version
!!
!! ARGUMENTS
!!
!!   y       - dependent variable, array of size y(1:n)
!!   dydx    - derivative of dependent variable, array of size dydx(1:n)
!!   n       - number of dependent variables
!!   x       - independent variable
!!   htry    - attempted stepsize 
!!   eps     - desired fractional accuracy
!!   yscal   - vector of size yscal(1:n) for scaling error
!!   hdid    - stepsize used
!!   hnext   - estimate for next stepsize
!!   derivs  - procedure(IN) name of the routine that contains the odes
!!   jakob   - procedure(IN) name of the routine that contains the jacobian of the odes
!!   bjakob  - procedure(IN) name of the routine that sets the pointers of the sparse jacobian
!!             Not used in this routine
!!
!! NOTES
!!  this file contains 2 routines to drive the integration of 
!!  nuclear reaction networks with rosenbrock time stepper
!!  and 2 choices for the linear pchem_algebra.
!!  routine pchem_rosenMa28 drives a rosenbrock step with ma28 pchem_algebra
!!  routine pchem_rosenGift drives a rosenbrock step with gift pchem_algebra
!!
!!***

subroutine pchem_rosenGift(y,dydx,n,x,htry,eps,yscal,hdid,hnext,  & 
     &                      derivs,jakob,bjakob,jcounts) 

  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY: Logfile_stampMessage
  use pchem_interface, ONLY: pchem_gift
  use PrimordialChemistry_data
!  use pchem_integrateInterface, ONLY: derivs

  implicit none
  
#include "constants.h"

!! argument declarations
  ! Note that bjakob is not used in this routine, and jakob explanation found in 
  ! bnNetwork_interface
  external               derivs, jakob, bjakob
  integer, intent(IN) :: n
  real, intent(IN)    :: yscal(n), dydx(n), htry, eps
  real, intent(INOUT) :: x, y(n), jcounts
  real, intent(INOUT)   :: hdid, hnext


!! local variables
  integer, parameter  :: nmax=30, maxtry=400
  real, save          :: errmax,h,xsav,dysav(nmax),err(nmax),g1(nmax),  & 
       &                 g2(nmax),g3(nmax),g4(nmax),ysav(nmax),xx
  real, parameter     ::    safety=0.9e0, grow=1.5e0, pgrow=-0.25e0,  & 
       &           shrnk=0.5e0,  pshrnk=-1.0e0/3.0e0, errcon=0.1296e0 

  integer, save       :: i,j,jtry
!!  shampine parameter set 
  real, parameter     ::  gam =  1.0e0/2.0e0,    a21 =  2.0e0,  & 
       &                  a31 =  48.0e0/25.0e0,  a32 =  6.0e0/25.0e0,  & 
       &                  c21 = -8.0e0,          c31 =  372.0e0/25.0e0,  & 
       &                  c32 =  12.0e0/5.0e0,   c41 = -112.0e0/125.0e0,  & 
       &                  c42 = -54.0e0/125.0e0, c43 = -2.0e0/5.0e0,  & 
       &                  b1  =  19.0e0/9.0e0,   b2  =  1.0e0/2.0e0,  & 
       &                  b3  =  25.0e0/108.0e0, b4  =  125.0e0/108.0e0,  & 
       &                  e1  =  17.0e0/54.0e0,  e2  =  7.0e0/36.0e0,  & 
       &                  e3  =  0.0e0,          e4  =  125.0e0/108.0e0,  & 
       &                  c1x =  1.0e0/2.0e0,    c2x = -3.0e0/2.0e0,  & 
       &                  c3x =  121.0e0/50.0e0, c4x =  29.0e0/250.0e0,  & 
       &                  a2x =  1.0e0,          a3x = 3.0e0/5.0e0


!!  for the gift linear pchem_algebra
  integer, parameter   :: nmaxp1 = nmax + 1
  real, save           :: dfdy(nmax,nmax), dmat(nmax,nmax)
! , & 
!       &                 av(nmax,nmaxp1)
  real, save           :: av(nmax,nmax), bb(nmax), yysav(nmax)
  real                 :: d
  integer              :: indx(nmax)
!!-------------------------------------------------------------------------------

!!  store the initial values 
  xsav = x 
  do i=1,n 
     ysav(i)  = y(i) 
!     print *, 'YSAV(',i,')= ', ysav(i)
     yysav(i) = y(i)
     dysav(i) = dydx(i) 
  enddo

! print *, 'N is ', n



!!  get the dense jacobian in sparse_dfdy 
  call jakob(xsav,ysav,dfdy,n,nmax)


!!  main loop 
  h = htry 
  do jtry = 1,maxtry 


   !!   form the a matrix
     xx = 1.0e0/(gam*h)

     do j=1,n
        do i=1,n
           dmat(i,j) = -dfdy(i,j) 
        enddo
     enddo
     do i=1,n
        dmat(i,i) = xx + dmat(i,i) 
     end do

   !!  set up and solve the right hand side for g1 
     do i=1,n 
        g1(i) = dysav(i) 
     enddo


     do j = 1, n
        do i = 1, n
           ! i is the row
           ! j is the column
           av(i,j) = dmat(i,j)
        enddo
     enddo
     do i = 1, n
!        av(i,n+1) = g1(i)
        bb(i) = g1(i)
     enddo
     ! here we solve AX = B (where B is column n+1)
     ! it is returning the answer also in column n+1
!     call ludcmp(av,n-1,nmax,indx,d)

!     call pchem_gift(av,nmax,nmaxp1)
     call leqs(av,bb,n-1,nmax)
!     call ludcmp(av,n-1,nmax,indx,d)
!     call lubksb(av,n-1,nmax,indx,bb)
     do i = 1, n
        g1(i) = bb(i)  !!av(i,n+1)
     enddo
        g1(iELEC) = g1(iHP) + g1(iDP) + 2.0*g1(iHEPP) + g1(iHEP) + g1(iHDP) + g1(iH2P) - g1(iHM) - g1(iDM)


   !!  compute intermediate values of y,x and dydx 
     do i=1,n 
        y(i) = ysav(i) + a21 * g1(i) 
     enddo
     x = xsav + a2x * h 
     call derivs(x,y,dydx) 


   !!  set up and solve the right hand side for g2 i
     do i=1,n 
        g2(i) = dydx(i) + c21*g1(i)/h 
     enddo


     do j = 1, n
        do i = 1, n
           av(i,j) = dmat(i,j)
        enddo
     enddo
     do i = 1, n
!        av(i,n+1) = g2(i)
        bb(i) = g2(i)
     enddo
!     call pchem_gift(av,nmax,nmaxp1)
     call leqs(av,bb,n-1,nmax)
!     call ludcmp(av,n-1,nmax,indx,d)
!     call lubksb(av,n-1,nmax,indx,bb)
     do i = 1, n
        g2(i) = bb(i)  !av(i,n+1) 
     enddo
        g2(iELEC) = g2(iHP) + g2(iDP) + 2.0*g2(iHEPP) + g2(iHEP) + g2(iHDP) + g2(iH2P) - g2(iHM) - g2(iDM)



   !!  compute intermediate values of y,x and dydx 
     do i=1,n 
        y(i) = ysav(i) + a31*g1(i) + a32*g2(i) 
     enddo
     x = xsav + a3x*h 
     call derivs(x,y,dydx) 


   !!  set up and solve the right hand side for g3 
     do i=1,n 
        g3(i)  = dydx(i) + (c31*g1(i) + c32*g2(i))/h 
     enddo


     do j = 1, n
        do i = 1, n
           av(i,j) = dmat(i,j)
        enddo
     enddo
     do i = 1, n
!        av(i,n+1) = g3(i)
        bb(i) = g3(i)
     enddo
!     call pchem_gift(av,nmax,nmaxp1)
     call leqs(av,bb,n-1,nmax)
!     call ludcmp(av,n-1,nmax,indx,d)
!     call lubksb(av,n-1,nmax,indx,bb) 
    do i = 1, n
        g3(i) = bb(i)  !av(i,n+1) 
     enddo
        g3(iELEC) = g3(iHP) + g3(iDP) + 2.0*g3(iHEPP) + g3(iHEP) + g3(iHDP) + g3(iH2P) - g3(iHM) - g3(iDM)



   !!  set up and solve the right hand side for g4 
     do i=1,n 
        g4(i)  = dydx(i) + (c41*g1(i) + c42*g2(i) + c43*g3(i))/h 
     end do


     do j = 1, n
        do i = 1, n
           av(i,j) = dmat(i,j)
        enddo
     enddo
     do i = 1, n
!        av(i,n+1) = g4(i)
        bb(i) = g4(i)
     enddo
!     call pchem_gift(av,nmax,nmaxp1)
     call leqs(av,bb,n-1,nmax)
!     call ludcmp(av,n-1,nmax,indx,d)
!     call lubksb(av,n-1,nmax,indx,bb)
     do i = 1, n
        g4(i) = bb(i)  !av(i,n+1) 
     enddo
        g4(iELEC) = g4(iHP) + g4(iDP) + 2.0*g4(iHEPP) + g4(iHEP) + g4(iHDP) + g4(iH2P) - g4(iHM) - g4(iDM)



   !!  compute the third and fourth order estimates of y 
     do i=1,n 
        y(i)   = ysav(i) + b1*g1(i) + b2*g2(i) + b3*g3(i) + b4*g4(i) 
        err(i) = e1*g1(i) + e2*g2(i) + e3*g3(i) + e4*g4(i) 
     enddo
       y(iELEC) = y(iHP) + y(iDP) + 2.0*y(iHEPP) + y(iHEP) + y(iHDP) + y(iH2P) - y(iHM) - y(iDM)

     x = xsav + h 

     if (x .eq. xsav) then 
        if (pchem_meshMe .EQ. MASTER_PE) print *, 'step size not significant in pchem_rosenGift' 
        call Logfile_stampMessage('[pchem_rosenGift] step size not significant')
        call Driver_abortFlash('ERROR: step size not significant in pchem_rosenGift')
     end if


   !!  determine the scaled accuracy 
     errmax = 0.0e0 
     do i=1,n 
        errmax = max(errmax,abs(err(i)/yscal(i))) 
     enddo
     errmax = errmax/eps 

   !!  if the step succeded, compute the size of the next step and return 
     if (errmax .le. 1.0) then 
        hdid = h 
        if (errmax .gt. errcon) then 
           hnext = safety * h * errmax**pgrow 
        else 
           hnext = grow * h 
        end if
        jcounts = jcounts + jtry
!       print *, 'JCOUNTS: ', jcounts
!       print *, 'HNEXT: ', hnext
!       if(jcounts .gt. 10000) call Driver_abortFlash('JCOUNTS IS HUGE. YOU HAVE A PROBLEM!!')
        return 

      !!   if the step did not succeed cut the stepsize and try again 
     else 
        hnext = safety * h * errmax**pshrnk 
        h     = sign(max(abs(hnext),shrnk*abs(h)),h) 
     end if
  enddo

!     print *, 'END OF ROSEN'

!!  too many tries
  if (pchem_meshMe .EQ. MASTER_PE) then
      print *, 'ERROR exceeded maxtry in routine pchem_rosenGift' 
  endif
  
  do i=1,14
    print *, 'SPECIES AFTER(',i,')= ', y(i)
    print *, 'SPECIES BEFORE(',i,')= ', yysav(i)
  enddo


  call Logfile_stampMessage('[pchem_rosenGift] ERROR exceeded maxtry')
  call Driver_abortFlash('ERROR: exceeded maxtry in pchem_rosenGift')

end subroutine  pchem_rosenGift
!!---------------------------------------------------------------------------------

