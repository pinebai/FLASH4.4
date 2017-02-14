!..this file contains 4 routines to drive the integration of
!..the NEI equations
!
!..routine neiint drives the integration of the odes
!
!..routine stifbs4 drives a bader-deuflhard step with ma28 algebra
!..routine simpr3 is a bader-deuflhard stepper with ma28 algebra
!..routine pzextr does extrapolations for the simpr3 routine
!
!
!
!
!
!
subroutine neiint(start,stptry,stpmin,stopp,bc,                   &
     &                  eps,dxsav,kmax,                                 &
     &                  xrk,yrk,xphys,yphys,xlogi,ylogi,                &
     &                  nok,nbad,kount,odescal,iprint,                  &
     &                  derivs,jakob,bjakob,steper)
  !
  ! *****************************************************************
  !
  !  LAST AUTHOR UPDATE: CHICAGO AUGUST 24, 2001
  !  LAST UPDATE:       CHICAGO OCTOBER 15, 2010
  !
  !  NAME:    NEIINT
  !
  !  PURPOSE: ode integrator for stiff odes with a sparse jacobian
  !
  !  INPUTS:
  !        START    REAL      beginning integration point
  !        STPTRY   REAL      suggested first step size
  !        STPMIN   REAL      minimum allowable step size
  !        STOPP    REAL      ending integration point
  !        BC       REAL      initial conditions, 
  !                           array of physical dimension yphys
  !        EPS      REAL      desired fraction error during the integration
  !        DXSAV    REAL      incremental value of indep variable 
  !                           at which to store the sol.
  !                           if zero     -> sol. is stored at every 
  !                                          intermediate point
  !                           if not zero -> sol. is done and saved 
  !                                          at every dxsav point
  !        KMAX     INTEGER   maximum number of sol. points to store, 
  !                           kkmax < xphys
  !        ODESCAL  REAL      error scaling factor
  !        IPRINT   INTEGER   determines if the sol. is printed as it evolves
  !        DERIVS   EXTERNAL  name of the routine that contains the odes
  !        JAKOB    EXTERNAL  rout. that contains the jacobian of the odes
  !        BJAKOB   EXTERNAL  rout. that sets the pointers of the sparse jacobian
  !        STEPER   EXTERNAL  rout. that will take a single step
  !
  !  NOTE: At the end of the computation, BC is updated
  !
  !  OUTPUTS:
  !        XRK      REAL      the independent variable solution
  !                           array of physical dimension XPHYS, 
  !                           logical dimension XLOGI
  !        YRK      REAL      the dependent variable solution
  !                           array of physical dimension (YPHYS,XPHYS),
  !                           logical  dimension (YLOGI,XLOGI)
  !        NOK      INTEGER   number of succesful steps taken
  !        NBAD     INTEGER   number of bad steps taken, 
  !                           bad but retried and then succesful
  !        KOUNT    INTEGER   total number of steps stored in arrays xrk and yrk
  !
  !  SUBROUTINES:
  !
  !  REFERENCES:
  !           none
  !
  !  MODIFICATION HISTORY:
  !           ORIGINAL ROUTINE: NETINT of the burning module
  !           August 2001: modified for the NEI module - S. ORLANDO
  !           October 2010: added intent attributes -  N. Flocke, K. Weide
  !
  ! *****************************************************************
  !
  !
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  
  !
  !..declare  
  external         derivs,jakob,bjakob,steper

  integer, intent(in)  :: kmax
  integer, intent(in)  :: iprint
  integer, intent(in)  :: xphys
  integer, intent(in)  :: yphys
  integer, intent(in)  :: xlogi
  integer, intent(in)  :: ylogi
  integer, intent(out) :: nok
  integer, intent(out) :: nbad
  integer, intent(out) :: kount

  real, intent(in)  :: start
  real, intent(in)  :: stptry
  real, intent(in)  :: stpmin
  real, intent(in)  :: stopp
  real, intent(inout) :: bc(yphys)
  real, intent(in)  :: eps
  real, intent(in)  :: dxsav
  real, intent(in)  :: odescal
  real, intent(out) :: xrk(xphys)
  real, intent(out) :: yrk(yphys,xphys)

  integer, parameter :: nmax=30, stpmax=10000
  real, parameter :: zero=0.0e0, one=1.0e0, tiny=1.0e-15

  integer, save :: i,j,nstp
  real, save :: yscal(nmax),y(nmax),dydx(nmax),x,xsav,h,hdid,hnext

  !
  !..here are the format statements for printouts as we integrate
100 format(1x,i4,e10.2)
101 format(1x,12e10.2)
102 format(1x,5(a,' ',e10.2))
103 format(1x,5(a,' ',i6))
  !
  !..initialize   
  if (ylogi .gt. yphys) then
     write(*,*) 'ylogi > yphys in routine neiint'
     call Driver_abortFlash('ylogi > yphys in routine neiint')
     stop
  end if
  !
  if (yphys .gt. nmax)  then
     write(*,*) 'yphys > nmax in routine neiint'
     call Driver_abortFlash('yphys > nmax in routine neiint')
     stop
  end if
  !
  x     = start   
  h     = sign(stptry,stopp-start)
  nok   = 0 
  nbad  = 0
  kount = 0   
  !
  !..store the first step 
  do i = 1,ylogi
     y(i) = bc(i)  
  enddo
  xsav = x-2.0e0*dxsav
  !
  !..take at most stpmax steps
  do nstp = 1,stpmax
     !
     !..positive definite population fractions
     do i = 1,ylogi
        y(i) = max(y(i),1.0e-30)
     enddo
     call derivs(x,y,dydx)
     !
     !..scaling vector used to monitor accuracy
     do i = 1,ylogi
        yscal(i) = max(odescal,abs(y(i)))
     enddo
     !
     !..store intermediate results
     if (kmax .gt. 0) then
        if ( (abs(dxsav) - abs(x-xsav)) .le. tiny) then
           if ( kount .lt. (kmax-1) ) then
              kount         = kount+1
              xrk(kount)    = x
              do i = 1,ylogi
                 yrk(i,kount) = y(i)
              enddo
              if (iprint .eq. 1) then
                 write(6,100) kount,xrk(kount)
                 write(6,101) (yrk(j,kount), j=1,ylogi)
              end if
              xsav = x
           end if
        end if
     end if
     !
     !..if the step can overshoot the stop point or the dxsav increment then cut it
     if ((x+h-stopp)*(x+h-start).gt.zero) h = stopp - x
     if (dxsav.ne.zero .and. h.gt.(xsav-x+dxsav)) h = xsav + dxsav-x
     !
     !..do an integration step
     call steper(y,dydx,ylogi,x,h,eps,yscal,hdid,hnext,              &
          &              derivs,jakob,bjakob)
     if (hdid.eq.h) then
        nok = nok+1
     else
        nbad = nbad+1
     end if
     !
     !..this is the normal exit point, save the final step
     if (nstp.eq.stpmax.or.(x-stopp)*(stopp-start).ge.zero) then
        do i = 1,ylogi
           bc(i) = y(i)
        enddo
        !
        if (kmax.ne.0) then
           kount         = kount+1
           xrk(kount)    = x
           do i = 1,ylogi
              yrk(i,kount) = y(i)
           end do
           if (iprint .eq. 1) then
              write(6,100) kount,xrk(kount)
              write(6,101) (yrk(j,kount), j=1,ylogi)
           end if
        end if
        return
     end if
     !
     !
     !..set the step size for the next iteration; stay above stpmin
     h = hnext
     if (abs(hnext).lt.stpmin) then
        write(6,*) ' '
        write(6,102) 'hnext=',hnext,' stpmin=',stpmin
        write(6,*) 'hnext < stpmin in netint'
        write(6,*) ' '
        call Driver_abortFlash('hnext < stpmin in netint')
     end if
     !
     !..back for another iteration or death
  enddo
  !
  write(6,*)
  write(6,103)'stpmax=',stpmax
  write(6,*)  'more than stpmax steps required in netint'
  write(6,*)
  call Driver_abortFlash('more than stpmax steps required in netint')
  !
  !
end subroutine neiint
!
!
!
!*****************************************************************
!*****************************************************************
!
!
!
subroutine stifbs4(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,         &
     &                   derivs,jakob,bjakob)
  !
  use Driver_interface, ONLY : Driver_abortFlash
  use Ionize_data, ONLY : ifirst => ion_ifirst
  implicit none

  !.. 
  !..for sparse analytic jacobians, ma28 linear algebra 
  !..
  !..semi-implicit extrapolation step for stiff ODEs with monitoring 
  !..of local truncation error to adjust stepsize. Inputs are the dependent  
  !..variable vector y(nv) and its derivative dydx(nv) at the starting of the 
  !..independent variable x. Also inputs are the stepsize to be attempted htry, 
  !..the required accuracy eps, and the vector yscal against which the error is 
  !..scaled. on output, y and x are replaced by their new values, hdid is the  
  !..stepsize actually accomplished, and hnext is the estimated next stepsize. 
  !..dervs is a user supplied function that computes the right hand side of 
  !..the equations.
  !..
  !..1/scalmx is the maximum increase in the step size allowed
  !
  !
  !..declare  

  integer, intent(in) :: nv
  real, intent(in) :: dydx(nv), htry, eps, yscal(nv)
  real, intent(inout) :: y(nv), x
  real, intent(out) :: hdid, hnext
  external         derivs,jakob,bjakob
  
  logical, save :: first,reduct
  integer, parameter :: nmax = 30, kmaxx=7, imax=kmaxx+1
  integer, parameter :: naij=200, n5 = 5*nmax, n8=8*nmax
  integer, save :: iloc(naij),jloc(naij),ivect(naij),jvect(naij),&
                   ikeep(n5),iw(n8),flag,nzo,ii
  integer, save :: i,iq,k,kk,km,kmax,kopt,nvold,nseq(imax)
  real, save :: eps1,epsold,errmax,fact,h,red,scale,xwork,wrkmin,&
                xest,xnew,a(imax),alf(kmaxx,kmaxx),err(kmaxx),   &
                yerr(nmax),ysav(nmax),yseq(nmax),&
                dfdy(naij),amat(naij),w(nmax),u
  real, parameter :: safe1 = 0.25e0, safe2 = 0.7e0, redmax=1.0e-5,&
                     redmin = 0.7e0, tiny = 1.0e-30, scalmx = 0.1e0
  !
  !! integer          ifirst
  !! common /neibjac/ ifirst  !Commented out by Chris 3 July 08, as ion_ifirst 
  !! is set in subroutine neimn.  We never refer to neibjac common block.

  !
  !
  data             first/.true./, epsold/-1.0e0/, nvold/-1/
  data             nseq /2, 6, 10, 14, 22, 34, 50, 70/
  !      data             ifirst/0/ 
  !
  !
  !..get and copy the nonzero locations
  if (ifirst .eq. 0) then
     !
     ifirst = 1
     !
     call bjakob(iloc,jloc,nzo,naij)
     !
     do i=1,nzo
        ivect(i) = iloc(i)
        jvect(i) = jloc(i)
     enddo
     !
     !..force the diagonal to be the pivot elements
     do i=1,nzo
        amat(i) = 1.0e-10
        if (ivect(i) .eq. jvect(i)) amat(i) = 1.0e0
     enddo
     !
     u  = 0.1e0
     call ma28ad(nv,nzo,amat,naij,iloc,naij,jloc,u,ikeep,iw,w,flag)
     !
     if (flag .lt. 0) then
        write(6,*) 'error in ma28ad flag',flag
        write(*,*) 'more than stpmax steps required in netint'  
        call Driver_abortFlash('more than stpmax steps required in netint')
     end if
     !
  end if
  !
  !..a new tolerance or a new number , so reinitialize
  if (eps .ne. epsold  .or.  nv .ne. nvold) then
     hnext = -1.0e29
     xnew  = -1.0e29
     eps1  = safe1 * eps
     !
     !..compute the work coefficients a_k
     a(1)  = nseq(1) + 1
     do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
     enddo
     !
     !..compute alf(k,q)
     do iq=2,kmaxx
        do k=1,iq-1
           alf(k,iq) = eps1**((a(k+1) - a(iq+1)) /                        &
                &               ((a(iq+1) - a(1) + 1.0e0) * (2*k + 1)))
        enddo
     enddo
     epsold = eps
     nvold  = nv
     !
     !..add cost of jacobians to work coefficients
     a(1)   = nv + a(1)
     do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
     enddo
     !
     !..determine optimal row number for convergence
     do kopt=2,kmaxx-1
        if (a(kopt+1) .gt. a(kopt)*alf(kopt-1,kopt)) go to 01
     enddo
01   kmax = kopt
  end if
  !
  !
  !..save the starting values 
  h    = htry
  do i=1,nv  
     ysav(i)  = y(i)
  enddo
  !
  !..get the sparse jacobian in dfdy
  call jakob(x,y,dfdy,nzo)
  !
  !..a new stepsize or a new integration, re-establish the order window
  if (h .ne. hnext  .or.  x .ne. xnew) then
     first = .true.
     kopt = kmax
  end if
  reduct = .false.
  !
  !..evaluate the sequence of semi implicit midpoint rules
02 do k=1,kmax
     xnew = x + h
     !
     if (xnew .eq. x) then
        write(*,*) 'stepsize too small in routine stiffbs'
        write(6,110) xnew,x,k
110     format(1x,2e11.3,i6)
        !
        !        do ii=1,nv
        !         write(6,111) y(ii),y(ii)*aion(ii),aion(ii),ii
        ! 111     format(1x,1p3e11.3,i6)
        !        enddo
        !
        call Driver_abortFlash('stepsize too small in routine stiffbs')
     end if
     !
     !
     call simpr3(ysav,dydx,dfdy,nmax,nv,x,h,nseq(k),yseq,             &
          &             nzo,amat,naij,ivect,jvect,jloc,ikeep,iw,w,flag,      &
          &             derivs) 
     !
     xest = (h/nseq(k))**2 
     call pzextr(k,xest,yseq,y,yerr,nv) 
     !
     !
     !..compute normalized error estimate
     if (k .ne. 1) then
        errmax = tiny
        do i=1,nv
           errmax = max(errmax,abs(yerr(i)/yscal(i)))
        enddo
        errmax   = errmax/eps   
        km = k - 1
        err(km) = (errmax/safe1)**(1.0e0/(2*km+1))
     end if
     !
     !..in order window
     if (k .ne. 1  .and. (k .ge. kopt-1  .or. first)) then
        !
        !..converged
        if (errmax .lt. 1.0) go to 04
        !
        !
        !..possible step size reductions
        if (k .eq. kmax  .or.  k .eq. kopt + 1) then
           red = safe2/err(km)
           go to 3
        else if (k .eq. kopt) then
           if (alf(kopt-1,kopt) .lt. err(km)) then
              red = 1.0e0/err(km)
              go to 03
           end if
        else if (kopt .eq. kmax) then
           if (alf(km,kmax-1) .lt. err(km)) then
              red = alf(km,kmax-1) * safe2/err(km)
              go to 03
           end if
        else if (alf(km,kopt) .lt. err(km)) then
           red = alf(km,kopt-1)/err(km)
           go to 03
        end if
     end if
     !
  enddo
  !
  !
  !..reduce stepsize by atleast redmin and at mosr redmax
03 red    = min(red,redmin)
  red    = max(red,redmax)
  h      = h * red
  reduct = .true.
  go to 2
  !
  !..successful step; get optimal row for convergence and corresponding stepsize
04 x = xnew
  hdid = h
  first = .false.
  wrkmin = 1.0e35
  do kk=1,km
     fact = max(err(kk),scalmx)
     xwork = fact * a(kk+1)
     if (xwork .lt. wrkmin) then
        scale  = fact
        wrkmin = xwork
        kopt   = kk + 1
     end if
  enddo
  !
  !..check for possible order increase, but not if stepsize was just reduced
  hnext = h/scale
  if (kopt .ge. k  .and.  kopt .ne. kmax  .and.  .not.reduct) then
     fact = max(scale/alf(kopt-1,kopt),scalmx)
     if (a(kopt+1)*fact .le. wrkmin) then
        hnext = h/fact
        kopt = kopt + 1
     end if
  end if
  return
end subroutine stifbs4
!
!
!*****************************************************************
!*****************************************************************
!
!
!
subroutine simpr3(y,dydx,dfdy,nmax,n,xs,htot,nnstep,yout,         &
     &                  nzo,a,naij,ivect,jvect,jloc,ikeep,iw,w,flag,    &
     &                  derivs) 
  !
  use Driver_interface, ONLY : Driver_abortFlash
  
  !
  !..
  !..an implicit midpoint stepper, for ma28 sparse linear algebra. 
  !.. 
  !..declare 
  implicit none
  integer, intent(in) :: n,nmax,nnstep,nzo,naij,ivect(naij),jvect(naij),&
                         jloc(naij),ikeep(1),flag
  integer, intent(inout) :: iw(1)
  real, intent(in) :: y(n),dydx(n),dfdy(naij),xs,htot
  real, intent(inout) :: a(naij),w(1)
  real, intent(out) :: yout(n)
  external         derivs
  
  integer, parameter :: nmaxx=30
  integer            :: i,nn 
  real, save :: h,x,del(nmaxx),ytemp(nmaxx)

  !
  !
  !..stepsize this trip, and make the a matrix 
  h = htot/nnstep 
  do i=1,nzo 
     a(i) = -h * dfdy(i) 
     if (ivect(i) .eq. jvect(i)) a(i) = 1.0e0 + a(i) 
  enddo
  !
  !..numeric decomp 
  call ma28bd(n,nzo,a,naij,ivect,jvect,jloc,ikeep,iw,w,flag) 
  !
  !
  if (flag .lt. 0) then 
     write(6,*) 'error in ma28bd flag',flag 
     call Driver_abortFlash('error in ma28bd flag')
  end if
  !
  !..use yout as temporary storage; the first step 
  do i=1,n 
     yout(i) = h * dydx(i) 
  enddo
  !
  call ma28cd(n,a,naij,jloc,ikeep,yout,w,1) 
  !
  do i=1,n 
     del(i)   = yout(i) 
     ytemp(i) = y(i) + del(i) 
  enddo
  x = xs + h 
  !
  call derivs(x,ytemp,yout) 
  !
  !
  !..use yout as temporary storage; general step 
  do nn=2,nnstep 
     do i=1,n 
        yout(i) = h*yout(i) - del(i) 
     enddo
     !
     call ma28cd(n,a,naij,jloc,ikeep,yout,w,1) 
     !
     do i=1,n 
        del(i)   = del(i) + 2.0e0 * yout(i) 
        ytemp(i) = ytemp(i) + del(i) 
     enddo
     x = x + h 
     !
     call derivs(x,ytemp,yout) 
     !
  enddo
  !
  !..take the last step 
  do i=1,n 
     yout(i) = h * yout(i) - del(i)  
  enddo
  !
  call ma28cd(n,a,naij,jloc,ikeep,yout,w,1) 
  !
  do i=1,n 
     yout(i) = ytemp(i) + yout(i) 
  enddo
  return 
end subroutine simpr3
!
!
!*****************************************************************
!*****************************************************************
!
!
!
subroutine pzextr(iest,xest,yest,yz,dy,nv) 

  implicit none
  !
  !..use polynomial extrapolation to evaluate nv functions at x=0 by fitting 
  !..a polynomial to a sequence of estimates with progressively smaller values 
  !..x=xest, and corresponding function vectors yest(1:nv). the call is number  
  !..iest in the sequence of calls. extrapolated function values are output as  
  !..yz(1:nv), and their estimated error is output as dy(1:nv) 
  !.. 
  !..declare 
  integer, intent(in) :: iest,nv
  real, intent(in) :: xest,yest(nv)
  real, intent(out) :: yz(nv),dy(nv)
  
  integer, parameter :: nmax=50, imax=13
  
  integer, save :: j,k1
  real, save :: delta,f1,f2,q,d(nmax),qcol(nmax,imax),x(imax) 
  !
  !
  !..save current independent variables 
  x(iest) = xest 
  do j=1,nv 
     dy(j) = yest(j) 
     yz(j) = yest(j) 
  enddo
  !
  !..store first estimate in first column 
  if (iest .eq. 1) then 
     do j=1,nv 
        qcol(j,1) = yest(j) 
     enddo
     !
  else 
     do j=1,nv 
        d(j) = yest(j) 
     enddo
     do k1=1,iest-1 
        delta = 1.0e0/(x(iest-k1) - xest) 
        f1    = xest * delta 
        f2    = x(iest-k1) * delta 
        !
        !..propagate tableu 1 diagonal more 
        do j=1,nv 
           q          = qcol(j,k1) 
           qcol(j,k1) = dy(j) 
           delta      = d(j) - q 
           dy(j)      = f1*delta 
           d(j)       = f2*delta 
           yz(j)      = yz(j) + dy(j) 
        enddo
     enddo
     do j=1,nv 
        qcol(j,iest) = dy(j) 
     enddo
  end if

  return 

end subroutine pzextr
