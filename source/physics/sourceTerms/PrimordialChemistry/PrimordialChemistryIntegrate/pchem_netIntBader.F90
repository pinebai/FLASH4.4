!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryIntegrate/pchem_netIntBader
!!
!! NAME
!!   pchem_baderMa28  
!!
!! SYNOPSIS
!!  subroutine pchem_baderMa28(real(IN)       ::y(:),
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
!!  for sparse analytic jacobians, ma28 linear pchem_algebra 
!!  
!!  semi-implicit extrapolation step for stiff ode's with monitoring 
!!  of local truncation error to adjust stepsize. inputs are the dependent  
!!  variable vector y(nv) and its derivative dydx(nv) at the starting of the 
!!  independent variable x. also input are the stepsize to be attempted htry, 
!!  the required accuracy eps, and the vector yscal against which the error is 
!!  scaled. on output, y and x are replaced by their new values, hdid is the  
!!  stepsize actually accomplished, and hnext is the estimated next stepsize. 
!!  derivs is a user supplied function that computes the right hand side of 
!!  the equations.
!!  
!!  1/scalmx is the maximum increase in the step size allowed
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
!!  this file contains 5 routines to work with bader-duflhard integration steps
!!  for the integration nuclear reaction networks with 2 choices for the linear pchem_algebra.
!!
!!  routine pchem_baderMa28 drives a bader-deuflhard step with ma28 pchem_algebra
!!  routine pchem_baderStepMa28 is a bader-deuflhard stepper with ma28 pchem_algebra
!!  routine pchem_baderGift drives a bader-deuflhard step with gift pchem_algebra
!!  routine pchem_baderStepGift is a bader-deuflhard stepper with ma28 pchem_algebra
!!  routine pchem_pzExtr does extrapolations for any of the *Step* routines
!!
!!***
!!---------------------------------------------------------------------------------




subroutine pchem_baderGift(y,dydx,nv,x,htry,eps,yscal,hdid,hnext, & 
     &                       derivs,jakob,bjakob)

  use PrimordialChemistry_data
  use pchem_interface
  use Driver_interface, ONLY : Driver_abortFlash
  !  Ummm.... bit of a mystery why I can use the interfaces in pchem_netIntRosen but not here.
  use pchem_integrateInterface, ONLY: pchem_baderStepGift, pchem_pzExtr
  use pchem_data
 
  implicit none

!!  for dense analytic jacobians, gift pchem_algebra
!!  
!!  semi-implicit extrapolation step for stiff ode's with monitoring 
!!  of local truncation error to adjust stepsize. inputs are the dependent  
!!  variable vector y(nv) and its derivative dydx(nv) at the starting of the 
!!  independent variable x. also input are the stepsize to be attempted htry, 
!!  the required accuracy eps, and the vector yscal against which the error is 
!!  scaled. on output, y and x are replaced by their new values, hdid is the  
!!  stepsize actually accomplished, and hnext is the estimated next stepsize. 
!!  dervs is a user supplied function that computes the right hand side of 
!!  the equations.
!!  
!!  1/scalmx is the maximum increase in the step size allowed

  interface   ! = pchem_network
     subroutine derivs(tt,y,dydt)   !! == pchem_network
       implicit none
       real, intent(IN) :: tt
       real, intent(INOUT), dimension(*)  :: y
       real, intent(OUT), dimension(*) :: dydt
     end subroutine derivs

     subroutine bjakob(iloc,jloc,nzo,np)
       implicit none
       integer, intent(IN)  ::   iloc(*),jloc(*),np
       integer, intent(OUT) ::   nzo
     end subroutine bjakob
  end interface

!!  declare arguments
  external               jakob
  integer, intent(IN) :: nv
  real, intent(IN)    :: dydx(nv), yscal(nv), htry, eps
  real, intent(INOUT) :: x, y(nv)
  real, intent(INOUT)   :: hdid, hnext


!!  declare local variables; you  need to save all of them or all hell breaks loose.
  logical, save       ::   first,reduct
  integer, parameter  ::   nmax = 30, kmaxx=7, imax=kmaxx+1   
  integer, save       ::   i,iq,k,kk,km,kmax,kopt
  integer, save       :: nvold,nseq(imax),ii
  real, save          :: eps1,epsold,errmax,fact,h,red,scale,xwork,wrkmin, & 
       &                 xest,xnew,a(imax),alf(kmaxx,kmaxx),err(kmaxx), & 
       &                 yerr(nmax),ysav(nmax),yseq(nmax),& 
       &                 dfdy(nmax,nmax)
  real, parameter    :: safe1 = 0.25e0, safe2 = 0.7e0, redmax=1.0e-5, & 
       &                  redmin = 0.7e0, tiny = 1.0e-30, scalmx = 0.1e0


  data             first/.true./, epsold/-1.0e0/, nvold/-1/
  data             nseq /2, 6, 10, 14, 22, 34, 50, 70/

  !couple of counters
  integer             :: l,j,z

!   print *,'BaderGift'
!   print *,'y: ', y
!!  a new tolerance or a new number , so reinitialize
  if (eps .ne. epsold  .or.  nv .ne. nvold) then
     hnext = -1.0e29
     xnew  = -1.0e29
     eps1  = safe1 * eps

   !!  compute the work coefficients a_k
     a(1)  = nseq(1) + 1
     do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
     enddo

   !!  compute alf(k,q)
     do iq=2,kmaxx
        do k=1,iq-1
           alf(k,iq) = eps1**((a(k+1) - a(iq+1)) / & 
                &               ((a(iq+1) - a(1) + 1.0e0) * (2*k + 1)))
        enddo
     enddo
     epsold = eps
     nvold  = nv

   !!  add cost of jacobians to work coefficients
     a(1)   = nv + a(1)
     do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
     enddo

   !!  determine optimal row number for convergence
     do kopt=2,kmaxx-1
        if (a(kopt+1) .gt. a(kopt)*alf(kopt-1,kopt)) go to 01
     enddo
01   kmax = kopt
  end if


!!  save the starting values 
  h    = htry
  do i=1,nv  
     ysav(i)  = y(i)
  enddo

!!  get the dense jacobian in dfdy
  !      call Timers_start ("jacobian (burn)")

 ! call jakob(x,y,dfdy,nv,nmax)
  call jakob(x,y,dfdy,nv,nmax) 
 
!   print *,'BaderGift2'
!   print *,'y: ', y

  !Making sure the matricies match
  !THEY DO!  
!  print *, 'INITALIZING DFDY'
!  do l=1,nv
!     do j=1,nv
!        print *, 'jakob(',l,j,')=',dfdy(l,j)
!     enddo
!  enddo

 !  print *, dfdy
 !  print *, x
 !  print *, y
 !  print *, nv
 !  print *, nmax
  

  !      call Timers_stop ("jacobian (burn)")

!!  a new stepsize or a new integration, re-establish the order window
  if (h .ne. hnext  .or.  x .ne. xnew) then
     first = .true.
     kopt = kmax
  end if
  reduct = .false.

!!  evaluate the sequence of semi implicit midpoint rules
02 do k=1,kmax
     xnew = x + h

     if (xnew .eq. x) then
        write(*,*) 'stepsize too small in routine baderGift'
        write(*,110) xnew,x,k
110     format(1x,1p2e11.3,i6)

        do ii=1,nv
           write(*,111) y(ii),y(ii)*aion(ii),aion(ii),ii
111        format(1x,1p3e11.3,i6)
        enddo

        call Driver_abortFlash('ERROR in baderGift: stepsize too small')
     end if



   !!  the semi implicit midpoint rules for this sequence
     call pchem_baderStepGift(ysav,dydx,dfdy,nmax,nv,x,h,nseq(k),yseq,derivs) 

   !!  extrapolate the error to zero
     xest = (h/nseq(k))**2 
!     print *, 'k: ', k
!     print *, 'xest: ', xest
!     print *, 'yseq: ', yseq
!     print *, 'y: '   , y
!     print *, 'yerr: ', yerr
!     print *, 'nv:   ', nv
!     print *, 'dfdy ', dfdy

!   print *,'BaderGift3'
!   print *,'y: ', y
!   print *,'yseq: ', yseq

     call pchem_pzExtr(k,xest,yseq,y,yerr,nv)
     

   !!  compute normalized error estimate
     if (k .ne. 1) then
        errmax = tiny
        do i=1,nv
           if( i .le. 3) then
      !     print *,k,'errmax:',errmax,'xest: ',xest
      !     print *,i,yerr(i),yscal(i)
          endif
           errmax = max(errmax,abs(yerr(i)/yscal(i)))
        enddo
        errmax   = errmax/eps   
        km = k - 1
        err(km) = (errmax/safe1)**(1.0e0/(2*km+1))
     end if

   !!  in order window
     if (k .ne. 1  .and. (k .ge. kopt-1  .or. first)) then

      !!   converged
        if (errmax .lt. 1.0) go to 04


      !!   possible step size reductions
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

  enddo


!!  reduce stepsize by atleast redmin and at mosr redmax
03 red    = min(red,redmin)
  red    = max(red,redmax) 
  h      = h * red
  reduct = .true.
  go to 2

!!  successful step; get optimal row for convergence and corresponding stepsize
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

!!  check for possible order increase, but not if stepsize was just reduced
  hnext = h/scale
  if (kopt .ge. k  .and.  kopt .ne. kmax  .and.  .not.reduct) then
     fact = max(scale/alf(kopt-1,kopt),scalmx)
     if (a(kopt+1)*fact .le. wrkmin) then
        hnext = h/fact
        kopt = kopt + 1
     end if
  end if
  return
end subroutine pchem_baderGift

!!---------------------------------------------------------------------------------


subroutine pchem_baderStepGift(y,dydx,dfdy,nmax,n,xs,htot,nnstep,yout,  & 
     &                      derivs) 

  use pchem_interface
  use PrimordialChemistry_data
  implicit none
  
!!   
!!  an implicit midpoint stepper, for GIFT sparse linear pchem_algebra. 
!!   
!!  declare arguments
  external               derivs 
  integer, intent(IN) :: n, nmax
  integer, intent(IN) :: nnstep
  real, intent(IN)    :: xs, htot, y(n), dydx(n), dfdy(nmax,nmax)
  real, intent(INOUT) :: yout(n)

!!  declare local variables
  integer, parameter  :: nmaxx=30 
  integer, save       ::  i,j,nn
  real, save          ::  h,x,del(nmaxx),ytemp(nmaxx)

!!  for the gift linear pchem_algebra
  integer, parameter  :: nmaxp1 = nmaxx + 1
  real, save          :: dmat(nmaxx,nmaxx)!, av(nmaxx,nmaxp1)
  real, save          :: av(nmaxx,nmaxx)
!  real, save         :: av(n,n)
  real, save          :: bb(nmaxx) !This is the solution
  integer             :: indx(nmaxx)
  real                :: d


!!  stepsize this trip, and make the a matrix 
  h = htot/nnstep 

  !      call Timers_start ("pchem_algebra (burn)")
  ! print *, 'this is n inside netIntBader: ', n
  do j=1,n
     do i=1,n
        dmat(i,j) = -h * dfdy(i,j) 
     enddo
  enddo
  do i=1,n
     dmat(i,i) = 1.0e0 + dmat(i,i) 
  end do



!!  use yout as temporary storage; the first step 
  do i=1,n 
     yout(i) = h * dydx(i) 
  enddo



  !Adding something new here, going to call one of my new functions to
  !turn the jacobian into a LU decomposed matrix


  do j = 1, n
     do i = 1, n
        av(i,j) = dmat(i,j)
!       print *,i,j,av(i,j)
     enddo
  enddo
  
!  print *, 'SETUP TEST MATRIX'
!  do i = 1,n
!     do j = 1,n
!        if (i .eq. j) then 
!          av(i,i) = i
!        else
!          av(i,j) = 0.
!        endif
!        av(i,j) = (i+j)**(1.5)
!        print *, i,j,av(i,j)
!     enddo
!  enddo


!  call ludcmp(av,n-1,nmaxx,indx,d) !!output is a lu decomposition of av

!  print *, 'AFTER LUDCMP'
!  do i = 1,n
!     do j = 1,n
!        print *, i,j,av(i,j)
!     enddo
!  enddo

!  print *, 'called ludcmp!'
  do i = 1,n 
    ! av(i,n+1) = yout(i)
      bb(i) = yout(i)  !NEW, this is the solution
!      print *, 'bb(',i,')=', bb(i)
  enddo


!  call pchem_gift(av,nmaxx,nmaxp1)
!  call lubksb(av,n-1,nmaxx,indx,bb)
  call leqs(av,bb,n-1,nmaxx) 
!!TRYING GUASSIAN ELIMINATION
!   print *, 'calling leqs'
!   call leqs(av,bb,n-1,nmaxx)
!   print *, 'called leqs'


  do i = 1, n
    ! yout(i) = av(i,n+1) 
    yout(i) = bb(i)
!    print *, i,bb(i)
  enddo
    bb(iELEC) = bb(iHP) + bb(iDP) + bb(iHEP) + 2.0*bb(iHEPP) + bb(iH2P) + bb(iHDP) - bb(iHM) - bb(iDM)
  !      call Timers_stop ("pchem_algebra (burn)")

  do i=1,n 
     del(i)   = yout(i) 
     ytemp(i) = y(i) + del(i) 
  enddo
  x = xs + h 

  !      call Timers_start ("derivs (burn)")

  call derivs(x,ytemp,yout) 

  !      call Timers_stop ("derivs (burn)")


!!  use yout as temporary storage; general step 
  do nn=2,nnstep 
     do i=1,n 
        yout(i) = h*yout(i) - del(i) 
     enddo

     !       call Timers_start ("pchem_algebra (burn)")

!     do j = 1, n
!        do i = 1, n
!           av(i,j) = dmat(i,j)
!        enddo
!     enddo
     do i = 1,n 
        !av(i,n+1) = yout(i)
        bb(i) = yout(i)
     enddo


    ! call pchem_gift(av,nmaxx,nmaxp1)
    ! call ludcmp(av,n,nmaxx,indx,d) 
 !    call lubksb(av,n-1,nmaxx,indx,bb)
     call leqs(av,bb,n-1,nmaxx)
     do i = 1, n
      !  yout(i) = av(i,n+1) 
      yout(i) = bb(i)
     enddo
     bb(iELEC) = bb(iHP) + bb(iDP) + bb(iHEP) + 2.0*bb(iHEPP) + bb(iH2P) + bb(iHDP) - bb(iHM) - bb(iDM)
     !       call Timers_stop ("pchem_algebra (burn)")


     do i=1,n 
        del(i)   = del(i) + 2.0e0 * yout(i) 
        ytemp(i) = ytemp(i) + del(i) 
     enddo
     x = x + h 


     !       call Timers_start ("derivs (burn)")

     call derivs(x,ytemp,yout) 

     !       call Timers_stop ("derivs (burn)")

  enddo

!!  take the last step
  do i=1,n 
     yout(i) = h * yout(i) - del(i)  
  enddo


  !      call Timers_start ("pchem_algebra (burn)")

!  do j = 1, n
!     do i = 1, n
!        av(i,j) = dmat(i,j)
!     enddo
!  enddo
  do i = 1,n 
    ! av(i,n+1) = yout(i)
    bb(i) = yout(i)
  enddo

 ! call pchem_gift(av,nmaxx,nmaxp1)
 ! call ludcmp(av,n,nmaxx,indx,d)  
!  call lubksb(av,n-1,nmaxx,indx,bb)
  call leqs(av,bb,n-1,nmaxx)
  do i = 1, n
   !  yout(i) = av(i,n+1)
   yout(i) = bb(i) 
  enddo
   bb(iELEC) = bb(iHP) + bb(iDP) + bb(iHEP) + 2.0*bb(iHEPP) + bb(iH2P) + bb(iHDP) - bb(iHM) - bb(iDM)

  !      call Timers_stop ("pchem_algebra (burn)")

  do i=1,n 
     yout(i) = ytemp(i) + yout(i) 
  enddo

!  print*, 'At end of Bader'

  return 

end subroutine pchem_baderStepGift

!!---------------------------------------------------------------------------------
!! NAME
!!
!!  pchem_pzExtr
!!
!! SYNOPSIS
!! 
!!  subroutine pchem_pzExtr(integer(IN) :: iest,
!!                       real(IN)    :: xest,
!!                       real(IN)    :: yest(:),
!!                       real(INOUT)   :: yz(:),
!!                       real(OUT)   :: dy(x:),
!!                       integer(IN) :: nv)
!!
!! DESCRIPTION
!!
!!  use polynomial extrapolation to evaluate nv functions at x=0 by fitting 
!!  a polynomial to a sequence of estimates with progressively smaller values 
!!  x=xest, and corresponding function vectors yest(1:nv). the call is number  
!!  iest in the sequence of calls. extrapolated function values are output as  
!!  yz(1:nv), and their estimated error is output as dy(1:nv) 
!!   
!! ARGUMENTS
!! 
!!  iest - integer(IN)      number in the sequence of calls
!!  xest - real(IN)         fitting location
!!  yest - real(IN)(1:nv)   function vectors at location xest
!!  yz -   real(INOUT)(1:nv)  extrapolated function values
!!  dy -   real(OUT)(1:nv)  estimated error of extrapolated function values
!!  nv -   integer(IN)      number of functions to evaluate
!!

subroutine pchem_pzExtr(iest,xest,yest,yz,dy,nv) 

  implicit none

!!  declare arguments
  integer, intent(IN)  :: nv, iest
  real, intent(IN)     :: xest
  real, intent(IN), dimension(nv) :: yest
  real, intent(OUT), dimension(nv) :: dy, yz
!!  real, intent(INOUT), dimension(nv) :: yz

!! local variables
  integer,save       :: j,k1 
  integer, parameter :: nmax=50
  integer, parameter :: imax=13
  real, save         ::  delta,f1,f2,q
  real, dimension(nmax), save     :: d
  real, dimension(nmax,imax),save :: qcol 
  real, dimension(imax), save     :: x(imax) 


!!  save current independent variables 
  x(iest) = xest 
  do j=1,nv 
     dy(j) = yest(j) 
     yz(j) = yest(j) 
  enddo

!!  store first estimate in first column 
  if (iest .eq. 1) then 
     do j=1,nv 
        qcol(j,1) = yest(j) 
     enddo

  else 
     do j=1,nv 
        d(j) = yest(j) 
     enddo
     do k1=1,iest-1 
        delta = 1.0e0/(x(iest-k1) - xest) 
        f1    = xest * delta 
        f2    = x(iest-k1) * delta 
!        print *,'delta: ', delta
!        print *,'iest: ', iest
!        print *,'k1: ', k1
!       print *,'x(iest-k1)-xest): ', x(iest-k1)-xest
      !!   propagate tableu 1 diagonal more 
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
 
end subroutine pchem_pzExtr

