!!****if* source/physics/sourceTerms/Burn/BurnIntegrate/bn_baderMa28
!!
!! NAME
!!   bn_baderMa28  
!!
!! SYNOPSIS
!!  subroutine bn_baderMa28(real(IN)       ::y(:),
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
!!  for sparse analytic jacobians, ma28 linear algebra 
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
!!  for the integration nuclear reaction networks with 2 choices for the linear algebra.
!!
!!  routine bn_baderMa28 drives a bader-deuflhard step with ma28 algebra
!!  routine bn_baderStepMa28 is a bader-deuflhard stepper with ma28 algebra
!!  routine bn_baderGift drives a bader-deuflhard step with gift algebra
!!  routine bn_baderStepGift is a bader-deuflhard stepper with ma28 algebra
!!  routine bn_pzExtr does extrapolations for any of the *Step* routines
!!
!!***
!!---------------------------------------------------------------------------------



subroutine bn_baderMa28(y,dydx,nv,x,htry,eps,yscal,hdid,hnext, & 
     &                   derivs,jakob,bjakob)


  use Burn_data, ONLY: aion
  use Driver_interface, ONLY : Driver_abortFlash
  ! can't use jakob interface; see notes in bnNetwork_interface for the mystery.
  use bnIntegrate_interface, ONLY: bn_baderStepMa28, bn_pzExtr

  implicit none

  interface   ! = bn_network
     subroutine derivs(tt,y,dydt)   !! == bn_network
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


  !!  argument declarations
  external               jakob
  integer, intent(IN) :: nv
  real, intent(IN)    :: dydx(nv), yscal(nv), htry, eps
  real, intent(INOUT) :: x, y(nv)
  real, intent(OUT)   :: hdid, hnext



  !! local variables 
  logical, save    :: first = .true.
  logical, save    :: reduct
  integer, parameter  ::    nmax = 30, kmaxx=7, imax=kmaxx+1   
  integer, parameter  ::    naij=300, n5 = 5*nmax, n8=8*nmax
  integer, save :: iloc(naij),jloc(naij), & 
       &                 ivect(naij),jvect(naij),ikeep(n5),iw(n8), & 
       &                 flag,nzo,ifirst,ii, &
       &            nDummy ! used to make number of arguments of jakob consistent
  integer, save :: i,iq,k,kk,km,kmax,kopt,nvold,nseq(imax)
  real, save ::    eps1,epsold,errmax,fact,h,red,scale,xwork,wrkmin, & 
       &                 xest,xnew,a(imax),alf(kmaxx,kmaxx),err(kmaxx), & 
       &                 yerr(nmax),ysav(nmax),yseq(nmax),& 
       &                 dfdy(naij),amat(naij),w(nmax),u

  real, parameter     ::   safe1 = 0.25e0, safe2 = 0.7e0, redmax=1.0e-5, & 
       &                  redmin = 0.7e0, tiny = 1.0e-30, scalmx = 0.1e0

  data             epsold/-1.0e0/, nvold/-1/
  data             nseq /2, 6, 10, 14, 22, 34, 50, 70/
  data             ifirst/0/ 

  !!-----------------------------------------------------------------------------
  !!  get and copy the nonzero locations
  if (ifirst .eq. 0) then

     !       call Timers_start ("algebra (burn)")

     ifirst = 1

     call bjakob(iloc,jloc,nzo,naij)

     do i=1,nzo
        ivect(i) = iloc(i)
        jvect(i) = jloc(i)
     enddo

     !!   force the diagonal to be the pivot elements
     do i=1,nzo
        amat(i) = 1.0e-10
        if (ivect(i) .eq. jvect(i)) amat(i) = 1.0e0
     enddo

     u  = 0.1e0

     call ma28ad(nv,nzo,amat,naij,iloc,naij,jloc,u,ikeep,iw,w,flag)

     if (flag .lt. 0) then
        write(*,*) 'error in ma28ad flag',flag
        write(*,*) 'more than stpmax steps required in bn_netInt'  
        call Driver_abortFlash('ERROR in bn_netInt: too many steps required')
     end if

     !       call Timers_stop ("algebra (burn)")

  end if

  !!  a new tolerance or a new number , so reinitialize
  if (eps .ne. epsold  .or.  nv .ne. nvold) then
     hnext = -1.0e29
     xnew  = -1.0e29
     eps1  = safe1 * eps

     !!   compute the work coefficients a_k
     a(1)  = nseq(1) + 1
     do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
     enddo

     !!   compute alf(k,q)
     do iq=2,kmaxx
        do k=1,iq-1
           alf(k,iq) = eps1**((a(k+1) - a(iq+1)) / & 
                &               ((a(iq+1) - a(1) + 1.0e0) * (2*k + 1)))
        enddo
     enddo
     epsold = eps
     nvold  = nv

     !!   add cost of jacobians to work coefficients
     a(1)   = nv + a(1)
     do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
     enddo

     !!   determine optimal row number for convergence
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

  !!  get the sparse jacobian in dfdy
  !      call Timers_start ("jacobian (burn)")

  call jakob(x,y,dfdy,nzo,nDummy)

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
        write(*,*) 'stepsize too small in routine bn_baderMa28'
        write(*,110) xnew,x,k
110     format(1x,2e11.3,i6)

        do ii=1,nv
           write(*,111) y(ii),y(ii)*aion(ii),aion(ii),ii
111        format(1x,3e11.3,i6)
        enddo

        call Driver_abortFlash('ERROR in bn_baderMa28: stepsize too small')
     end if

     call bn_baderStepMa28(ysav,dydx,dfdy,nmax,nv,x,h,nseq(k),yseq, & 
          &             nzo,amat,naij,ivect,jvect,jloc,ikeep,iw,w,flag, & 
          &             derivs) 

     xest = (h/nseq(k))**2 
     call bn_pzExtr(k,xest,yseq,y,yerr,nv) 


     !!   compute normalized error estimate
     if (k .ne. 1) then
        errmax = tiny
        do i=1,nv
           errmax = max(errmax,abs(yerr(i)/yscal(i)))
        enddo
        errmax   = errmax/eps   
        km = k - 1
        err(km) = (errmax/safe1)**(1.0e0/(2*km+1))
     end if

     !!   in order window
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

end subroutine bn_baderMa28

!!---------------------------------------------------------------------------------

subroutine bn_baderStepMa28(y,dydx,dfdy,nmax,n,xs,htot,nnstep,yout,  & 
     &                  nzo,a,naij,ivect,jvect,jloc,ikeep,iw,w,flag,  & 
     &                  derivs) 

  use Driver_interface, ONLY : Driver_abortFlash
!  use bnIntegrate_interface, ONLY: derivs

  implicit none

!! This routine was simpr3
!!  
!!  an implicit midpoint stepper, for ma28 sparse linear algebra. 
!!   
!!  declare arguments
!!  NOTE jloc and ikeep are not used anywhere.....

  external               derivs 
  integer, intent(IN) :: n, naij, nmax
  integer, intent(IN) :: nnstep, nzo, ivect(naij), jvect(naij)
  integer, intent(INOUT) :: ikeep(1), jloc(naij)  ! because LBR can't figure out what it does
  integer, intent(INOUT) :: iw(1), flag !used by ma28 LU decomposition
  real, intent(IN)    :: y(n), dydx(n), dfdy(naij), xs, htot
  real, intent(INOUT) :: yout(n), w(1)
  real, intent(OUT)   :: a(naij)

!!  declare local variables
  integer, parameter ::     nmaxx=30 
  integer, save :: i,nn
  real, save ::    h,x,del(nmaxx),ytemp(nmaxx)

!!  stepsize this trip, and make the a matrix 
  h = htot/nnstep 
  do i=1,nzo 
     a(i) = -h * dfdy(i) 
     if (ivect(i) .eq. jvect(i)) a(i) = 1.0e0 + a(i) 
  enddo

!!  numeric decomp 
  !      call Timers_start ("algebra (burn)")

  call ma28bd(n,nzo,a,naij,ivect,jvect,jloc,ikeep,iw,w,flag) 

  !      call Timers_stop ("algebra (burn)")


  if (flag .lt. 0) then 
     write(*,*) 'error in ma28bd flag',flag 
     call Driver_abortFlash('ERROR in return from ma28bd')
  end if

!!  use yout as temporary storage; the first step 
  do i=1,n 
     yout(i) = h * dydx(i) 
  enddo

  !      call Timers_start ("algebra (burn)")

  call ma28cd(n,a,naij,jloc,ikeep,yout,w,1) 

  !      call Timers_stop ("algebra (burn)")

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

     !       call Timers_start ("algebra (burn)")

     call ma28cd(n,a,naij,jloc,ikeep,yout,w,1) 

     !       call Timers_stop ("algebra (burn)")

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

  !      call Timers_start ("algebra (burn)")

  call ma28cd(n,a,naij,jloc,ikeep,yout,w,1) 

  !      call Timers_stop ("algebra (burn)")

  do i=1,n 
     yout(i) = ytemp(i) + yout(i) 
  enddo
  return 


end subroutine bn_baderStepMa28

!!---------------------------------------------------------------------------------

subroutine bn_baderGift(y,dydx,nv,x,htry,eps,yscal,hdid,hnext, & 
     &                       derivs,jakob,bjakob)

  use Burn_data, ONLY: aion
  use Driver_interface, ONLY : Driver_abortFlash
  !  Ummm.... bit of a mystery why I can use the interfaces in bn_netIntRosen but not here.
  use bnIntegrate_interface, ONLY: bn_baderStepGift, bn_pzExtr

  implicit none

!!  for dense analytic jacobians, gift algebra
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

  interface   ! = bn_network
     subroutine derivs(tt,y,dydt)   !! == bn_network
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
  real, intent(OUT)   :: hdid, hnext


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

  call jakob(x,y,dfdy,nv,nmax)

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
110     format(1x,2e11.3,i6)

        do ii=1,nv
           write(*,111) y(ii),y(ii)*aion(ii),aion(ii),ii
111        format(1x,3e11.3,i6)
        enddo

        call Driver_abortFlash('ERROR in baderGift: stepsize too small')
     end if

   !!  the semi implicit midpoint rules for this sequence
     call bn_baderStepGift(ysav,dydx,dfdy,nmax,nv,x,h,nseq(k),yseq,derivs) 

   !!  extrapolate the error to zero
     xest = (h/nseq(k))**2 
     call bn_pzExtr(k,xest,yseq,y,yerr,nv)


   !!  compute normalized error estimate
     if (k .ne. 1) then
        errmax = tiny
        do i=1,nv
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
end subroutine bn_baderGift

!!---------------------------------------------------------------------------------


subroutine bn_baderStepGift(y,dydx,dfdy,nmax,n,xs,htot,nnstep,yout,  & 
     &                      derivs) 

  use bn_interface, ONLY: bn_gift

  implicit none
  
!!   
!!  an implicit midpoint stepper, for ma28 sparse linear algebra. 
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

!!  for the gift linear algebra
  integer, parameter  :: nmaxp1 = nmaxx + 1
  real, save          :: dmat(nmaxx,nmaxx), av(nmaxx,nmaxp1)


!!  stepsize this trip, and make the a matrix 
  h = htot/nnstep 

  !      call Timers_start ("algebra (burn)")

  do j=1,n
     do i=1,n
        dmat(i,j) = -h * dfdy(i,j) 
     enddo
  enddo
  do i=1,n
     dmat(i,i) = 1.0e0 + dmat(i,i) 
  end do

  !      call Timers_stop ("algebra (burn)")

!!  use yout as temporary storage; the first step 
  do i=1,n 
     yout(i) = h * dydx(i) 
  enddo


  !      call Timers_start ("algebra (burn)")

  do j = 1, n
     do i = 1, n
        av(i,j) = dmat(i,j)
     enddo
  enddo
  do i = 1,n 
     av(i,n+1) = yout(i)
  enddo
  call bn_gift(av,nmaxx,nmaxp1)
  do i = 1, n
     yout(i) = av(i,n+1) 
  enddo

  !      call Timers_stop ("algebra (burn)")

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

     !       call Timers_start ("algebra (burn)")

     do j = 1, n
        do i = 1, n
           av(i,j) = dmat(i,j)
        enddo
     enddo
     do i = 1,n 
        av(i,n+1) = yout(i)
     enddo
     call bn_gift(av,nmaxx,nmaxp1)
     do i = 1, n
        yout(i) = av(i,n+1) 
     enddo

     !       call Timers_stop ("algebra (burn)")


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


  !      call Timers_start ("algebra (burn)")

  do j = 1, n
     do i = 1, n
        av(i,j) = dmat(i,j)
     enddo
  enddo
  do i = 1,n 
     av(i,n+1) = yout(i)
  enddo
  call bn_gift(av,nmaxx,nmaxp1)
  do i = 1, n
     yout(i) = av(i,n+1) 
  enddo

  !      call Timers_stop ("algebra (burn)")

  do i=1,n 
     yout(i) = ytemp(i) + yout(i) 
  enddo

  return 

end subroutine bn_baderStepGift

!!---------------------------------------------------------------------------------
!! NAME
!!
!!  bn_pzExtr
!!
!! SYNOPSIS
!! 
!!  subroutine bn_pzExtr(integer(IN) :: iest,
!!                       real(IN)    :: xest,
!!                       real(IN)    :: yest(:),
!!                       real(INOUT)   :: yz(:),
!!                       real(OUT)   :: dy(:),
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

subroutine bn_pzExtr(iest,xest,yest,yz,dy,nv) 

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
 
end subroutine bn_pzExtr

