!!****if* source/physics/Eos/EosNuclear/kernel/findtemp
!!
!! NAME
!!
!!  findtemp
!!
!! SYNOPSIS
!!
!!  call findtemp(:: lr,
!!                 :: lt0,
!!                 :: y,
!!                 :: epsin,
!!                 :: keyerrt,
!!                 :: rfeps)
!!
!! DESCRIPTION
!!
!! Finds the new temperature
!!
!! ARGUMENTS
!!
!!   lr : 
!!
!!   lt0 : 
!!
!!   y : 
!!
!!   epsin : 
!!
!!   keyerrt : 
!!
!!   rfeps : 
!!
!!
!!
!!***

subroutine findtemp(lr,lt0,y,epsin,keyerrt,rfeps)

  use eosmodule

  implicit none

  real :: lr,lt0,y,epsin
  real :: eps,lt,ldt
  real :: tol
  real :: d1,d2,d3
  real :: eps0,eps1,lt1

  real :: ltn,ltmax,ltmin
  real :: tinput,rfeps

  integer :: rl = 0
  integer :: itmax,i,keyerrt
  integer :: ii,jj,kk
  
  keyerrt=0

  tol=rfeps ! need to find energy to less than 1 in 10^-10
  itmax=20 ! use at most 20 iterations, then bomb

  lt=lt0
  lt1=lt 

  eps0=epsin
  eps1=eps0
  

  ltmax=logtemp(ntemp)
  ltmin=logtemp(1)

  ! Note: We are using Ewald's Lagrangian interpolator here!

  !preconditioning 1: do we already have the right temperature?
  call findthis(lr,lt,y,eps,alltables(:,:,:,2),d1,d2,d3)
  if (abs(eps-eps0).lt.tol*abs(eps0)) then
     return
  endif
  lt1=lt
  eps1=eps
 
!!$  write(*,"(i4,1P12E19.10)") 0,lr,lt0,y,lt,eps,eps0,abs(eps-eps0)/eps0,d2
!!$  write(*,"(i4,1P12E19.10)") 0,lr,lt0,y,ltmin,ltmax

  do i=1,itmax
!!$     print *, 'iteration', i
!!$     print *, eps0, eps
!!$     print *, lt0, lt, d2
     !d2 is the derivative deps/dlogtemp;
     ! let's try this:
     ldt = -(eps - eps0)/d2 
!     if(ldt.gt.0.d0) ltmin=lt
!     if(ldt.le.0.d0) ltmax=lt
     ltn = lt+ldt
!!$     if (ltn <= ltmin) then
!!$        print *, 'problem', i
!!$        print *, lt0
!!$        print *, lt, ltn, ldt
!!$        print *, eps, eps0
!!$        print *, abs(eps-eps0)/abs(eps0)
!!$        print *, d2
!!$        stop
!!$     endif
        
     ltn = min(ltn,ltmax)
     ltn = max(ltn,ltmin)
 
     lt1=lt
     lt=ltn
     eps1=eps
!!$     write(*,"(i4,1P12E19.10)") i,lr,lt0,y,lt,eps,eps0,abs(eps-eps0)/eps0,d2,ldt
     call findthis(lr,lt,y,eps,alltables(:,:,:,2),d1,d2,d3)
!     call findthis(lr,lt,y,eps,epst,d1,d2,d3)
     if (abs(eps - eps0).lt.tol*abs(eps0)) then
!!$        write(*,"(1P12E19.10)") tol,abs(eps-eps0)/eps0
        exit
     endif
     !setup new d2

     ! if we are closer than 10^-2  to the 
     ! root (eps-eps0)=0, we are switching to 
     ! the secant method, since the table is rather coarse and the
     ! derivatives may be garbage.
     if(abs(eps-eps0).lt.1.0d-3*abs(eps0)) then
        d2 = (eps-eps1)/(lt-lt1)
     endif
!     if(i.ge.10) then
!       write(*,*) "EOS: Did not converge in findtemp!"
!       write(*,*) "rl,logrho,logtemp0,ye,lt,eps,eps0,abs(eps-eps0)/eps0"
!     if(i.gt.5) then
!       write(*,"(i4,1P10E22.14)") i,lr,lt0,y,lt,eps,eps0,abs(eps-eps0)/eps0
!     endif
  enddo


 if(i.ge.itmax) then
    !keyerrt=667
    ! calling bisection
!!$    print *, 'calling bisection:'
!!$    write(*,*) "i,rl,logrho,logtemp0,lt,ye,eps,eps0,abs(eps-eps0)/eps0"
!!$    write(*,"(i4,i4,1P10E19.10)") i,rl,lr,lt0,lt,y,eps,eps0,abs(eps-eps0)/eps0
    call bisection(lr,lt0,y,eps0,lt,alltables(:,:,:,2),keyerrt,1)
!    call bisection(lr,lt0,y,eps0,alltables(:,:,:,2),keyerrt,1)

    if(keyerrt.eq.667) then
       if(lt.ge.log10(t_max_hack)) then
          ! handling for too high temperatures
          lt = log10(t_max_hack)
          keyerrt=0
          goto 12
       else if(abs(lt-log10(t_max_hack))/log10(t_max_hack).lt.0.025d0) then
          lt0 = min(lt,log10(t_max_hack))
          keyerrt=0
          goto 12
       else if(lt <= ltmin) then
          lt0 = max(lt,ltmin)
          keyerrt=0
          goto 12
       else
          ! total failure
          write(*,*) "EOS: Did not converge in findtemp!"
          write(*,*) "i,rl,logrho,logtemp0,lt,ye,eps,eps0,abs(eps-eps0)/eps0"
!!          write(*,"(i4,i4,1P10E19.10)") i,rl,lr,lt0,lt,y,eps,eps0,abs(eps-eps0)/eps0
          write(*,*) i,rl,lr,lt0,lt,y,eps,eps0,abs(eps-eps0)/eps0
          write(*,*) "d2=", d2
          write(*,*) "Tried calling bisection... didn't help... :-/"
          write(*,*) "Bisection error: ",keyerrt
       endif
    endif

    lt0=min(lt,log10(t_max_hack))
    return
 endif

12 continue
 
  lt0=min(lt,log10(t_max_hack))


end subroutine findtemp

subroutine findtemp_entropy(lr,lt0,y,sin,keyerrt,rfeps)

! This routine finds the new temperature based
! on rho, Y_e, entropy

  use eosmodule

  implicit none

  real :: lr,lt0,y,sin
  real :: s,lt,ldt
  real :: tol
  real :: d1,d2,d3
  real :: s0,s1,lt1

  real :: ltn,ltmax,ltmin
  real :: tinput,rfeps

  integer :: rl = 0
  integer itmax,i,keyerrt
  integer ii,jj,kk

  keyerrt=0

  tol=rfeps ! need to find energy to less than 1 in 10^-10
  itmax=20 ! use at most 20 iterations, then bomb

  lt=lt0
  lt1=lt 

  s0=sin
  s1=s0

  ltmax=logtemp(ntemp)
  ltmin=logtemp(1)

  !preconditioning 1: do we already have the right temperature?
  call findthis(lr,lt,y,s,alltables(:,:,:,3),d1,d2,d3)
  if (abs(s-s0).lt.tol*abs(s0)) then
     return
  endif
  lt1=lt
  s1=s
 

  do i=1,itmax
     !d2 is the derivative ds/dlogtemp;
     ldt = -(s - s0)/d2 
     ltn = lt+ldt
     ltn = min(ltn,ltmax)
     ltn = max(ltn,ltmin)
     lt1=lt
     lt=ltn
     s1=s
     call findthis(lr,lt,y,s,alltables(:,:,:,3),d1,d2,d3)
     if (abs(s - s0).lt.tol*abs(s0)) then
       exit
     endif
     !setup new d2

     ! if we are closer than 10^-2  to the 
     ! root (eps-eps0)=0, we are switching to 
     ! the secant method, since the table is rather coarse and the
     ! derivatives may be garbage.
     if(abs(s-s0).lt.1.0d-3*abs(s0)) then
        d2 = (s-s1)/(lt-lt1)
     endif
  enddo


 if(i.ge.itmax) then
    keyerrt=667
    call bisection(lr,lt0,y,s0,lt,alltables(:,:,:,3),keyerrt,2)
!    call bisection(lr,lt0,y,s0,alltables(:,:,:,3),keyerrt,2)
    if(keyerrt.eq.667) then
          write(*,*) "EOS: Did not converge in findtemp_entropy!"
          write(*,*) "rl,logrho,logtemp0,ye,lt,s,s0,abs(s-s0)/s0"
!!          write(*,"(i4,i4,1P10E19.10)") i,rl,lr,lt0,y,lt,s,s0,abs(s-s0)/s0
          write(*,*) i,rl,lr,lt0,y,lt,s,s0,abs(s-s0)/s0
          write(*,*) "Tried calling bisection... didn't help... :-/"
          write(*,*) "Bisection error: ",keyerrt
    endif
    
    lt0=lt
    return
 endif


  lt0=lt


end subroutine findtemp_entropy

subroutine findtemp_press(lr,lt0,y,pin,keyerrt,rfeps)

! This routine finds the new temperature based
! on rho, Y_e, entropy

  use eosmodule

  implicit none

  real :: lr,lt0,y,pin
  real :: p,lt,ldt
  real :: tol
  real :: d1,d2,d3
  real :: p0,p1,lt1

  real :: ltn,ltmax,ltmin
  real :: tinput,rfeps

  integer :: rl = 0
  integer itmax,i,keyerrt
  integer ii,jj,kk

  keyerrt=0

  tol=rfeps ! need to find energy to less than 1 in 10^-10
  itmax=20 ! use at most 20 iterations, then bomb

  lt=lt0
  lt1=lt 

  p0=pin
  p1=p0

  ltmax=logtemp(ntemp)
  ltmin=logtemp(1)

  !preconditioning 1: do we already have the right temperature?
  call findthis(lr,lt,y,p,alltables(:,:,:,1),d1,d2,d3)
  if (abs(p-p0).lt.tol*abs(p0)) then
     return
  endif
  lt1=lt
  p1=p
 
  do i=1,itmax
     !d2 is the derivative dp/dlogtemp;
     ldt = -(p - p0)/d2 
     ltn = lt+ldt
     ltn = min(ltn,ltmax)
     ltn = max(ltn,ltmin)
     lt1=lt
     lt=ltn
     p1=p
     call findthis(lr,lt,y,p,alltables(:,:,:,1),d1,d2,d3)
     if (abs(p - p0).lt.tol*abs(p0)) then
       exit
     endif
     !setup new d2

     ! if we are closer than 10^-2  to the 
     ! root (eps-eps0)=0, we are switching to 
     ! the secant method, since the table is rather coarse and the
     ! derivatives may be garbage.
     if(abs(p-p0).lt.1.0d-3*abs(p0)) then
        d2 = (p-p1)/(lt-lt1)
     endif
  enddo


 if(i.ge.itmax) then
    keyerrt=667
    call bisection(lr,lt0,y,p0,lt,alltables(:,:,:,1),keyerrt,2)
!    call bisection(lr,lt0,y,p0,alltables(:,:,:,1),keyerrt,2)
    if(keyerrt.eq.667) then
          write(*,*) "EOS: Did not converge in findtemp_press!"
          write(*,*) "rl,logrho,logtemp0,ye,lt,p,p0,abs(p-p0)/p0"
          write(*,"(i4,i4,10E19.10)") i,rl,lr,lt0,y,lt,p,p0,abs(p-p0)/p0
          write(*,*) "Tried calling bisection... didn't help... :-/"
          write(*,*) "Bisection error: ",keyerrt
    endif
    
    lt0=lt
    return
 endif


  lt0=lt


end subroutine findtemp_press
