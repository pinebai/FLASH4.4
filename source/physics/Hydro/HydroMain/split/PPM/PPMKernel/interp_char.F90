!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/interp_char
!!
!! NAME
!! 
!!  interp_char
!!
!! SYNOPSIS
!!
!!  call interp_char(integer(IN) :: numIntCells,
!!                   integer(IN) :: numCells,
!!                   real(OUT)   :: al(numCells), 
!!                   real(IN)    :: a(numCells), 
!!                   real(OUT)   :: ar(numCells), 
!!                   real(IN)    :: coeff1(numCells), 
!!                   real(IN)    :: coeff2(numCells), 
!!                   real(IN)    :: coeff3(numCells), 
!!                   real(IN)    :: coeff4(numCells), 
!!                   real(IN)    :: coeff5(numCells))
!!
!! DESCRIPTION
!!
!!  Interpolate interface values and monotonize.
!!  We solve eqns (7) and (8) from Colella and Sekora, or equivalantly, eqn (46) combined
!!  with eqn (38) and (39) from Stone et al.
!!
!! ARGUMENTS
!!
!! numIntCells :
!! numCells :
!! al :
!! a :
!! ar :
!! coeff1 :
!! coeff2 :
!! coeff3 :
!! coeff4 :
!! coeff5 :
!!
!!
!! SIDE EFFECTS
!!
!!  Modifies array hy_dela, which is referenced in subroutine detect.
!!
!!***

#define HY_FASTLEFT 1
#define HY_SLOWLEFT 2
#define HY_ENTROPY  3
#define HY_SLOWRGHT 4
#define HY_FASTRGHT 5

#define HY_DENS 1
#define HY_VELX 2
#define HY_VELY 3
#define HY_VELZ 4
#define HY_PRES 5
#define HY_GAMC 6
#define HY_GAME 7

subroutine interp_char(sweepDir,numIntCells, numCells, &
                       rhol, rho, rhor, &
                       ul,   u,   ur,   &
                       utl,  ut,  utr,  &
                       uttl, utt, uttr, &
                       pl,   p,   pr,   &
                       gamc, game,       &
                       coeff1, coeff2, coeff3, coeff4, coeff5)


  use Hydro_data, ONLY:  hy_dela
                     
  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: sweepDir,numIntCells, numCells
  real, intent(IN), DIMENSION(numCells)  ::  coeff1, coeff2, coeff3, coeff4, coeff5
  real, intent(IN), DIMENSION(numCells)  ::  rho,u,ut,utt,p,gamc,game
  real, intent(OUT), DIMENSION(numCells) ::  rhol,rhor,ul,ur,utl,utr,uttl,uttr,pl,pr
  real, dimension(numCells) :: scrch1, scrch2, scrch3, scrch4
  integer :: i, numIntCells5, numIntCells6, numIntCells8


  integer :: numVars,nWave
  real, dimension(5,5) :: reigL,reig0,reigR,leigL,leig0,leigR
  real, dimension(5)   :: lambdaL,lambda0,lambdaR,delbarL,delbar0,delbarR
  real, dimension(7)   :: wLL,wL,w0,wR,wRR


  numIntCells5 = numIntCells + 5
  numIntCells6 = numIntCells + 6
  numIntCells8 = numIntCells + 8

  numVars = 7 ! rho,u,ut,utt,p,gamc,game

  rhol=0.
  rhor=0.
  ul=0.
  ur=0.
  utl=0.
  utr=0.
  uttl=0.
  uttr=0.
  pl=0.
  pr=0.

  do i = 2, numIntCells8
     scrch1(i) = rho(i) - rho(i-1)
     scrch2(i) = abs ( scrch1(i) + scrch1(i) )
     scrch4(i) = sign (1.e00, scrch1(i))
  end do


  ! apply Eq. 1.8 of Colella & Woodward -- guarantee that a(i+1/2) lies
  ! between a(i) and a(i+1)
    
  do i = 2, numIntCells6+1
     hy_dela(i)   = coeff1(i) * scrch1(i+1) + coeff2(i) * scrch1(i)
       
     if (hy_dela(i) .LT. 0.e0) then
        scrch3(i) = -1.e0
     else
        scrch3(i) = 1.e0
     endif

     ! This hy_dela for density is used in detect.F90 for contact steepening
     ! and is initialized here.
     ! Note: This hy_dela should be initialized before calling detect.F90
     !       for contact steepening
     hy_dela(i)   = min(abs(hy_dela(i)), scrch2(i), scrch2(i+1))* scrch3(i)

     if (-scrch4(i)*scrch4(i+1) >= 0.e0) hy_dela(i) = 0.e0
       
  end do


  ! compute some common factors
  !!* We solve eqns (7) and (8) in Colella & Sekora, instead of solving eqn 1.6 from
  !!  Colella & Woodward (or eqn 12 from Colella & Sekora)
  !!* We apply limitings to characteristic variables, NOT to primitive variables

  do i = 3, numIntCells8-2

     select case(sweepDir)
     case(1)
        wLL=(/rho(i-2),u(i-2),ut(i-2),utt(i-2),p(i-2),gamc(i-2),game(i-2)/)
        wL =(/rho(i-1),u(i-1),ut(i-1),utt(i-1),p(i-1),gamc(i-1),game(i-1)/)
        w0 =(/rho( i ),u( i ),ut( i ),utt( i ),p( i ),gamc( i ),game( i )/)
        wR =(/rho(i+1),u(i+1),ut(i+1),utt(i+1),p(i+1),gamc(i+1),game(i+1)/)
        wRR=(/rho(i+2),u(i+2),ut(i+2),utt(i+2),p(i+2),gamc(i+2),game(i+2)/)

     case(2)
        wLL=(/rho(i-2),ut(i-2),u(i-2),utt(i-2),p(i-2),gamc(i-2),game(i-2)/)
        wL =(/rho(i-1),ut(i-1),u(i-1),utt(i-1),p(i-1),gamc(i-1),game(i-1)/)
        w0 =(/rho( i ),ut( i ),u( i ),utt( i ),p( i ),gamc( i ),game( i )/)
        wR =(/rho(i+1),ut(i+1),u(i+1),utt(i+1),p(i+1),gamc(i+1),game(i+1)/)
        wRR=(/rho(i+2),ut(i+2),u(i+2),utt(i+2),p(i+2),gamc(i+2),game(i+2)/)

     case(3)
        wLL=(/rho(i-2),ut(i-2),utt(i-2),u(i-2),p(i-2),gamc(i-2),game(i-2)/)
        wL =(/rho(i-1),ut(i-1),utt(i-1),u(i-1),p(i-1),gamc(i-1),game(i-1)/)
        w0 =(/rho( i ),ut( i ),utt( i ),u( i ),p( i ),gamc( i ),game( i )/)
        wR =(/rho(i+1),ut(i+1),utt(i+1),u(i+1),p(i+1),gamc(i+1),game(i+1)/)
        wRR=(/rho(i+2),ut(i+2),utt(i+2),u(i+2),p(i+2),gamc(i+2),game(i+2)/)
     end select

     call eigensystem(sweepDir,wL,reigL,leigL,lambdaL)
     call eigensystem(sweepDir,w0,reig0,leig0,lambda0)
     call eigensystem(sweepDir,wR,reigR,leigR,lambdaR)


     ! One can use a different slope limiter - default is MC limiter
     do nWave = HY_FASTLEFT,HY_FASTRGHT
        delbarL(nWave) = &
             mc(dot_product(leigL(nWave,HY_DENS:HY_PRES),wL(HY_DENS:HY_PRES)-wLL(HY_DENS:HY_PRES)),&
                dot_product(leigL(nWave,HY_DENS:HY_PRES),w0(HY_DENS:HY_PRES)-wL (HY_DENS:HY_PRES)))


        delbar0(nWave) = &
             mc(dot_product(leig0(nWave,HY_DENS:HY_PRES),w0(HY_DENS:HY_PRES)-wL(HY_DENS:HY_PRES)),&
                dot_product(leig0(nWave,HY_DENS:HY_PRES),wR(HY_DENS:HY_PRES)-w0(HY_DENS:HY_PRES)))


        delbarR(nWave) = &
             mc(dot_product(leigR(nWave,HY_DENS:HY_PRES),wR (HY_DENS:HY_PRES)-w0(HY_DENS:HY_PRES)),&
                dot_product(leigR(nWave,HY_DENS:HY_PRES),wRR(HY_DENS:HY_PRES)-wR(HY_DENS:HY_PRES)))
     enddo


     delbarL(HY_DENS:HY_PRES) = reigL(HY_DENS:HY_PRES,HY_FASTLEFT)*delbarL(HY_FASTLEFT)+&
                                reigL(HY_DENS:HY_PRES,HY_SLOWLEFT)*delbarL(HY_SLOWLEFT)+&
                                reigL(HY_DENS:HY_PRES,HY_ENTROPY) *delbarL(HY_ENTROPY) +&
                                reigL(HY_DENS:HY_PRES,HY_SLOWRGHT)*delbarL(HY_SLOWRGHT)+&
                                reigL(HY_DENS:HY_PRES,HY_FASTRGHT)*delbarL(HY_FASTRGHT)


     delbar0(HY_DENS:HY_PRES) = reig0(HY_DENS:HY_PRES,HY_FASTLEFT)*delbar0(HY_FASTLEFT)+&
                                reig0(HY_DENS:HY_PRES,HY_SLOWLEFT)*delbar0(HY_SLOWLEFT)+&
                                reig0(HY_DENS:HY_PRES,HY_ENTROPY) *delbar0(HY_ENTROPY) +&
                                reig0(HY_DENS:HY_PRES,HY_SLOWRGHT)*delbar0(HY_SLOWRGHT)+&
                                reig0(HY_DENS:HY_PRES,HY_FASTRGHT)*delbar0(HY_FASTRGHT)

     delbarR(HY_DENS:HY_PRES) = reigR(HY_DENS:HY_PRES,HY_FASTLEFT)*delbarR(HY_FASTLEFT)+&
                                reigR(HY_DENS:HY_PRES,HY_SLOWLEFT)*delbarR(HY_SLOWLEFT)+&
                                reigR(HY_DENS:HY_PRES,HY_ENTROPY) *delbarR(HY_ENTROPY) +&
                                reigR(HY_DENS:HY_PRES,HY_SLOWRGHT)*delbarR(HY_SLOWRGHT)+&
                                reigR(HY_DENS:HY_PRES,HY_FASTRGHT)*delbarR(HY_FASTRGHT)



     ! Parabolic interpolation at the left and right cell interfaces
     ! Colella-Woodward Eqn 1.9, Sekora-Colella Eqn 7, Stone et al Eqn 46
     wLL(HY_DENS:HY_PRES) = 0.5*(w0(HY_DENS:HY_PRES)+wL(HY_DENS:HY_PRES)) &
          - (delbar0(HY_DENS:HY_PRES)-delbarL(HY_DENS:HY_PRES))/6.

     wRR(HY_DENS:HY_PRES) = 0.5*(w0(HY_DENS:HY_PRES)+wR(HY_DENS:HY_PRES)) &
          - (delbarR(HY_DENS:HY_PRES)-delbar0(HY_DENS:HY_PRES))/6.


     ! Ensure that the interpolated values lie between the cell-centered values
     do nWave=HY_DENS,HY_PRES
        wLL(nWave) = max(min(w0(nWave),wL(nWave)),wLL(nWave))
        wLL(nWave) = min(max(w0(nWave),wL(nWave)),wLL(nWave))
        wRR(nWave) = max(min(w0(nWave),wR(nWave)),wRR(nWave))
        wRR(nWave) = min(max(w0(nWave),wR(nWave)),wRR(nWave))
     enddo


     ! Prepare outputs
     rhol(i) = wLL(HY_DENS)
     pl(i)   = wLL(HY_PRES)

     rhor(i) = wRR(HY_DENS)
     pr(i)   = wRR(HY_PRES)

     select case(sweepDir)
     case(1)
        ul(i)   = wLL(HY_VELX)
        utl(i)  = wLL(HY_VELY)
        uttl(i) = wLL(HY_VELZ)

        ur(i)   = wRR(HY_VELX)
        utr(i)  = wRR(HY_VELY)
        uttr(i) = wRR(HY_VELZ)

     case(2)
        ul(i)   = wLL(HY_VELY)
        utl(i)  = wLL(HY_VELX)
        uttl(i) = wLL(HY_VELZ)

        ur(i)   = wRR(HY_VELY)
        utr(i)  = wRR(HY_VELX)
        uttr(i) = wRR(HY_VELZ)
     case(3)
        ul(i)   = wLL(HY_VELZ)
        utl(i)  = wLL(HY_VELX)
        uttl(i) = wLL(HY_VELY)

        ur(i)   = wRR(HY_VELZ)
        utr(i)  = wRR(HY_VELX)
        uttr(i) = wRR(HY_VELY)
     end select

  enddo



  contains

    subroutine eigensystem(sweepDir,W,reig,leig,lambda)

      implicit none
      integer, intent(IN) :: sweepDir
      real, intent(IN),dimension(7) :: W
      real, intent(OUT), dimension(5,5) :: reig,leig
      real, intent(OUT), dimension(5)   :: lambda


      real :: a  ! sound speed
      real :: uN ! normal velocity
      real :: dinv

      select case(sweepDir)
      case(1)
         uN = W(HY_VELX)
      case(2)
         uN = W(HY_VELY)
      case(3)
         uN = W(HY_VELZ)
      end select

      ! sound velocity
      a = sqrt(W(HY_GAMC)*W(HY_PRES)/W(HY_DENS))

      ! eigen value
      lambda(HY_FASTLEFT:HY_FASTRGHT) = (/uN-a,uN,uN,uN,uN+a/)


      dinv = 1./W(HY_DENS)

      reig = 0.
      leig = 0.

      ! right eigenvectors: reig(1:nVars,1:nWave)
      reig(HY_DENS,HY_FASTLEFT) = W(HY_DENS)
      reig(HY_DENS,HY_FASTRGHT) = W(HY_DENS)

      select case(sweepDir)
      case(1)
         reig(HY_VELX,HY_FASTLEFT)= -a
         reig(HY_VELX,HY_FASTRGHT)=  a

         reig(HY_VELY,HY_SLOWLEFT)= -dinv
         reig(HY_VELZ,HY_SLOWRGHT)=  dinv

      case(2)
         reig(HY_VELY,HY_FASTLEFT)= -a
         reig(HY_VELY,HY_FASTRGHT)=  a

         reig(HY_VELX,HY_SLOWLEFT)=  dinv
         reig(HY_VELZ,HY_SLOWRGHT)= -dinv

      case(3)
         reig(HY_VELZ,HY_FASTLEFT)= -a
         reig(HY_VELZ,HY_FASTRGHT)=  a

         reig(HY_VELX,HY_SLOWLEFT)= -dinv
         reig(HY_VELY,HY_SLOWRGHT)=  dinv
      end select
      reig(HY_PRES,HY_FASTLEFT) =  W(HY_DENS)*a*a
      reig(HY_PRES,HY_FASTRGHT) =  W(HY_DENS)*a*a


      ! entropy wave
      reig(HY_DENS,HY_ENTROPY) = 1.


      ! left eigenvectors: leig(1:nWave,1:nVars)
      select case(sweepDir)
      case(1)
         leig(HY_FASTLEFT,HY_VELX)= -a
         leig(HY_FASTRGHT,HY_VELX)=  a

         leig(HY_SLOWLEFT,HY_VELY)= -W(HY_DENS)
         leig(HY_SLOWRGHT,HY_VELZ)=  W(HY_DENS)


      case(2)
         leig(HY_FASTLEFT,HY_VELY)= -a
         leig(HY_FASTRGHT,HY_VELY)=  a

         leig(HY_SLOWLEFT,HY_VELX)=  W(HY_DENS)
         leig(HY_SLOWRGHT,HY_VELZ)= -W(HY_DENS)
      case(3)
         leig(HY_FASTLEFT,HY_VELZ)= -a
         leig(HY_FASTRGHT,HY_VELZ)=  a

         leig(HY_SLOWLEFT,HY_VELX)= -W(HY_DENS)
         leig(HY_SLOWRGHT,HY_VELY)=  W(HY_DENS)
      end select

      leig(HY_FASTLEFT, HY_PRES) =  dinv
      leig(HY_FASTRGHT, HY_PRES) =  dinv


      ! scale fast waves with 1/2*a^2
      leig(HY_FASTLEFT,:) = 0.5/a**2*leig(HY_FASTLEFT,:)
      leig(HY_FASTRGHT,:) = 0.5/a**2*leig(HY_FASTRGHT,:)

      ! entropy wave
      leig(HY_ENTROPY,HY_DENS) = 1.
      leig(HY_ENTROPY,HY_PRES) =-1./a**2



    end subroutine eigensystem


    function vanLeer(a,b)
      implicit none
      real :: a,b,vanLeer
      if (a*b <=0.) then
         vanLeer=0.
      else
         vanLeer=2.*a*b/(a+b)
      endif
    end function vanLeer


    function mc(a,b)
      implicit none
      real :: a,b,mc
      mc = (sign(1.,a)+sign(1.,b))*min(abs(a),.25*abs(a+b),abs(b))
    end function mc


    function minmod(a,b)
      implicit none
      real :: a,b,minmod
      minmod=.5 * (sign(1.,a) + sign(1.,b))*min(abs(a),abs(b))
    end function minmod


    function signum(x)
      implicit none
      real :: x,signum
      signum=sign(.5,x)-sign(.5,-x)
    end function signum



  end subroutine interp_char






