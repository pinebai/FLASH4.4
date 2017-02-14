!!****ih* source/physics/sourceTerms/Burn/BurnIntegrate/bnIntegrate_interface
!!
!!  SYNOPSIS
!!
!!   use bnIntegrate_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the BurnIntegrate module that defines its
!! public interfaces.
!!
!!***
Module bnIntegrate_interface 

  !! These routines are found in the physics/sourceTerms/Burn/BurnIntegrate
  !! They are all used for integrating the ODE of nuclear burning

  interface
     subroutine bn_netIntegrate(start,stptry,stpmin,stopp,bc,  & 
          &                  eps,dxsav,kmax,   & 
          &                  xrk,yrk,xphys,yphys,xlogi,ylogi,  & 
          &                  nok,nbad,kount,odescal,iprint,  & 
          &                  derivs,jakob,bjakob,steper) 

       implicit none

       !! These first 3 interfaces taken from bnNetwork_interface.  Can't apparently have circular calls
       interface 
          !! derivs routine as a dummy argument expands to routine bn_network
          subroutine derivs(tt,y,dydt)
            real, intent(IN)  :: tt
            real, intent(INOUT) ::y(*)
            real, intent(OUT) :: dydt(*)
          end subroutine derivs
       end interface
       interface
          subroutine jakob(tt,y,dfdy,nzo,ndummy)
            implicit none
            integer, intent(IN) :: nzo, ndummy
            real, intent(IN)    :: tt
            real, intent(INOUT) :: y(*)
            real, intent(OUT)   :: dfdy(*)
          end subroutine jakob
       end interface
       interface
          subroutine bjakob(iloc,jloc,nzo,np)
            implicit none
            integer, intent(IN)  ::   iloc(*),jloc(*),np
            integer, intent(OUT) ::   nzo
          end subroutine bjakob
       end interface
       interface
          ! This interface is from below.....
          subroutine steper(y,dydx,nv,x,htry,eps,yscal,hdid,hnext, & 
               &                       derivs,jakob,bjakob)
            implicit none
            interface
               subroutine derivs(tt,y,dydt)
                 real, intent(IN)  :: tt
                 real, intent(INOUT) ::y(*)
                 real, intent(OUT) :: dydt(*)
               end subroutine derivs
            end interface
            interface
               subroutine jakob(tt,y,dfdy,nzo,ndummy)
                 implicit none
                 integer, intent(IN) :: nzo, ndummy
                 real, intent(IN)    :: tt
                 real, intent(INOUT) :: y(*)
                 real, intent(OUT)   :: dfdy(*)
               end subroutine jakob
            end interface
            interface
               subroutine bjakob(iloc,jloc,nzo,np)
                 implicit none
                 integer, intent(IN)  ::   iloc(*),jloc(*),np
                 integer, intent(OUT) ::   nzo
               end subroutine bjakob
            end interface
            integer, intent(IN) :: nv
            real, intent(INOUT) :: y(nv)
            real, intent(IN)    :: dydx(nv), yscal(nv), htry, eps
            real, intent(OUT)   :: hdid, hnext
            real, intent(INOUT) :: x
          end subroutine steper
       end interface

       integer, intent(IN)  :: xphys,yphys,xlogi,ylogi
       integer, intent(IN)  :: kmax, iprint
       real, intent(IN)     :: odescal, dxsav, eps
       real, intent(IN)     :: start, stptry, stpmin, stopp
       real, intent(INOUT), dimension(yphys) :: bc
       integer, intent(OUT) :: nok, nbad, kount
       real, intent(OUT), dimension(xphys)       :: xrk
       real, intent(OUT), dimension(yphys,xphys) :: yrk



       !! bn_netIntegrate calls four external functions derivs, jakob, bjakob, steper
       !! Their interfaces are given below.  They should correspond to (in the case of Aprox13)
       !! derivs = bn_network = aprox13  -- note these were the approximate names in flash2
       !! jakob  = bn_        = saprox13
       !! bjakob = bn_        = baprox13
       !! steper = bn_baderMa28, for example

     end subroutine bn_netIntegrate

  end interface

  interface
     subroutine bn_baderStepMa28(y,dydx,dfdy,nmax,n,xs,htot,nnstep,yout,  & 
          &                  nzo,a,naij,ivect,jvect,jloc,ikeep,iw,w,flag,  & 
          &                  derivs) 
       implicit none
       interface
          subroutine derivs(tt,y,dydt)
            real, intent(IN)  :: tt
            real, intent(INOUT) ::y(*)
            real, intent(OUT) :: dydt(*)
          end subroutine derivs
       end interface
       integer, intent(IN) :: n, naij, nmax
       integer, intent(INOUT) :: ikeep(1), jloc(naij)  ! because LBR can't figure out what it does
       integer, intent(INOUT) :: iw(1), flag !used by ma28 LU decomposition
       real, intent(OUT)   :: a(naij)
       real, intent(IN)    :: y(n), dydx(n), dfdy(naij), xs, htot
       integer, intent(IN) :: nnstep, nzo, ivect(naij), jvect(naij)
       real, intent(INOUT) :: yout(n), w(1)
     end subroutine bn_baderStepMa28
  end interface

  interface
     subroutine bn_baderStepGift(y,dydx,dfdy,nmax,n,xs,htot,nnstep,yout,  & 
          &                      derivs) 
       implicit none
       interface
          subroutine derivs(tt,y,dydt)
            real, intent(IN)  :: tt
            real, intent(INOUT) ::y(*)
            real, intent(OUT) :: dydt(*)
          end subroutine derivs
       end interface
       integer, intent(IN) :: n, nmax
       integer, intent(IN) :: nnstep
       real, intent(INOUT) :: yout(n)
       real, intent(IN)    :: xs, htot, y(n), dydx(n), dfdy(nmax,nmax)
     end subroutine bn_baderStepGift
  end interface

  interface
     subroutine bn_pzExtr(iest,xest,yest,yz,dy,nv) 
       implicit none
       integer, intent(IN)  :: nv, iest
       real, intent(IN)     :: xest
       real, intent(IN), dimension(nv) :: yest
       real, intent(OUT), dimension(nv) :: dy, yz
     end subroutine bn_pzExtr
  end interface

!-------------------------------------------------------------------------
  interface
     subroutine bn_baderMa28(y,dydx,nv,x,htry,eps,yscal,hdid,hnext, & 
          &                   derivs,jakob,bjakob)
       implicit none
       interface
          subroutine derivs(tt,y,dydt)
            real, intent(IN)  :: tt
            real, intent(INOUT) ::y(*)
            real, intent(OUT) :: dydt(*)
          end subroutine derivs
       end interface
       interface
          subroutine jakob(tt,y,dfdy,nzo,ndummy)
            implicit none
            integer, intent(IN) :: nzo, ndummy
            real, intent(IN)    :: tt
            real, intent(INOUT) :: y(*)
            real, intent(OUT)   :: dfdy(*)
          end subroutine jakob
       end interface
       interface
          subroutine bjakob(iloc,jloc,nzo,np)
            implicit none
            integer, intent(IN)  ::   iloc(*),jloc(*),np
            integer, intent(OUT) ::   nzo
          end subroutine bjakob
       end interface
       integer, intent(IN) :: nv
       real, intent(IN)    :: dydx(nv), yscal(nv), htry, eps
       real, intent(INOUT) :: x, y(nv)
       real, intent(OUT)   :: hdid, hnext
     end subroutine bn_baderMa28
  end interface

  interface
     subroutine bn_baderGift(y,dydx,nv,x,htry,eps,yscal,hdid,hnext, & 
          &                       derivs,jakob,bjakob)
       implicit none
       interface
          subroutine derivs(tt,y,dydt)
            real, intent(IN)  :: tt
            real, intent(INOUT) ::y(*)
            real, intent(OUT) :: dydt(*)
          end subroutine derivs
       end interface
       interface
          subroutine jakob(tt,y,dfdy,nzo,ndummy)
            implicit none
            integer, intent(IN) :: nzo, ndummy
            real, intent(IN)    :: tt
            real, intent(INOUT) :: y(*)
            real, intent(OUT)   :: dfdy(*)
          end subroutine jakob
       end interface
       interface
          subroutine bjakob(iloc,jloc,nzo,np)
            implicit none
            integer, intent(IN)  ::   iloc(*),jloc(*),np
            integer, intent(OUT) ::   nzo
          end subroutine bjakob
       end interface
       integer, intent(IN) :: nv
       real, intent(IN)    :: dydx(nv), yscal(nv), htry, eps
       real, intent(INOUT) :: x, y(nv)
       real, intent(OUT)   :: hdid, hnext
     end subroutine bn_baderGift
  end interface


  interface
     subroutine bn_rosenMa28(y,dydx,n,x,htry,eps,yscal,hdid,hnext,  & 
          &                      derivs,jakob,bjakob) 
       implicit none
       interface
          subroutine derivs(tt,y,dydt)
            real, intent(IN)  :: tt
            real, intent(INOUT) ::y(*)
            real, intent(OUT) :: dydt(*)
          end subroutine derivs
       end interface
       interface
          subroutine jakob(tt,y,dfdy,nzo,ndummy)
            implicit none
            integer, intent(IN) :: nzo, ndummy
            real, intent(IN)    :: tt
            real, intent(INOUT) :: y(*)
            real, intent(OUT)   :: dfdy(*)
          end subroutine jakob
       end interface
       interface
          subroutine bjakob(iloc,jloc,nzo,np)
            implicit none
            integer, intent(IN)  ::   iloc(*),jloc(*),np
            integer, intent(OUT) ::   nzo
          end subroutine bjakob
       end interface
       integer, intent(IN) :: n
       real, intent(IN)    :: yscal(n), dydx(n), htry, eps
       real, intent(INOUT) :: x, y(n)
       real, intent(OUT)   :: hdid, hnext
     end subroutine bn_rosenMa28
  end interface

  interface
     subroutine bn_rosenGift(y,dydx,n,x,htry,eps,yscal,hdid,hnext,  & 
          &                      derivs,jakob,bjakob) 
       implicit none
       interface
          subroutine derivs(tt,y,dydt)
            real, intent(IN)  :: tt
            real, intent(INOUT) ::y(*)
            real, intent(OUT) :: dydt(*)
          end subroutine derivs
       end interface
       interface
          subroutine jakob(tt,y,dfdy,nzo,ndummy)
            implicit none
            integer, intent(IN) :: nzo, ndummy
            real, intent(IN)    :: tt
            real, intent(INOUT) :: y(*)
            real, intent(OUT)   :: dfdy(*)
          end subroutine jakob
       end interface
       interface
          subroutine bjakob(iloc,jloc,nzo,np)
            implicit none
            integer, intent(IN)  ::   iloc(*),jloc(*),np
            integer, intent(OUT) ::   nzo
          end subroutine bjakob
       end interface
       integer, intent(IN) :: n
       real, intent(IN)    :: yscal(n), dydx(n), htry, eps
       real, intent(INOUT) :: x, y(n)
       real, intent(OUT)   :: hdid, hnext
     end subroutine bn_rosenGift
  end interface

  !!------------------- This is the general dummy name
  !!   it expands to either bn_baderGift or bn_baderMa28 or bn_rosenGift or bn_rosenMa28
  interface
     subroutine steper(y,dydx,nv,x,htry,eps,yscal,hdid,hnext, & 
          &                       derivs,jakob,bjakob)
       implicit none
       interface
          subroutine derivs(tt,y,dydt)
            real, intent(IN)  :: tt
            real, intent(INOUT) ::y(*)
            real, intent(OUT) :: dydt(*)
          end subroutine derivs
       end interface
       interface
          subroutine jakob(tt,y,dfdy,nzo,ndummy)
            implicit none
            integer, intent(IN) :: nzo, ndummy
            real, intent(IN)    :: tt
            real, intent(INOUT) :: y(*)
            real, intent(OUT)   :: dfdy(*)
          end subroutine jakob
       end interface
       interface
          subroutine bjakob(iloc,jloc,nzo,np)
            implicit none
            integer, intent(IN)  ::   iloc(*),jloc(*),np
            integer, intent(OUT) ::   nzo
          end subroutine bjakob
       end interface
       integer, intent(IN) :: nv
       real, intent(IN)    :: dydx(nv), yscal(nv), htry, eps
       real, intent(INOUT) :: x, y(nv)
       real, intent(OUT)   :: hdid, hnext
     end subroutine steper
  end interface



end Module bnIntegrate_interface
