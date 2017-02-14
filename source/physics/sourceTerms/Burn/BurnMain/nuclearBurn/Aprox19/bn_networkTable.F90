!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19/bn_networkTable
!!
!! NAME
!!  
!!  bn_networkTable
!!
!!
!! SYNOPSIS
!! 
!!  call bn_networkTable()
!!
!!  
!! DESCRIPTION
!!
!!  routine networkTable generates the raw reaction rates for routine aprox13
!!  by using table interpolation instead of analytical expressions.
!!  a cubic polynomial is hardwired for speed
!!
!!***

subroutine bn_networkTable()

#include "Flash.h"

   
  use Burn_dataEOS, ONLY: btemp, den
  use Burn_data
  use bn_dataAprox19

  implicit none

  !..uses tables instead of analytical expressions to evaluate the 
  !..raw reaction rates. a cubic polynomial is hardwired for speed.

  integer          i,j,imax,iat,mp,ifirst
  parameter        (mp = 4)
  real             tlo,thi,tstp,den_sav,btemp_sav,     &
       &                 x,x1,x2,x3,x4,a,b,c,d,e,f,g,h,p,q,  &
       &                 alfa,beta,gama,delt
  data             ifirst/0/


  !..make the table
  real :: entropy, dst, dsd
  if (ifirst .eq. 0) then
     ifirst = 1

     !..set the log temperature loop limits
     !..use 120 points per decade
     imax = 481
     if (imax .gt. nrattab) stop 'imax too small in bn_networkTable'
     tlo  = 6.0e0
     thi  = 10.0e0
     tstp = (thi - tlo)/float(imax-1)

     !..save the input
     btemp_sav = btemp
     den_sav   = den

     !..form the table
     den = 1.0e0
     do i=1,imax
        btemp = tlo + float(i-1)*tstp
        btemp = 10.0e0**(btemp)
        call bn_networkRates
        ttab(i) = btemp
        do j=1,nrat
           rattab(j,i) = ratraw(j)
        enddo
     enddo

     !..restore the input
     den  = den_sav
     btemp = btemp_sav
  end if


  !..normal execution starts here
  !..set the density dependence vector
  dtab(irpp)    = den
  dtab(irheng)  = den
  dtab(irhegn)  = 1.0e0
  dtab(irhng)   = den
  dtab(irdgn)   = 1.0e0
  dtab(irdpg)   = den
  dtab(irhegp)  = 1.0e0
  dtab(ir33)    = den
  dtab(ir34)    = den
  dtab(ircpg)   = den 
  dtab(irnpg)   = den 
  dtab(ifa)     = 1.0e0
  dtab(ifg)     = 1.0e0
  dtab(iropg)   = den
  dtab(irnag)   = den
  dtab(ircag)   = den
  dtab(iroga)   = 1.0e0
  dtab(ir3a)    = den*den
  dtab(irg3a)   = 1.0e0
  dtab(ir1212)  = den
  dtab(ir1216)  = den
  dtab(ir1616)  = den
  dtab(iroag)   = den
  dtab(irnega)  = 1.0e0
  dtab(irneag)  = den
  dtab(irmgga)  = 1.0e0
  dtab(irmgag)  = den
  dtab(irsiga)  = 1.0e0
  dtab(irmgap)  = den 
  dtab(iralpa)  = den
  dtab(iralpg)  = den
  dtab(irsigp)  = 1.0e0
  dtab(irsiag)  = den
  dtab(irsga)   = 1.0e0
  dtab(irppa)   = den
  dtab(irsiap)  = den
  dtab(irppg)   = den
  dtab(irsgp)   = 1.0e0
  dtab(irsag)   = den
  dtab(irarga)  = 1.0e0
  dtab(irsap)   = den
  dtab(irclpa)  = den
  dtab(irclpg)  = den
  dtab(irargp)  = 1.0e0
  dtab(irarag)  = den
  dtab(ircaga)  = 1.0e0
  dtab(irarap)  = den
  dtab(irkpa)   = den
  dtab(irkpg)   = den
  dtab(ircagp)  = 1.0e0
  dtab(ircaag)  = den
  dtab(irtiga)  = 1.0e0
  dtab(ircaap)  = den
  dtab(irscpa)  = den
  dtab(irscpg)  = den
  dtab(irtigp)  = 1.0e0
  dtab(irtiag)  = den
  dtab(ircrga)  = 1.0e0
  dtab(irtiap)  = den
  dtab(irvpa)   = den
  dtab(irvpg)   = den
  dtab(ircrgp)  = 1.0e0
  dtab(ircrag)  = den
  dtab(irfega)  = 1.0e0
  dtab(ircrap)  = den
  dtab(irmnpa)  = den
  dtab(irmnpg)  = den
  dtab(irfegp)  = 1.0e0
  dtab(irfeag)  = den
  dtab(irniga)  = 1.0e0
  dtab(irfeap)  = den
  dtab(ircopa)  = den
  dtab(ircopg)  = den
  dtab(irnigp)  = 1.0e0
  dtab(irfepg)  = den
  dtab(ircogp)  = 1.0e0
  dtab(ir52ng)  = den
  dtab(ir53gn)  = 1.0e0
  dtab(ir53ng)  = den
  dtab(ir54gn)  = 1.0e0
  dtab(irpen)   = 1.0e0
  dtab(ispen)   = 1.0e0
  dtab(irnep)   = 1.0e0
  dtab(isnep)   = 1.0e0
  dtab(irn56ec) = 1.0e0
  dtab(isn56ec) = 1.0e0



  !..hash locate the temperature
  iat = int((log10(btemp) - tlo)/tstp) + 1
  iat = max(1,min(iat - mp/2 + 1,imax - mp + 1))

  !..setup the lagrange interpolation coefficients for a cubic
  x  = btemp
  x1 = ttab(iat)
  x2 = ttab(iat+1)
  x3 = ttab(iat+2)
  x4 = ttab(iat+3)
  a  = x - x1
  b  = x - x2
  c  = x - x3
  d  = x - x4
  e  = x1 - x2
  f  = x1 - x3
  g  = x1 - x4
  h  = x2 - x3
  p  = x2 - x4
  q  = x3 - x4
  alfa =  b*c*d/(e*f*g)
  beta = -a*c*d/(e*h*p)
  gama =  a*b*d/(f*h*q)
  delt = -a*b*c/(g*p*q)

  !..crank off the raw reaction rates
  do j=1,nrat
     ratraw(j) = (alfa*rattab(j,iat)    &
          &            + beta*rattab(j,iat+1)  &
          &            + gama*rattab(j,iat+2)  &
          &            + delt*rattab(j,iat+3)  &
          &              ) * dtab(j)
  enddo
  return
end subroutine bn_networkTable







