!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox13/bn_networkDenseJakob
!!
!! NAME
!!  
!!  bn_networkDenseJakob
!!
!! SYNOPSIS
!! 
!!  call bn_networkDenseJakob (real, intent(IN) :: tt,
!!                          real, intent(INOUT) :: y(:),
!!                            real, intent(OUT) :: dfdy(nphys,nphys),
!!                          integer, intent(IN) :: nlog,
!!                          integer, intent(IN) :: nphys)
!!
!!  
!! DESCRIPTION
!!
!!  routine networkDenseJakob sets up the dense aprox13 jacobian 
!!  
!!
!! ARGUMENTS
!!
!! tt :
!! y :
!! dfdy :
!! nlog :
!! nphys :
!!
!!***


subroutine bn_networkDenseJakob(tt,y,dfdy,nlog,nphys)   

!! This is the only one LBR will do to include by name.  Essentially all
!!   of the bn_dataAprox13 is used...
  use Burn_data, ONLY: ratdum, ife52, icr48, iti44, ica40, iar36, is32, &
       isi28, img24, ine20, io16, ic12, ihe4, ini56
  use bn_dataAprox13, ONLY: &
       irscpg,irvpg, irvpa, irmnpg, irmnpa, ircopg, ircopa, irfeag, ircrag,  &
       irtiag, ircaag, irarag, irsag, irsiag, irmgag, irneag, iroag, ircag,  &
       ir3a, irfeap, ircrap, irtiap, ircaap, irarap, irsap, irsiap, irmgap,  &
       irg3a, ir1216, ir1212, iroga, ir1616, irnega, irmgga, irsigp, irsiga, &
       irsgp, irsga, irargp, irarga, ircagp, ircaga, irtigp, irtiga, ircrgp, &
       ircrga, irfegp, irfega, irnigp, irniga, iralpg, iralpa, irppg, irppa, &
       irclpg, irclpa, irkpg, irkpa, irscpa
 

  implicit none

  !    save

#include "constants.h"
#include "Flash.h"


!!  declare
  integer, intent(IN)  :: nlog, nphys
  real, intent(IN)     :: tt
  real, intent(INOUT), dimension(*)           :: y
  real, intent(OUT), dimension(nphys,nphys) :: dfdy
  integer          i,j
  real         r1,s1,t1,u1,v1,w1,x1,y1


  !!  zero the jacobian
  do j=1,nlog
     do i=1,nlog
        dfdy(i,j) = 0.0e0
     enddo
  enddo


  !!  positive definite mass fractions
  do i=1,NSPECIES
     y(i) = min(1.0e0,max(y(i),1.0e-30))
  enddo


  !!  branching ratios for dummy proton links
  r1     = ratdum(iralpa)/(ratdum(iralpa) + ratdum(iralpg))
  s1     = ratdum(irppa)/(ratdum(irppa) + ratdum(irppg))
  t1     = ratdum(irclpa)/(ratdum(irclpa) + ratdum(irclpg))
  u1     = ratdum(irkpa)/(ratdum(irkpa) + ratdum(irkpg))
  v1     = ratdum(irscpa)/(ratdum(irscpa) + ratdum(irscpg))
  w1     = ratdum(irvpa)/(ratdum(irvpa) + ratdum(irvpg))
  x1     = ratdum(irmnpa)/(ratdum(irmnpa) + ratdum(irmnpg))
  y1     = ratdum(ircopa)/(ratdum(ircopa) + ratdum(ircopg))


  !!  he4 jacobian elements
  dfdy(ihe4,ihe4)  = -9.0e0 * y(ihe4) * y(ihe4) * ratdum(ir3a)    &
       &                   - y(ic12)  * ratdum(ircag)       &
       &                   - y(io16)  * ratdum(iroag)       &
       &                   - y(ine20) * ratdum(irneag)  &
       &                   - y(img24) * ratdum(irmgag) &
       &                   - y(isi28) * ratdum(irsiag)  &
       &                   - y(is32)  * ratdum(irsag)  &
       &                   - y(iar36) * ratdum(irarag) &
       &                   - y(ica40) * ratdum(ircaag) &
       &                   - y(iti44) * ratdum(irtiag) &
       &                   - y(icr48) * ratdum(ircrag) &
       &                   - y(ife52) * ratdum(irfeag) 

  dfdy(ihe4,ihe4)  = dfdy(ihe4,ihe4) &
       &                   - y(img24) * ratdum(irmgap) * (1.0e0-r1) &
       &                   - y(isi28) * ratdum(irsiap) * (1.0e0-s1) &
       &                   - y(is32) * ratdum(irsap)   * (1.0e0-t1) &
       &                   - y(iar36) * ratdum(irarap) * (1.0e0-u1) &
       &                   - y(ica40) * ratdum(ircaap) * (1.0e0-v1)  &
       &                   - y(iti44) * ratdum(irtiap) * (1.0e0-w1) &
       &                   - y(icr48) * ratdum(ircrap) * (1.0e0-x1)  &
       &                   - y(ife52) * ratdum(irfeap) * (1.0e0-y1) 


  dfdy(ihe4,ic12)  = 2.0e0 * y(ic12) * ratdum(ir1212)    &
       &                   + 0.5e0 * y(io16) * ratdum(ir1216)    &
       &                   + 3.0e0 * ratdum(irg3a)   &
       &                   - y(ihe4) * ratdum(ircag) 

  dfdy(ihe4,io16)  = 0.5e0 * y(ic12) * ratdum(ir1216) &
       &                   + 1.12e0 * y(io16) * ratdum(ir1616)  &
       &                   + 0.68e0 * s1 * y(io16) * ratdum(ir1616)  &
       &                   + ratdum(iroga) &
       &                   - y(ihe4) * ratdum(iroag)          

  dfdy(ihe4,ine20) =  ratdum(irnega) &
       &                  - y(ihe4) * ratdum(irneag)

  dfdy(ihe4,img24) =   ratdum(irmgga) &
       &                   - y(ihe4) * ratdum(irmgag)  &
       &                   - y(ihe4) * ratdum(irmgap) * (1.0e0-r1)

  dfdy(ihe4,isi28) =   ratdum(irsiga) &
       &                   - y(ihe4) * ratdum(irsiag) &
       &                   - y(ihe4) * ratdum(irsiap) * (1.0e0-s1) &
       &                   + r1 * ratdum(irsigp) 


  dfdy(ihe4,is32)  =   ratdum(irsga) &
       &                   - y(ihe4) * ratdum(irsag) &
       &                   - y(ihe4) * ratdum(irsap) * (1.0e0-t1) &
       &                   + s1 * ratdum(irsgp)

  dfdy(ihe4,iar36) =   ratdum(irarga) &
       &                   - y(ihe4) * ratdum(irarag) &
       &                   - y(ihe4) * ratdum(irarap) * (1.0e0-u1) &
       &                   + t1 * ratdum(irargp)

  dfdy(ihe4,ica40) =   ratdum(ircaga) &
       &                   - y(ihe4) * ratdum(ircaag) &
       &                   - y(ihe4) * ratdum(ircaap) * (1.0e0-v1) &
       &                   + u1 * ratdum(ircagp)

  dfdy(ihe4,iti44) =   ratdum(irtiga) &
       &                   - y(ihe4) * ratdum(irtiag) &
       &                   - y(ihe4) * ratdum(irtiap) * (1.0e0-w1) &
       &                   + v1 * ratdum(irtigp)

  dfdy(ihe4,icr48) =   ratdum(ircrga) &
       &                   - y(ihe4) * ratdum(ircrag) &
       &                   - y(ihe4) * ratdum(ircrap) * (1.0e0-x1) &
       &                   + w1 * ratdum(ircrgp)

  dfdy(ihe4,ife52) =   ratdum(irfega) &
       &                   - y(ihe4) * ratdum(irfeag)  &
       &                   - y(ihe4) * ratdum(irfeap) * (1.0e0-y1) &
       &                   + x1 * ratdum(irfegp) 

  dfdy(ihe4,ini56) =   ratdum(irniga) &
       &                   + y1 * ratdum(irnigp)  


  !!  c12 jacobian elements
  dfdy(ic12,ihe4) = 3.0e0 * y(ihe4) * y(ihe4) * ratdum(ir3a) &
       &                 - y(ic12) * ratdum(ircag)   

  dfdy(ic12,ic12) = -4.0e0 * y(ic12) * ratdum(ir1212)  &
       &                  - y(io16) * ratdum(ir1216)   &
       &                  - ratdum(irg3a) &
       &                  - y(ihe4) * ratdum(ircag)

  dfdy(ic12,io16) = -y(ic12) * ratdum(ir1216)    &
       &                 + ratdum(iroga)



  !!  16o jacobian elements
  dfdy(io16,ihe4) = y(ic12)*ratdum(ircag) &
       &                - y(io16)*ratdum(iroag) 

  dfdy(io16,ic12) = -y(io16)*ratdum(ir1216)  &
       &                 + y(ihe4)*ratdum(ircag)

  dfdy(io16,io16) = - y(ic12) * ratdum(ir1216)   &
       &                  - 4.0e0 * y(io16) * ratdum(ir1616)  &
       &                  - y(ihe4) * ratdum(iroag) &
       &                  - ratdum(iroga) 

  dfdy(io16,ine20) = ratdum(irnega)



  !!  20ne jacobian elements
  dfdy(ine20,ihe4)  = y(io16) * ratdum(iroag)   &
       &                  - y(ine20) * ratdum(irneag)       

  dfdy(ine20,ic12)  = 2.0e0 * y(ic12) * ratdum(ir1212)  

  dfdy(ine20,io16)  = y(ihe4) * ratdum(iroag)

  dfdy(ine20,ine20) = -y(ihe4) * ratdum(irneag)  &
       &                   - ratdum(irnega)

  dfdy(ine20,img24) = ratdum(irmgga)



  !!  24mg jacobian elements
  dfdy(img24,ihe4)  = y(ine20) * ratdum(irneag) &
       &                   -y(img24) * ratdum(irmgag)  &
       &                   -y(img24) * ratdum(irmgap) * (1.0e0-r1) 

  dfdy(img24,ic12)  = 0.5e0 * y(io16) * ratdum(ir1216)

  dfdy(img24,io16)  = 0.5e0 * y(ic12) * ratdum(ir1216)

  dfdy(img24,ine20) = y(ihe4) * ratdum(irneag)

  dfdy(img24,img24) = -y(ihe4) * ratdum(irmgag)  &
       &                   - ratdum(irmgga)  &
       &                   - y(ihe4) * ratdum(irmgap) * (1.0e0-r1)  

  dfdy(img24,isi28) = ratdum(irsiga)  &
       &                  + r1 * ratdum(irsigp)



  !!  28si jacobian elements
  dfdy(isi28,ihe4)  = y(img24) * ratdum(irmgag) &
       &                  - y(isi28) * ratdum(irsiag)  &
       &                  + y(img24) * ratdum(irmgap) * (1.0e0-r1) &
       &                  - y(isi28) * ratdum(irsiap) * (1.0e0-s1)

  dfdy(isi28,ic12)  = 0.5e0 * y(io16) * ratdum(ir1216) 

  dfdy(isi28,io16)  =   0.5e0 * y(ic12) * ratdum(ir1216)  &
       &                    + 1.12e0 * y(io16) * ratdum(ir1616)   &
       &                    + 0.68e0 * y(io16) * s1 * ratdum(ir1616)     

  dfdy(isi28,img24) = y(ihe4) * ratdum(irmgag)  &
       &                  + y(ihe4) * ratdum(irmgap) * (1.0e0-r1)

  dfdy(isi28,isi28) = -y(ihe4) * ratdum(irsiag)  &
       &                   - ratdum(irsiga) &
       &                   - r1 * ratdum(irsigp)  &
       &                   - y(ihe4) * ratdum(irsiap) * (1.0e0-s1)

  dfdy(isi28,is32)  = ratdum(irsga)  &
       &                  + s1 * ratdum(irsgp)



  !!  32s jacobian elements
  dfdy(is32,ihe4)  = y(isi28) * ratdum(irsiag) &
       &                 - y(is32) * ratdum(irsag) &
       &                 + y(isi28) * ratdum(irsiap) * (1.0e0-s1) &
       &                 - y(is32) * ratdum(irsap) * (1.0e0-t1)

  dfdy(is32,io16)  = 0.68e0 * y(io16) * ratdum(ir1616) * (1.0e0-s1) &
       &                   + 0.2e0 * y(io16) * ratdum(ir1616) 

  dfdy(is32,isi28) = y(ihe4) * ratdum(irsiag)  &
       &                  + y(ihe4) * ratdum(irsiap) * (1.0e0-s1)

  dfdy(is32,is32)  = -y(ihe4) * ratdum(irsag) &
       &                  - ratdum(irsga)  &
       &                  - s1 * ratdum(irsgp) &
       &                  - y(ihe4) * ratdum(irsap) * (1.0e0-t1)

  dfdy(is32,iar36) = ratdum(irarga)  &
       &                 + t1 * ratdum(irargp)



  !!  36ar jacobian elements
  dfdy(iar36,ihe4)  = y(is32)  * ratdum(irsag) &
       &                  - y(iar36) * ratdum(irarag) &
       &                  + y(is32)  * ratdum(irsap) * (1.0e0-t1)  &
       &                  - y(iar36) * ratdum(irarap) * (1.0e0-u1)

  dfdy(iar36,is32)  = y(ihe4) * ratdum(irsag) &
       &                   + y(ihe4) * ratdum(irsap) * (1.0e0-t1)

  dfdy(iar36,iar36) = -y(ihe4) * ratdum(irarag)  &
       &                   - ratdum(irarga)  &
       &                   - t1 * ratdum(irargp) &
       &                   - y(ihe4) * ratdum(irarap) * (1.0e0-u1)

  dfdy(iar36,ica40) = ratdum(ircaga)  &
       &                  + ratdum(ircagp) * u1



  !!  40ca jacobian elements
  dfdy(ica40,ihe4)   = y(iar36) * ratdum(irarag) &
       &                   - y(ica40) * ratdum(ircaag)  &
       &                   + y(iar36) * ratdum(irarap)*(1.0e0-u1) &
       &                   - y(ica40) * ratdum(ircaap)*(1.0e0-v1)

  dfdy(ica40,iar36)  = y(ihe4) * ratdum(irarag) &
       &                   + y(ihe4) * ratdum(irarap)*(1.0e0-u1)

  dfdy(ica40,ica40)  = -y(ihe4) * ratdum(ircaag) &
       &                    - ratdum(ircaga)  &
       &                    - ratdum(ircagp) * u1 &
       &                    - y(ihe4) * ratdum(ircaap)*(1.0e0-v1)

  dfdy(ica40,iti44)  = ratdum(irtiga)  &
       &                    + ratdum(irtigp) * v1



  !!  44ti jacobian elements
  dfdy(iti44,ihe4)   = y(ica40) * ratdum(ircaag) &
       &                   - y(iti44) * ratdum(irtiag) &
       &                   + y(ica40) * ratdum(ircaap)*(1.0e0-v1) &
       &                   - y(iti44) * ratdum(irtiap)*(1.0e0-w1)

  dfdy(iti44,ica40)  = y(ihe4) * ratdum(ircaag) &
       &                   + y(ihe4) * ratdum(ircaap)*(1.0e0-v1)

  dfdy(iti44,iti44)  = -y(ihe4) * ratdum(irtiag) &
       &                    - ratdum(irtiga)  &
       &                    - v1 * ratdum(irtigp) &
       &                    - y(ihe4) * ratdum(irtiap)*(1.0e0-w1)

  dfdy(iti44,icr48)  = ratdum(ircrga)  &
       &                   + w1 * ratdum(ircrgp)



  !!  48cr jacobian elements
  dfdy(icr48,ihe4)  = y(iti44) * ratdum(irtiag) &
       &                  - y(icr48) * ratdum(ircrag) &
       &                  + y(iti44) * ratdum(irtiap)*(1.0e0-w1) &
       &                  - y(icr48) * ratdum(ircrap)*(1.0e0-x1)

  dfdy(icr48,iti44) = y(ihe4) * ratdum(irtiag) &
       &                  + y(ihe4) * ratdum(irtiap)*(1.0e0-w1)

  dfdy(icr48,icr48) = -y(ihe4) * ratdum(ircrag)  &
       &                   - ratdum(ircrga)  &
       &                   - w1 * ratdum(ircrgp) &
       &                   - y(ihe4) * ratdum(ircrap)*(1.0e0-x1)

  dfdy(icr48,ife52) = ratdum(irfega)  &
       &                  + x1 * ratdum(irfegp)



  !!  52fe jacobian elements
  dfdy(ife52,ihe4)  = y(icr48) * ratdum(ircrag) &
       &                  - y(ife52) * ratdum(irfeag) &
       &                  + y(icr48) * ratdum(ircrap) * (1.0e0-x1)  &
       &                  - y(ife52) * ratdum(irfeap) * (1.0e0-y1)

  dfdy(ife52,icr48) = y(ihe4) * ratdum(ircrag) &
       &                  + y(ihe4) * ratdum(ircrap) * (1.0e0-x1) 

  dfdy(ife52,ife52) = - y(ihe4) * ratdum(irfeag) &
       &                    - ratdum(irfega) &
       &                    - x1 * ratdum(irfegp) &
       &                    - y(ihe4) * ratdum(irfeap) * (1.0e0-y1)

  dfdy(ife52,ini56) = ratdum(irniga) &
       &                  + y1 * ratdum(irnigp)



  !!  56ni jacobian elements
  dfdy(ini56,ihe4)  = y(ife52) * ratdum(irfeag) &
       &                  + y(ife52) * ratdum(irfeap) * (1.0e0-y1)

  dfdy(ini56,ife52) = y(ihe4) * ratdum(irfeag) &
       &                  + y(ihe4) * ratdum(irfeap) * (1.0e0-y1)

  dfdy(ini56,ini56) = -ratdum(irniga) &
       &                   - y1 * ratdum(irnigp)


  return
end   subroutine bn_networkDenseJakob
