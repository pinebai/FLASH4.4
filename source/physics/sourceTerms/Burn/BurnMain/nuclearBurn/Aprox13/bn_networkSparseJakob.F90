!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox13/bn_networkSparseJakob
!!
!! NAME
!!
!!  bn_networkSparseJakob
!!
!! SYNOPSIS
!! 
!!  call bn_networkSparseJakob(real, intent(IN)    :: tt,
!!                             real, intent(IN)    :: y(:),
!!                             real, intent(OUT)   :: dfdy(nphys,nphys),
!!                             integer, intent(IN) :: nzo,
!!                             integer, intent(IN) :: nDummy)
!!
!!  
!! DESCRIPTION
!!
!!  routine networkSparseJakob sets up the sparse aprox13 jacobian 
!!  input is tt (irrelevant here) and the abundances y(*). 
!!  output is the jacobian dfdy(nzo).
!!
!!  This routine (and bn_networkDenseJakob) can be called as an external
!!    'jakob'
!!
!! ARGUMENTS
!!
!!   tt -- not used in this subroutine
!!   y  -- abundances
!!   dfdy -- jacobian
!!   nzo -- size of jacobian
!!   nDummy -- a fake argument to make the number of arguments the same
!!
!!***


subroutine bn_networkSparseJakob(tt,y,dfdy,nzo,nDummy)

  use Burn_data
  use Driver_interface, ONLY : Driver_abortFlash
  ! for communication with networkSparsePointers
  use bn_dataNetworkSize, ONLY : neloc, eloc, nterms
  use bn_dataAprox13
 
  implicit none

  !    save


  !!  declare arguments and local variables
  ! nDummy brings number of variables up to be the same as bn_networkDenseJakob
  integer, intent(IN) :: nzo, nDummy  
  real, intent(IN)    :: tt
  real, intent(INOUT) :: y(*)
  real, intent(OUT)   :: dfdy(*)

  integer          i,nt,iat
  real             r1,s1,t1,u1,v1,w1,x1,y1,a1

  !!  communicate with the jacobian builder networkSparsePointers, now in bn_dataNetworkSize
!!  integer          neloc
!!  parameter        (neloc=65)
!!  integer          eloc(neloc),nterms
!!  common /elca13/  eloc,nterms

  !!  initialize
  nt = 0
  do i=1,nzo
     dfdy(i) = 0.0e0
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
  !!  d(he4)/d(he4)
  a1 = -9.0e0 * y(ihe4) * y(ihe4) * ratdum(ir3a)                     &
       &    - y(ic12)  * ratdum(ircag)                                     &
       &    - y(io16)  * ratdum(iroag)                                     &
       &    - y(ine20) * ratdum(irneag)                                    &
       &    - y(img24) * ratdum(irmgag)                                    &
       &    - y(isi28) * ratdum(irsiag)                                    &
       &    - y(is32)  * ratdum(irsag)                                     &
       &    - y(iar36) * ratdum(irarag)                                    &
       &    - y(ica40) * ratdum(ircaag)                                    &
       &    - y(iti44) * ratdum(irtiag)                                    &
       &    - y(icr48) * ratdum(ircrag)                                    &
       &    - y(ife52) * ratdum(irfeag) 
  a1 =  a1                                                          &
       &    - y(img24) * ratdum(irmgap) * (1.0e0-r1)                       &
       &    - y(isi28) * ratdum(irsiap) * (1.0e0-s1)                       &
       &    - y(is32) * ratdum(irsap)   * (1.0e0-t1)                       &
       &    - y(iar36) * ratdum(irarap) * (1.0e0-u1)                       &
       &    - y(ica40) * ratdum(ircaap) * (1.0e0-v1)                       &
       &    - y(iti44) * ratdum(irtiap) * (1.0e0-w1)                       &
       &    - y(icr48) * ratdum(ircrap) * (1.0e0-x1)                       &
       &    - y(ife52) * ratdum(irfeap) * (1.0e0-y1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(c12)
  a1 =  2.0e0 * y(ic12) * ratdum(ir1212)                             &
       &    + 0.5e0 * y(io16) * ratdum(ir1216)                             &
       &    + 3.0e0 * ratdum(irg3a)                                        &
       &    - y(ihe4) * ratdum(ircag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(o16)
  a1 =  0.5e0 * y(ic12) * ratdum(ir1216)                             &
       &    + 1.12e0 * y(io16) * ratdum(ir1616)                            &
       &    + 0.68e0 * s1 * y(io16) * ratdum(ir1616)                       &
       &    + ratdum(iroga)                                                &
       &    - y(ihe4) * ratdum(iroag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(ne20)
  a1  = ratdum(irnega)                                               &
       &   - y(ihe4) * ratdum(irneag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(mg24)
  a1  =  ratdum(irmgga)                                              &
       &     - y(ihe4) * ratdum(irmgag)                                    &
       &     - y(ihe4) * ratdum(irmgap) * (1.0e0-r1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(si28)
  a1  =  ratdum(irsiga)                                              &
       &    - y(ihe4) * ratdum(irsiag)                                     &
       &    - y(ihe4) * ratdum(irsiap) * (1.0e0-s1)                        &
       &    + r1 * ratdum(irsigp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(s32)
  a1  =  ratdum(irsga)                                               &
       &    - y(ihe4) * ratdum(irsag)                                      &
       &    - y(ihe4) * ratdum(irsap) * (1.0e0-t1)                         &
       &    + s1 * ratdum(irsgp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(ar36)
  a1  =  ratdum(irarga)                                              &
       &    - y(ihe4) * ratdum(irarag)                                     &
       &    - y(ihe4) * ratdum(irarap) * (1.0e0-u1)                        &
       &    + t1 * ratdum(irargp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(ca40)
  a1  =  ratdum(ircaga)                                              &
       &    - y(ihe4) * ratdum(ircaag)                                     &
       &    - y(ihe4) * ratdum(ircaap) * (1.0e0-v1)                        &
       &    + u1 * ratdum(ircagp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(ti44)
  a1  =  ratdum(irtiga)                                              &
       &    - y(ihe4) * ratdum(irtiag)                                     &
       &    - y(ihe4) * ratdum(irtiap) * (1.0e0-w1)                        &
       &    + v1 * ratdum(irtigp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(cr48)
  a1  =  ratdum(ircrga)                                              &
       &    - y(ihe4) * ratdum(ircrag)                                     &
       &    - y(ihe4) * ratdum(ircrap) * (1.0e0-x1)                        &
       &    + w1 * ratdum(ircrgp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(fe52)
  a1  = ratdum(irfega)                                               &
       &    - y(ihe4) * ratdum(irfeag)                                     &
       &    - y(ihe4) * ratdum(irfeap) * (1.0e0-y1)                        &
       &    + x1 * ratdum(irfegp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(he4)/d(ni56)
  a1  =  ratdum(irniga)                                              &
       &     + y1 * ratdum(irnigp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  c12 jacobian elements
  !!  d(c12)/d(he4)
  a1  =  3.0e0 * y(ihe4) * y(ihe4) * ratdum(ir3a)                    &
       &     - y(ic12) * ratdum(ircag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(c12)/d(c12)
  a1  = -4.0e0 * y(ic12) * ratdum(ir1212)                            &
       &     - y(io16) * ratdum(ir1216)                                    &
       &     - ratdum(irg3a)                                               &
       &     - y(ihe4) * ratdum(ircag) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(c12)/d(o16)
  a1  = -y(ic12) * ratdum(ir1216)                                    &
       &     + ratdum(iroga)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  16o jacobian elements
  !!  d(o16)/d(he4)
  a1  = y(ic12)*ratdum(ircag)                                        &
       &    - y(io16)*ratdum(iroag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(o16)/d(c12)
  a1  = -y(io16)*ratdum(ir1216)                                      &
       &    + y(ihe4)*ratdum(ircag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(o16)/d(o16)
  a1  = - y(ic12) * ratdum(ir1216)                                   &
       &      - 4.0e0 * y(io16) * ratdum(ir1616)                           &
       &      - y(ihe4) * ratdum(iroag)                                    &
       &      - ratdum(iroga)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(o16)/d(n20)
  a1  = ratdum(irnega) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  20ne jacobian elements
  !!  d(ne20)/d(he4)
  a1  = y(io16) * ratdum(iroag)                                      &
       &    - y(ine20) * ratdum(irneag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ne20)/d(c12)
  a1  = 2.0e0 * y(ic12) * ratdum(ir1212)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ne20)/d(o16)
  a1  = y(ihe4) * ratdum(iroag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ne20)/d(ne20)
  a1  = -y(ihe4) * ratdum(irneag)                                    &
       &     - ratdum(irnega)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ne20)/d(mg24)
  a1  = ratdum(irmgga)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  24mg jacobian elements
  !!  d(mg24)/d(he4)
  a1   = y(ine20) * ratdum(irneag)                                   &
       &      -y(img24) * ratdum(irmgag)                                   &
       &      -y(img24) * ratdum(irmgap) * (1.0e0-r1) 
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(mg24)/d(c12)
  a1  = 0.5e0 * y(io16) * ratdum(ir1216) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(mg24)/d(o16)
  a1  = 0.5e0 * y(ic12) * ratdum(ir1216) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(mg24)/d(ne20)
  a1  = y(ihe4) * ratdum(irneag) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(mg24)/d(mg24)
  a1  = -y(ihe4) * ratdum(irmgag)                                    &
       &     - ratdum(irmgga)                                              &
       &     - y(ihe4) * ratdum(irmgap) * (1.0e0-r1)  
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(mg24)/d(si28)
  a1  = ratdum(irsiga)                                               &
       &    + r1 * ratdum(irsigp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  28si jacobian elements
  !!  d(si28)/d(he4)
  a1   =  y(img24) * ratdum(irmgag)                                  &
       &      - y(isi28) * ratdum(irsiag)                                  &
       &      + y(img24) * ratdum(irmgap) * (1.0e0-r1)                     &
       &      - y(isi28) * ratdum(irsiap) * (1.0e0-s1)
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(si28)/d(c12)
  a1  =  0.5e0 * y(io16) * ratdum(ir1216)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(si28)/d(o16)
  a1  =  0.5e0 * y(ic12) * ratdum(ir1216)                            &
       &     + 1.12e0 * y(io16) * ratdum(ir1616)                           &
       &     + 0.68e0 * y(io16) * s1 * ratdum(ir1616)     
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(si28)/d(mg24)
  a1  =  y(ihe4) * ratdum(irmgag)                                    &
       &      + y(ihe4) * ratdum(irmgap) * (1.0e0-r1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(si28)/d(si28)
  a1  = -y(ihe4) * ratdum(irsiag)                                    &
       &     - ratdum(irsiga)                                              &
       &     - r1 * ratdum(irsigp)                                         &
       &     - y(ihe4) * ratdum(irsiap) * (1.0e0-s1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(si28)/d(s32)
  a1  = ratdum(irsga)                                                &
       &    + s1 * ratdum(irsgp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  32s jacobian elements
  !!  d(s32)/d(he4)
  a1   = y(isi28) * ratdum(irsiag)                                   &
       &     - y(is32) * ratdum(irsag)                                     &
       &     + y(isi28) * ratdum(irsiap) * (1.0e0-s1)                      &
       &     - y(is32) * ratdum(irsap) * (1.0e0-t1)
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(s32)/d(o16)
  a1  =  0.68e0 * y(io16) * ratdum(ir1616) * (1.0e0-s1)              &
       &      + 0.2e0 * y(io16) * ratdum(ir1616) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(s32)/d(si28)
  a1  =  y(ihe4) * ratdum(irsiag)                                    &
       &     + y(ihe4) * ratdum(irsiap) * (1.0e0-s1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(s32)/d(s32)
  a1  = -y(ihe4) * ratdum(irsag)                                     &
       &     - ratdum(irsga)                                               &
       &     - s1 * ratdum(irsgp)                                          &
       &     - y(ihe4) * ratdum(irsap) * (1.0e0-t1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(s32)/d(ar36)
  a1  =  ratdum(irarga)                                              &
       &     + t1 * ratdum(irargp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  36ar jacobian elements
  !!  d(ar36)/d(he4)
  a1   =  y(is32)  * ratdum(irsag) &
       &      - y(iar36) * ratdum(irarag) &
       &      + y(is32)  * ratdum(irsap) * (1.0e0-t1)  &
       &      - y(iar36) * ratdum(irarap) * (1.0e0-u1)
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ar36)/d(s32)
  a1  = y(ihe4) * ratdum(irsag) &
       &    + y(ihe4) * ratdum(irsap) * (1.0e0-t1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ar36)/d(ar36)
  a1  = -y(ihe4) * ratdum(irarag) &
       &     - ratdum(irarga) &
       &     - t1 * ratdum(irargp)&
       &     - y(ihe4) * ratdum(irarap) * (1.0e0-u1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ar36)/d(ca40)
  a1  = ratdum(ircaga) &
       &    + ratdum(ircagp) * u1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  40ca jacobian elements
  !!  d(ca40)/d(he4)
  a1   = y(iar36) * ratdum(irarag)&
       &     - y(ica40) * ratdum(ircaag) &
       &     + y(iar36) * ratdum(irarap)*(1.0e0-u1)&
       &     - y(ica40) * ratdum(ircaap)*(1.0e0-v1) 
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ca40)/d(ar36)
  a1  =  y(ihe4) * ratdum(irarag)&
       &      + y(ihe4) * ratdum(irarap)*(1.0e0-u1)

  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ca40)/d(ca40)
  a1  = -y(ihe4) * ratdum(ircaag)&
       &     - ratdum(ircaga) &
       &     - ratdum(ircagp) * u1&
       &     - y(ihe4) * ratdum(ircaap)*(1.0e0-v1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ca40)/d(ti44)
  a1  = ratdum(irtiga) &
       &    + ratdum(irtigp) * v1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  44ti jacobian elements
  !!  d(ti44)/d(he4)
  a1   =  y(ica40) * ratdum(ircaag)&
       &      - y(iti44) * ratdum(irtiag)&
       &      + y(ica40) * ratdum(ircaap)*(1.0e0-v1)&
       &      - y(iti44) * ratdum(irtiap)*(1.0e0-w1)
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ti44)/d(ca40)
  a1  =  y(ihe4) * ratdum(ircaag)&
       &     + y(ihe4) * ratdum(ircaap)*(1.0e0-v1)

  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ti44)/d(ti44)
  a1  = -y(ihe4) * ratdum(irtiag)&
       &     - ratdum(irtiga) &
       &     - v1 * ratdum(irtigp)&
       &     - y(ihe4) * ratdum(irtiap)*(1.0e0-w1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ti44)/d(cr48)
  a1  = ratdum(ircrga) &
       &    + w1 * ratdum(ircrgp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  48cr jacobian elements
  !!  d(cr48)/d(he4)
  a1   =  y(iti44) * ratdum(irtiag)                                  &
       &      - y(icr48) * ratdum(ircrag)                                  &
       &      + y(iti44) * ratdum(irtiap)*(1.0e0-w1)                       &
       &      - y(icr48) * ratdum(ircrap)*(1.0e0-x1)
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(cr48)/d(ti44)
  a1  =  y(ihe4) * ratdum(irtiag)                                    &
       &      + y(ihe4) * ratdum(irtiap)*(1.0e0-w1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(cr48)/d(cr48)
  a1  = -y(ihe4) * ratdum(ircrag)                                    &
       &     - ratdum(ircrga)                                              &
       &     - w1 * ratdum(ircrgp)                                         &
       &     - y(ihe4) * ratdum(ircrap)*(1.0e0-x1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(cr48)/d(fe52)
  a1  = ratdum(irfega)                                               &
       &    + x1 * ratdum(irfegp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  52fe jacobian elements
  !!  d(fe52)/d(he4)
  a1   =  y(icr48) * ratdum(ircrag) &
       &      - y(ife52) * ratdum(irfeag) &
       &      + y(icr48) * ratdum(ircrap) * (1.0e0-x1)  &
       &      - y(ife52) * ratdum(irfeap) * (1.0e0-y1)
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(fe52)/d(cr48)
  a1  =  y(ihe4) * ratdum(ircrag) &
       &     + y(ihe4) * ratdum(ircrap) * (1.0e0-x1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(fe52)/d(fe52)
  a1  = -y(ihe4) * ratdum(irfeag)                                    &
       &     - ratdum(irfega)                                              &
       &     - x1 * ratdum(irfegp)                                         &
       &     - y(ihe4) * ratdum(irfeap) * (1.0e0-y1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(fe52)/d(ni56)
  a1  = ratdum(irniga)                                               &
       &    + y1 * ratdum(irnigp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  56ni jacobian elements
  !!  d(ni56)/d(he4)
  a1  = y(ife52) * ratdum(irfeag)                                    &
       &      + y(ife52) * ratdum(irfeap) * (1.0e0-y1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ni56)/d(fe52)
  a1  = y(ihe4) * ratdum(irfeag)                                     &
       &    + y(ihe4) * ratdum(irfeap) * (1.0e0-y1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !!  d(ni56)/d(ni56)
  a1  = -ratdum(irniga)                                              &
       &     - y1 * ratdum(irnigp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !!  bullet check the counting
  if (nt .ne. nterms) then
     write(6,*) 'nt =',nt,'  nterms =',nterms
     write(6,*) 'error in routine networkSparseJakob: nt .ne. nterms'
     call Driver_abortFlash('ERROR in networkSparseJakob: nt /= nterms')
  end if
  return

end   subroutine bn_networkSparseJakob

