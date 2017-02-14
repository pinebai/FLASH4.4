!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Iso7/bn_networkSparseJakob
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
!!  This routine is one of (a few) that can be called as an external
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

  use Burn_dataEOS, ONLY:  btemp, den
  use Driver_interface, ONLY : Driver_abortFlash
  use Burn_data
  use bn_dataIso7
  use bn_dataNetworkSize, ONLY:  eloc, nterms

!!  include 'mpif.h'

  implicit none
  !      save

#include "Flash.h"

  !
  !..this routine sets up the sparse bn_network jacobian. 
  !..input is tt (irrelevant here) and the abundances y(*). 
  !..output is the jacobian dfdy(nzo).

  !..declare arguments
  integer, intent(IN) :: nzo, nDummy
  real, intent(IN)    :: tt
  real, intent(INOUT) :: y(*)
  real, intent(OUT)   :: dfdy(*)

  !..declare locals
  integer :: i,nt,iat
  real    :: t9,yeff_ca40,yeff_ti44,                   &
       &                 t9i,t932,t9i32,rsi2ni,rni2si,rsi2nida,           &
       &                 rsi2nidsi,rni2sida,a1


  !..MPI error information
  integer   FAIL, ierr
  parameter (FAIL = -1)


  !..communicate with the jacobian builder through bn_datIso7

  !..initialize
  nt   = 0
  do i=1,nzo
     dfdy(i) = 0.0e0
  enddo


  !..positive definite mass fractions
  do i=1,NSPECIES
     y(i) = min(1.0e0,max(y(i),1.0e-30))
  enddo


  !..some common factors and branching ratios
  t9    = btemp * 1.0e-9
  t9i   = 1.0e0/t9
  t932  = t9 * sqrt(t9)
  t9i32 = 1.0e0/t932

  !..rsi2ni is the rate for silicon to nickel
  !..rni2si is the rate for nickel to silicon
  rsi2ni    = 0.0e0
  rsi2nida  = 0.0e0
  rsi2nidsi = 0.0e0
  rni2si    = 0.0e0
  rni2sida  = 0.0e0

  !      if (t9 .gt. 20.0) then
  if (t9 .gt. 2.5 .and. y(ic12)+y(io16) .le. 4.0e-3) then
     yeff_ca40 = (t9i32**3) * exp(239.42*t9i-74.741)
     yeff_ti44 = (t932**3) * exp(-274.12*t9i+74.914)
     rsi2ni    = yeff_ca40*den**3*y(ihe4)**3*ratdum(ircaag)*y(isi28)
     rsi2nida  = 3.0e0 * rsi2ni/y(ihe4)
     rsi2nidsi = rsi2ni/y(isi28)
     rni2si = min(1.0e20,yeff_ti44*ratdum(irtiga)/(den**3*y(ihe4)**3))
     rni2sida  = -3.0e0 * rni2si/y(ihe4)
     if (rni2si .eq. 1.0e20) rni2sida = 0.0e0
  end if


  !..set up the jacobian
  !..4he jacobian elements
  !..d(he4)/d(he4)
  a1 = -6.0 * y(ihe4) * y(ihe4) * ratdum(ir3a)                      &
       &    - y(ic12) * ratdum(ircag)                                     &
       &    - y(io16) * ratdum(iroag)                                     &
       &    - y(ine20) * ratdum(irneag)                                   &
       &    - y(img24) * ratdum(irmgag)                                   &
       &    - 7.0e0 * rsi2ni                                              &
       &    - 7.0e0 * rsi2nida * y(ihe4)                                  &
       &    + 7.0e0 * rni2sida * y(ini56)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(c12)
  a1 =  3.0e0 * ratdum(irg3a)                                       &
       &    - y(ihe4) * ratdum(ircag)                                     &
       &    + 2.0e0 * y(ic12) * ratdum(ir1212)                            &
       &    + 0.5e0 * y(io16) * ratdum(ir1216)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(o16)
  a1 =  ratdum(iroga)                                               &
       &    + 0.5e0 * y(ic12) * ratdum(ir1216)                            &
       &    + 2.0e0 * y(io16) * ratdum(ir1616)                            &
       &    - y(ihe4) * ratdum(iroag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(ne20)
  a1 =  ratdum(irnega)                                              &
       &    - y(ihe4) * ratdum(irneag) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(mg24)
  a1 =  ratdum(irmgga)                                              &
       &    - y(ihe4) * ratdum(irmgag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(si28)
  a1 =  ratdum(irsiga)                                              &
       &    - 7.0e0 * rsi2nidsi * y(ihe4)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(ni56)
  a1 =  7.0e0 * rni2si
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..12c jacobian elements 
  !..d(c12)/d(he4)
  a1 =  3.0e0 * y(ihe4) * y(ihe4) * ratdum(ir3a)                    &
       &    - y(ic12) * ratdum(ircag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(c12)/d(c12)
  a1 = -ratdum(irg3a)                                               &
       &    - y(ihe4) * ratdum(ircag)                                     &
       &    - 4.0e0 * y(ic12) * ratdum(ir1212)                            &
       &    - y(io16) * ratdum(ir1216)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(c12)/d(o16)
  a1 =  ratdum(iroga)                                               &
       &    - y(ic12) * ratdum(ir1216)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..16o jacobian elements 
  !..d(o16)/d(he4)
  a1 =  y(ic12) * ratdum(ircag)                                     &
       &    - y(io16) * ratdum(iroag) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(o16)/d(c12)
  a1 =  y(ihe4) * ratdum(ircag)                                     &
       &    - y(io16) * ratdum(ir1216) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(o16)/d(o16)
  a1 = -ratdum(iroga)                                               &
       &    - y(ic12) * ratdum(ir1216)                                    &
       &    - 4.0e0 * y(io16) * ratdum(ir1616)                            &
       &    - y(ihe4) * ratdum(iroag) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(o16)/d(ne20)
  a1 =  ratdum(irnega)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..20ne jacobian elements 
  !..d(ne20)/d(he4)
  a1 =  y(io16) * ratdum(iroag) - y(ine20) * ratdum(irneag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ne20)/d(c12)
  a1 =  2.0e0 * y(ic12) * ratdum(ir1212)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ne20)/d(o16)
  a1 =  y(ihe4) * ratdum(iroag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ne20)/d(ne20)
  a1 = -ratdum(irnega) - y(ihe4) * ratdum(irneag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ne20)/d(mg24)
  a1 =  ratdum(irmgga)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..24mg jacobian elements 
  !..d(mg24)/d(he4)
  a1 =  y(ine20) * ratdum(irneag)                                   &
       &    - y(img24) * ratdum(irmgag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(mg24)/d(c12)
  a1 =  0.5e0 * y(io16) * ratdum(ir1216)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(mg24)/d(o16)
  a1 =  0.5e0 * y(ic12) * ratdum(ir1216)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(mg24)/d(ne20)
  a1 =  y(ihe4) * ratdum(irneag)  
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(mg24)/d(mg24)
  a1 = -ratdum(irmgga)                                              &
       &    - y(ihe4) * ratdum(irmgag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(mg24)/d(si28)
  a1 =  ratdum(irsiga)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..28si jacobian elements 
  !..d(si28)/d(he4)
  a1 =  y(img24) * ratdum(irmgag)                                   &
       &    - rsi2ni                                                      &
       &    - rsi2nida * y(ihe4)                                          &
       &    + rni2sida * y(ini56)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(si28)/d(c12)
  a1 =  0.5e0 * y(io16) * ratdum(ir1216) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(si28)/d(o16)
  a1 =  2.0e0 * y(io16) * ratdum(ir1616)                            &
       &    + 0.5e0 * y(ic12) * ratdum(ir1216) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(si28)/d(mg24)
  a1 =  y(ihe4) * ratdum(irmgag) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(si28)/d(si28)
  a1 = -ratdum(irsiga)                                              &
       &    - rsi2nidsi * y(ihe4)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(si28)/d(ni56)
  a1 =  rni2si 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..ni56 jacobian elements 
  !..d(ni56)/d(he4)
  a1 =  rsi2ni                                                      &
       &    + rsi2nida * y(ihe4)                                          &
       &    - rni2sida * y(ini56)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ni56)/d(si28)
  a1 = rsi2nidsi * y(ihe4)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ni56)/d(ni56)
  a1 = -rni2si
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..bullet check the counting
  if (nt .ne. nterms) then
     write(6,*) 'nt =',nt,'  nterms =',nterms
     write(6,*) 'error in routine bn_networkSparseJakob: nt .ne. nterms'
     call Driver_abortFlash("Error in routine bn_networkSparseJakob: nt .ne. nterms")
     stop
  end if
  return
end subroutine bn_networkSparseJakob

