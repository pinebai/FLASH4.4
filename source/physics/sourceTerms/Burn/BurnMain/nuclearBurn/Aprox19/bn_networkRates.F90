!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19/bn_networkRates
!!
!! NAME
!!  
!!  bn_networkRates
!!
!!
!! SYNOPSIS
!! 
!!  bn_networkRates()
!!
!!  
!! DESCRIPTION
!!
!!  this routine generates nuclear reaction rates for the aprox19 network.
!!  
!!  reactions are
!!       r1616 = o16+o16      rcag  = c12(ag)     rarap = ar36(ap)
!!       r1212 = c12+c12      rclpa = cl35(pa)    r3a   = 3a-c12
!!       rg3a  = c12-3a       rkpa  = k39(pa)     rsigp = si28(gp)
!!       roga  = o16(ga)      roag  = o16(ag)     rsgp  = s32(gp)
!!       rnega = ne20(ga)     rneag = ne20(ag)    rargp = ar36(gp)
!!       rmgga = mg24(ga)     rmgag = mg24(ag)    rcagp = ca40(gp)
!!       rsiga = si28(ga)     rmgap = mg24(ap)    r1216 = o16+c12
!!       rsga  = s32(ga)      rsiag = si28(ag)    ratdum(irppg)  = p31(pg)
!!       rarga = ar36(ga)     rsiap = si28(ap)    ralpg = al27(pg)
!!       rcaga = ca40(ga)     rsag  = s32(ag)     rclpg = cl35(pg)
!!       ralpa = al27(pa)     rsap  = s32(ap)     rkpg  = k39(pg)
!!       rppa  = p31(pa)      rarag = ar36(ag)
!!       rcaag = ca40(ag)     rtiga = ti44(ga)    rscpa = sc43(pa)
!!       rcaap = ca40(ap)     rtigp = ti44(gp)    rscpg = sc43(pg)
!!       rtiag = ti44(ag)     rcrga = cr48(ga)    rvpa  = v 47(pa)
!!       rtiap = ti44(ap)     rcrgp = cr48(gp)    rvpg  = v 47(pg)
!!       rcrag = cr48(ag)     rfega = fe52(ga)    rmnpa = mn51(pa)
!!       rcrap = cr48(ap)     rfegp = fe52(gp)    rmnpg = mn51(pg)
!!       rfeag = fe52(ag)     rniga = ni56(ga)    rcopa = co55(pa)
!!       rfeap = fe52(ap)     rnigp = ni56(gp)    rcopg = co55(pg)
!!       rcogp = co55(gp)     rfepg = fe54(pg)
!!       rcpg  = c12(pg)      rnpg  = n14(pg)     ropg  = o16(pg)
!!       rnag  = n14(ag)      rpp   = p+p
!!       r34   = he3+he4      r33   = he3+he3
!!       rheng = he3(ng)      rhegn = he4(gn)     rhng  = p(ng)
!!       rdgn  = 2d(gn)       rdpg  = 2d(pg)      rhegp = 3he(gp)
!!       r52ng = fe52(ng)     r53gn = fe53(gn)    r53ng = fe53(ng)
!!       r54gn = fe54(gn)
!!  
!!  the lighter nuclei reaction rates are from cf88, some from fczii,
!!  and all the heavy nuclei rates are from sew76
!!
!!***
subroutine bn_networkRates
   
  use Burn_dataEOS
  use Burn_data
  use bn_dataAprox19

  implicit none

  !..declare  
  integer          i
  real tt9,t9r,t9,t912,t913,t923,t943,t953,t932,                    &
       &     t985,t92,t93,t972,t9r32,t9i,t9i13,t9i23,t9i32,               &
       &     t9i12,t9ri,term,term1,term2,term3,term4,rev,                 &
       &     t9ax,t9am13x,ddd,gt9cd,tq1,tq2,revx,eta,                     &
       &     t9a56x,r2abe,t9a,t9a13,t9a56,t9a23,gt9h,                     &
       &     rfepgx,tq10,rbeac,oneth,fivsix,xeta
  parameter        (oneth  = 1.0e0/3.0e0, fivsix = 5.0e0/6.0e0)


  !..zero all the rates
  real :: entropy, dst, dsd
  do i=1,nrat
     ratraw(i) = 0.0e0
  enddo

  !..limit t9 to 10, except for the reverse ratios
  tt9   = btemp * 1.0e-9
  t9r   = tt9
  t9    = min(tt9,10.0e0)
  t912  = sqrt(t9)
  t913  = t9**oneth
  t923  = t913*t913
  t943  = t9*t913
  t953  = t9*t923
  t932  = t9*t912
  t985  = t9**1.6
  t92   = t9*t9
  t93   = t9*t92
  t972  = t912*t93
  t9r32 = t9r*sqrt(t9r)
  t9i   = 1.0e0/t9
  t9i13 = 1.0e0/t913
  t9i23 = 1.0e0/t923
  t9i32 = 1.0e0/t932
  t9i12 = 1.0e0/t912
  t9ri  = 1.0e0/t9r


  !..p(p,e+nu)d(p,g)3he  
  ratraw(irpp) = 0.5e0 * 4.01e-15 * t9i23 * exp(-3.380e0*t9i13) *   &
       &               (1.0e0 + 0.123e0*t913 + 1.09e0*t923 + 0.938e0*t9) *&
       &               den


  !..3he(ng)4he  
  term           = 6.62 * (1.0e0 + 905.0*t9)
  rev            = 2.61e+10 * t9r32 * exp(-238.81*t9ri)
  ratraw(irheng) = den*term
  ratraw(irhegn) = rev*term


  !..p(ng)2d 
  term           = 4.4e+04 * (1.0e0 - 0.86*t912 + 0.429*t9)
  rev            = 4.71e+09 * t9r32 * exp(-25.82*t9ri)
  ratraw(irhng)  = den*term
  ratraw(irdgn)  = rev*term


  !..2d(pg)3he
  term           = 2.24e+03 * t9i23 * exp(-3.720*t9i13) *           &
       &                 (1.0e0 + 0.112*t913 + 3.38*t923 + 2.65*t9)
  rev            = 1.63e+10*t9r32*exp(-63.750*t9ri)
  ratraw(irdpg)  = den*term
  ratraw(irhegp) = rev*term


  !..he3(he3,2p)he4
  ratraw(ir33)   = 0.5 * den * 6.04e+10*t9i23 * exp(-12.276*t9i13) *&
       &                 (1.0e0 + 0.034*t913 - 0.522*t923 - 0.124*t9      &
       &                  + 0.353*t943 + 0.213*t953)


  !..3he(4he,g)7be; rate is beta limited (in burn) by decay of b8
  t9ax         = t9/(1.0e0 + 0.0495*t9)
  t9a56x       = t9ax**fivsix
  t9am13x      = 1.0e0/(t9ax**oneth)
  ratraw(ir34) = den*5.61e+6 * t9a56x * t9i32 * exp(-12.826*t9am13x)


  !..12c(p,g)13n ; assume c12+2p goes to n14; no beta limiting, can go to o14
  term1         = 2.04e+07 * t9i23 *                                &
       &                exp(-13.690*t9i13 - t92/2.25) *                   &
       &                (1.0e0 + 0.030*t913 + 1.19*t923 + 0.254*t9        &
       &                 + 2.06*t943 + 1.12*t953)
  term2         = 1.08e+05 * t9i32 * exp(-4.925*t9i)                &
       &                + 2.15e+05 * t9i32 * exp(-18.179*t9i)
  ratraw(ircpg) = den * (term1 + term2)


  !..14n(pg)15o; gives rate of formation of (dummy nucleus) n15. rate is limited
  !..(in burn) by the beta decay of o14 and o15
  term1          = 4.90e+07 * t9i23 *                               &
       &                 exp(-15.228*t9i13 - t92/10.850436) *             &
       &                 (1.0e0 + 0.027*t913 - 0.778*t923 - 0.149*t9      &
       &                   + 0.261*t943 + 0.127*t953)
  term2          = 2.37e+03 * t9i32 * exp(-3.011*t9i)               &
       &                 + 2.19e+04 * exp(-12.530*t9i)
  ratraw(irnpg)  = den * (term1+term2)

  !..fraction fg of n15 goes (pg) to o16 while 
  !..fraction fa of n15 goes (pa) to c12
  !..result is  n14+2p  goes to o16    rate= rnpg*fg
  !..                   goes to c12+a  rate=rnpg*fa
  !..terms 1 and 2 are for 15n(pg); terms 3 and 4 are for 15n(pa)
  term1       = 9.78e+08 * t9i23 *                                  &
       &              exp(-15.251*t9i13- t92/0.2025) *                    &
       &              (1.0e0 + 0.027*t913 + 0.219*t923 + 0.042*t9         &
       &              + 6.83*t943 + 3.32*t953)
  term2       = 1.11e+04 * t9i32 * exp(-3.328*t9i) +                &
       &              1.49e+04 * t9i32 * exp(-4.665*t9i) +                &
       &              3.80e+06 * t9i32 * exp(-11.048*t9i)
  term3       = 1.08e+12 * t9i23 *                                  &
       &              exp(-15.251*t9i13- t92/0.272484) *                  &
       &              (1.0e0 + 0.027*t913 + 2.62*t923 + 0.501*t9          &
       &              + 5.36*t943 + 2.60*t953)
  term4       = 1.19e+08 * t9i32 * exp(-3.676*t9i) +                &
       &              5.41e+08 * t9i12 * exp(-8.926*t9i) +                &
       &              0.1*(4.72e+08 * t9i32 * exp(-7.721*t9i) +           &
       &              2.20e+09 * t9i32 * exp(-11.418*t9i))
  ddd         = term1 + term2 + term3 + term4
  ratraw(ifa) = (term3 + term4)/ddd
  ratraw(ifg) = 1.0e0 - ratraw(ifa)


  !..16o(pg)17f; o16+2p goes to n14+alpha; beta limited (in burn) by f17 decay 
  ratraw(iropg) = den *                                             &
       &                1.50e+08/(t923*(1.+2.13*(1.-exp(-0.728*t923)))) * &
       &                exp(-16.692*t9i13)


  !..14n(ag)18f; n14+1.5*a goes to ne20; no beta limiting
  term1         = 7.78e+09 * t9i23 *                                &
       &                exp(-36.031*t9i13- t92/0.776161) *                &
       &                (1.0e0 + 0.012*t913 + 1.45*t923 + 0.117*t9        &
       &                 + 1.97*t943 + 0.406*t953)
  term2         = 2.36e-10 * t9i32 * exp(-2.798*t9i) +              &
       &                2.03 * t9i32 * exp(-5.054*t9i) +                  &
       &                1.15e+04 * t9i23 * exp(-12.310*t9i)
  ratraw(irnag) = den*(term1+term2)


  !..12c(ag)16o; note 1.7 times cf88 value
  ratraw(ircag) = 1.04e+08/(t92*(1.0e0 + 0.0489e0*t9i23)**2) *      &
       &                  exp(-32.120*t9i13 - t92/12.222016)              &
       &                + 1.76e+08/(t92*(1.+0.2654*t9i23)**2) *           &
       &                  exp(-32.120*t9i13)                              &
       &                + 1.25e+03*t9i32*exp(-27.499*t9i)                 &
       &                + 1.43e-02*t92*t93*exp(-15.541*t9i)
  ratraw(ircag) = ratraw(ircag) * den * 1.7e0
  rev           = 5.13e+10 * t9r32 * exp(-83.111*t9ri)
  ratraw(iroga) = rev*ratraw(ircag)/den


  !..triple alpha to c12 (has been divided by six)
  !..rate revised according to caughlan and fowler (nomoto ref.) 1988 q = -0.092
  r2abe = (7.40e+05 * t9i32)* exp(-1.0663*t9i) +                    &
       &        4.164e+09 * t9i23 * exp(-13.49*t9i13 - t92/0.009604) *    &
       &        (1.0e0 + 0.031*t913 + 8.009*t923 + 1.732*t9 +             &
       &        49.883*t943 + 27.426*t953)
  !..q = 7.367
  rbeac = (130.*t9i32) * exp(-3.3364*t9i) +                         &
       &        2.510e+07 * t9i23 * exp(-23.57*t9i13 - t92/0.055225) *    &
       &        (1.0e0 + 0.018*t913 + 5.249*t923 + 0.650*t9 +             &
       &         19.176*t943 + 6.034*t953)
  !..q = 7.275
  if (t9.gt.0.08) then
     ratraw(ir3a) = 2.90e-16 * (r2abe*rbeac) +                      &
          &                  0.1 * 1.35e-07 * t9i32 * exp(-24.811*t9i)
  else
     ratraw(ir3a) = 2.90e-16*(r2abe*rbeac) *                        &
          &         (0.01 + 0.2*(1.0e0 + 4.0e0*exp(-(0.025*t9i)**3.263)) /   &
          &         (1.0e0 + 4.0e0*exp(-(t9/0.025)**9.227))) +               &
          &         0.1 * 1.35e-07 * t9i32 * exp(-24.811*t9i)
  end if
  ratraw(ir3a)  = ratraw(ir3a) * (den*den)/6.0e0
  rev           = 2.00e+20*exp(-84.424*t9ri)
  ratraw(irg3a) = rev*(t9r*t9r*t9r) * 6.0e0 *ratraw(ir3a)/(den*den)



  !..c12 + c12 reaction; see cf88 references 47-49
  t9a            = t9/(1.0e0 + 0.0396*t9)
  t9a13          = t9a**oneth
  t9a56          = t9a**fivsix
  ratraw(ir1212) = 4.27e+26 * t9a56/t932 *                          &
       &                 exp(-84.165/t9a13-2.12e-03*t9*t9*t9)
  ratraw(ir1212) = 0.5*den*ratraw(ir1212)


  !..c12 + o16 reaction;  q = 16.755; valid for t9 .gt. 0.5
  !..y(nc12)*y(no16)*rc28 is the rate of formation of the si28 compound nucleus
  if (t9.ge.0.5) then
     t9a   = t9/(1.+0.055*t9)
     t9a13 = t9a**oneth
     t9a23 = t9a13*t9a13
     t9a56 = t9a**fivsix
     ratraw(ir1216) = 1.72e+31 * t9a56 * t9i32 * exp(-106.594/t9a13) /&
          &                  (exp(-0.18*t9a*t9a) + 1.06e-03*exp(2.562*t9a23))
     ratraw(ir1216) = den*ratraw(ir1216)
  else
     ratraw(ir1216) = 0.0
  endif


  !..16o+16o rate; q = 16.542; 
  !..y16*y16*r32 is rate of formation of 32s compound nucleus
  ratraw(ir1616) = 7.10e+36 * t9i23 *                               &
       &                 exp(-135.93 * t9i13 - 0.629*t923 -               &
       &                     0.445*t943 + 0.0103*t9*t9)
  ratraw(ir1616) = den*ratraw(ir1616)*0.5


  !..16o(ag)20ne + inverse
  term1          = 9.37e+09 * t9i23 *                               &
       &                 exp(-39.757*t9i13- t92/2.515396)
  term2 =          62.1 * t9i32 * exp(-10.297*t9i) +                &
       &                 538.0e0 * t9i32 * exp(-12.226*t9i) +             &
       &                 13.0e0 * t92 * exp(-20.093*t9i)
  ratraw(iroag)  = den*(term1+term2)
  rev            = 5.65e+10*t9r32*exp(-54.937*t9ri)
  ratraw(irnega) = rev*(term1+term2)


  !..20ne(ag)24mg + inverse
  term1          = 4.11e+11 * t9i23 *                               &
       &                 exp(-46.766*t9i13- t92/4.923961) *               &
       &                 (1.0e0 + 0.009*t913 + 0.882*t923 + 0.055*t9 +    &
       &                 0.749*t943 + 0.119*t953)
  term2          = 5.27e+03 * t9i32 * exp(-15.869*t9i) +            &
       &                 6.51e+03 * t912 * exp(-16.223*t9i)
  term3          = 0.1e0 * (42.1 * t9i32 * exp(-9.115*t9i) +        &
       &                 32.0 * t9i23 * exp(-9.383*t9i))
  ratraw(irneag) = den*(term1+term2+term3)/                         &
       &                 (1.0e0 + 5.0e0*exp(-18.960*t9i))
  rev            = 6.01e+10 * t9r32 * exp(-108.059*t9ri)
  ratraw(irmgga) = (ratraw(irneag)/den)*rev


  !..24mg(ag)28si + inverse
  gt9cd          = (1.0e0 + 5.0e0*exp(-15.882*t9i))
  ratraw(irmgag) = den * (4.78e+01 * t9i32 * exp(-13.506*t9i) +     &
       &                 2.38e+03 * t9i32 * exp(-15.218*t9i) +            &
       &                 2.47e+02 * t932 * exp(-15.147*t9i) +             &
       &                 0.1 * (1.72e-09 * t9i32 * exp(-5.028*t9i) +      &
       &                 1.25e-03 * t9i32 * exp(-7.929*t9i) +             &
       &                 2.43e+01 * t9i * exp(-11.523*t9i)))/gt9cd
  rev            = 6.27e+10 * t9r32 * exp(-115.862*t9ri)
  ratraw(irsiga) = rev * ratraw(irmgag)/den


  !..24mg(ap)27al + inverse
  gt9h           = 1.0e0 + exp(-9.792*t9i)/3.0e0 +                  &
       &                 2.0e0 * exp(-11.773*t9i)/3.0e0
  term1          = 1.10e+08 * t9i23 *                               &
       &                 exp(-23.261*t9i13 - t92/0.024649) *              &
       &                 (1.0e0 + 0.018*t913 + 12.85*t923 + 1.61*t9 +     &
       &                  89.87*t943 + 28.66*t953)
  term2          = 129.0e0 * t9i32 * exp(-2.517*t9i) +              &
       &                 5660.0e0 * t972 * exp(-3.421*t9i) +              &
       &                 0.1 * (3.89e-08 * t9i32 * exp(-0.853*t9i) +      &
       &                 8.18e-09 * t9i32 * exp(-1.001*t9i))
  rev            = 1.81 * exp(-18.572*t9ri)
  ratraw(irmgap) = rev * den * (term1+term2)/gt9h
  ratraw(iralpa) = den * (term1+term2)/gt9h


  !..27al(pg)28si
  ratraw(iralpg) = den * (1.67e+08 *t9i23 *                         &
       &                 exp(-23.261*t9i13- t92/0.024025) *               &
       &                 (1.0e0 + 0.018*t913 + 5.81*t923 + 0.728*t9       &
       &                  + 27.31*t943 + 8.71*t953) +                     &
       &                 2.20e+00 * t9i32 * exp(-2.269*t9i) +             &
       &                 1.22e+01 * t9i32 * exp(-2.491*t9i) +             &
       &                 1.50e+04 * t9 * exp(-4.112*t9i) +                &
       &                 0.1 * (6.50e-10 * t9i32 * exp(-0.853*t9i) +      &
       &                 1.63e-10 * t9i32 * exp(-1.001*t9i)))/gt9h
  rev            = 1.13e+11*t9r32*exp(-134.434*t9ri)
  ratraw(irsigp) = ratraw(iralpg)*rev/den


  !..28si(ag)32s + inverse
  term           = 4.82e+22 * t9i23 * exp(-61.015 * t9i13 *         &
       &                (1.0e0 + 6.340e-02*t9 + 2.541e-03*t92 -           &
       &                 2.900e-04*t93))
  rev            = 6.461e+10 * t9r32 * exp(-80.643*t9ri)
  ratraw(irsiag) = den * term
  ratraw(irsga)  = rev * term


  !..28si(ap)31p + inverse; (given as 31p(pa)28si)
  term           = 4.16e+13 * t9i23 * exp(-25.631 * t9i13 *         &
       &                 (1.0e0 + 2.798e-03*t9 + 2.763e-03*t92 -          &
       &                 2.341e-04*t93))
  rev            = 0.5825e0 * exp(-22.224*t9ri)
  ratraw(irppa)  = den*term
  ratraw(irsiap) = ratraw(irppa)*rev


  !..31p(pg)32s + inverse
  term          = 1.08e+16 * t9i23 * exp(-27.042 * t9i13 *          &
       &                (1.0e0 + 1.928e-01*t9 - 1.540e-02*t92 +           &
       &                6.444e-04*t93))
  rev           = 3.764e+10 * t9r32 * exp(-102.865*t9ri)
  ratraw(irppg) = den*term
  ratraw(irsgp) = rev*term


  !..32s(ag)36ar + inverse
  term           = 1.16e+24 * t9i23 * exp(-66.690 * t9i13 *         &
       &                 (1.0e0 + 4.913e-02*t9 + 4.637e-03*t92 -          &
       &                 4.067e-04*t93))
  rev            = 6.616e+10 * t9r32 * exp(-77.080*t9ri)
  ratraw(irsag)  = den*term
  ratraw(irarga) = rev*term


  !..32s(ap)35cl + inverse (given as 35cl(pa)32s)
  term           = 1.27e+16 * t9i23 * exp(-31.044 * t9i13 *         &
       &                 (1.0e0 + 1.041e-01*t9 - 1.368e-02*t92 +          &
       &                 6.969e-04*t93))
  rev            = 1.144 * exp(-21.643*t9ri)
  ratraw(irsap)  = den*rev*term
  ratraw(irclpa) = den*term


  !..35cl(pg)36ar
  term           = 4.48e+16 * t9i23 * exp(-29.483 * t9i13 *         &
       &                 (1.0e0 + 1.761e-01*t9 - 1.322e-02*t92 +          &
       &                 5.245e-04*t93))
  rev            = 7.568e+10*t9r32*exp(-98.722*t9ri)
  ratraw(irclpg) = den*term
  ratraw(irargp) = rev*term


  !..36ar(ag)40ca + inverse
  term           = 2.81e+30 * t9i23 * exp(-78.271 * t9i13 *         &
       &                (1.0e0 + 1.458e-01*t9 - 1.069e-02*t92 +           &
       &                 3.790e-04*t93))
  rev            = 6.740e+10 * t9r32 * exp(-81.711*t9ri)
  ratraw(irarag) = den*term
  ratraw(ircaga) = rev*term


  !..36ar(ap)39k + inverse (given as 39k(pa)36ar)
  term           = 2.76e+13 * t9i23 * exp(-34.922 * t9i13 *         &
       &                 (1.0e0 + 4.826e-03*t9 - 5.534e-03*t92 +          &
       &                 4.021e-04*t93))
  rev            = 1.128*exp(-14.959*t9ri)
  ratraw(irarap) = den*term*rev
  ratraw(irkpa)  = den*term


  !..39k(pg)40ca + inverse 
  term           = 4.09e+16 * t9i23 * exp(-31.727 * t9i13 *         &
       &                 (1.0e0 + 1.622e-01*t9 - 1.119e-02*t92 +          &
       &                 3.910e-04*t93))
  rev            = 7.600e+10 * t9r32 * exp(-96.657*t9ri)
  ratraw(irkpg)  = den*term
  ratraw(ircagp) = rev*term


  !..40ca(ag)44ti + inverse 
  term           = 4.66e+24 * t9i23 * exp(-76.435 * t9i13 *         &
       &                 (1.0e0 + 1.650e-02*t9 + 5.973e-03*t92 -          &
       &                 3.889e-04*t93))
  rev            = 6.843e+10 * t9r32 * exp(-59.510*t9ri)
  ratraw(ircaag) = den*term
  ratraw(irtiga) = rev*term


  !..40ca(ap)43sc + inverse (given as 43sc(pa)40ca)
  term           = 4.54e+14 * t9i23 * exp(-32.177 * t9i13 *         &
       &                 (1.0e0 - 1.206e-02*t9 + 7.753e-03*t92 -          &
       &                 5.071e-04*t93))
  rev            = 2.229 * exp(-40.966*t9ri)
  ratraw(ircaap)= den*rev*term
  ratraw(irscpa) = den*term


  !..43sc(pg)44ti + inverse 
  term           = 3.85e+16 * t9i23 * exp(-33.234 * t9i13 *         &
       &                 (1.0e0 + 1.023e-01*t9 - 2.242e-03*t92 -          &
       &                 5.463e-05*t93))
  rev            = 1.525e+11 * t9r32 * exp(-100.475*t9ri)
  ratraw(irscpg) = den*term
  ratraw(irtigp) = rev*term


  !..44ti(ag)48cr + inverse 
  term           = 1.37e+26 * t9i23 * exp(-81.227 * t9i13 *         &
       &                 (1.0e0 + 1.066e-01*t9 - 1.102e-02*t92 +          &
       &                 5.324e-04*t93))
  rev            = 6.928e+10*t9r32*exp(-89.289*t9ri)
  ratraw(irtiag) = den*term
  ratraw(ircrga) = rev*term


  !..44ti(ap)47v + inverse (given as 47v(pa)44ti)
  term           = 6.54e+20 * t9i23 * exp(-66.678 * t9i13 *         &
       &                 (1.0e0 + 2.655e-02*t9 - 3.947e-03*t92 +          &
       &                 2.522e-04*t93))
  rev            = 1.104 * exp(-4.723*t9ri)
  ratraw(irtiap) = rev*den*term
  ratraw(irvpa)  = den*term


  !..47v(pg)48cr + inverse
  term           = 2.05e+17 * t9i23 * exp(-35.568 * t9i13 *         &
       &                 (1.0e0 + 9.979e-02*t9 - 2.269e-03*t92 -          &
       &                 6.662e-05*t93))
  rev            = 7.649e+10*t9r32*exp(-93.999*t9ri)
  ratraw(irvpg)  = den*term
  ratraw(ircrgp) = rev*term


  !..48cr(ag)52fe + inverse
  term           = 1.04e+23 * t9i23 * exp(-81.420 * t9i13 *         &
       &                 (1.0e0 + 6.325e-02*t9 - 5.671e-03*t92 +          &
       &                 2.848e-04*t93))
  rev            = 7.001e+10 * t9r32 * exp(-92.177*t9ri)
  ratraw(ircrag) = den*term
  ratraw(irfega) = rev*term


  !..48cr(ap)51mn + inverse
  term           = 1.83e+26 * t9i23 * exp(-86.741 * t9i13 *         &
       &                 (1.0e0 + 1.384e-02*t9 + 1.081e-03*t92 -          &
       &                 5.933e-05*t93))
  rev            = 0.6087*exp(-6.510*t9ri)
  ratraw(ircrap) = den*term
  ratraw(irmnpa) = rev*den*term


  !..51mn(pg)52fe + inverse 
  term           = 3.77e+17 * t9i23 * exp(-37.516 * t9i13 *         &
       &                 (1.0e0 + 8.922e-02*t9 - 1.256e-03*t92 -          &
       &                  9.453e-05*t93))
  rev            = 1.150e+11*t9r32*exp(-85.667*t9ri)
  ratraw(irmnpg) = den*term
  ratraw(irfegp) = rev*term


  !..52fe(ag)56ni + inverse 
  term           = 1.05e+27 * t9i23 * exp(-91.674 * t9i13 *         &
       &                 (1.0e0 + 7.846e-02*t9 - 7.430e-03*t92 +          &
       &                 3.723e-04*t93))
  rev            = 7.064e+10*t9r32*exp(-92.850*t9ri)
  ratraw(irfeag) = den*term
  ratraw(irniga) = rev*term


  !..52fe(ap)55co + inverse
  term           = 1.30e+27 * t9i23 * exp(-91.674 * t9i13 *         &
       &                 (1.0e0 + 1.367e-02*t9 + 7.428e-04*t92 -          &
       &                 3.050e-05*t93))
  rev            = 0.4597*exp(-9.470*t9ri)
  ratraw(irfeap) = den*term
  ratraw(ircopa) = rev*den*term


  !..55co(pg)56ni + inverse 
  term           = 1.21e+18 * t9i23 * exp(-39.604 * t9i13 *         &
       &                 (1.0e0 + 9.894e-02*t9 - 3.131e-03*t92 -          &
       &                 2.160e-05*t93))
  rev            = 1.537e+11*t9r32*exp(-83.382*t9ri)
  ratraw(ircopg) = den*term
  ratraw(irnigp) = rev*term


  !..54fe(pg)55co + inverse 
  term           = 4.51e+17 * t9i23 * exp(-38.483 * t9i13 *         &
       &                 (1.0e0 + 9.593e-02*t9 - 3.445e-03*t92 +          &
       &                 8.594e-05*t93))
  rev            = 2.400e+09 * t9r32 * exp(-58.605*t9ri)
  ratraw(irfepg) = den*term
  ratraw(ircogp) = rev*term
  revx  = rev
  rfepgx= ratraw(irfepg)


  !..52fe(ng)53fe + inverse
  tq1            = t9/0.348
  tq2            = t9-0.348e0
  tq10           = tq1**0.10
  term           = 9.604e+05 * exp(-0.0626*tq2)
  ratraw(ir52ng) = den*term
  ratraw(ir53gn) = 2.43e+09 * t9r32 * exp(-123.951*t9ri) * term


  !..53fe(ng)54fe + inverse 
  term           = 1.817e+06 * tq10 * exp(-0.06319*tq2)
  ratraw(ir53ng) = den*term
  ratraw(ir54gn) = 1.56e+11 * t9r32 * term * exp(-155.284*t9ri)


  !..no electron capture rates for now
  ratraw(irpen)   = 0.0e0
  ratraw(irnep)   = 0.0e0
  ratraw(ispen)   = 0.0e0
  ratraw(isnep)   = 0.0e0
  ratraw(irn56ec) = 0.0e0
  ratraw(isn56ec) = 0.0e0

  return
end subroutine bn_networkRates





