!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bn_sneutx
!!
!! NAME
!!
!!  bn_sneutx
!!
!! SYNOPSIS
!!
!!  call bn_sneutx()
!!
!! ARGUMENTS
!!
!! DESCRIPTION
!!
!!    this routine computes neutrino losses from the analytic fits of
!!    itoh et al. apjs 102, 411, 1996
!!    
!!    input is the temperature temp, density den, 
!!    mean number of nucleons abar, and mean charge zbar.
!!    the input comes from data structure Burn_dataEOS
!!    
!!    output into Burn_data is the 
!!    pair neutrino contributions spair,
!!    plasma neutrino contributions splas, 
!!    photoneutrino contributions sphot, 
!!    bremstrahlung neutrino contributions sbrem, 
!!    recombination neutrino contributions srecomb
!!    and the sum total ?? sneut,
!!    all in erg/g/s.
!!
!!  NOTES
!!    Called by bn_ifermi12
!!
!!***

subroutine bn_sneutx

  use Burn_dataEOS, ONLY:  btemp,den,abar,zbar,z2bar,ytot1,ye
  use Burn_data, ONLY: sneut,sphot,spair,splas,sbrem,srecomb
  use bn_interface, ONLY: bn_ifermi12

  implicit none
  
!!    
!!    declare local variables
!!  real :: bn_ifermi12    ! external function not needed after interfaces added
  real :: xmue,t9,xl,xlp5,xl2,xl3,xl4,xl5,xl6,xl8,xl9, & 
       &                 xlm1,xlm2,xlm3,rm,xnum,xden,fpair,fphoto,gl, & 
       &                 zeta,zeta2,zeta3,qpair,gl2,gl12,gl32,gl72,gl6, & 
       &                 ft,fl,cc,xlnt,fxy,qv,tau,qphoto,den6,tfermi,t8, & 
       &                 t832,t83,t86,t8m2,t8m5,eta,etam1,etam2,fbrem, & 
       &                 gbrem,a0,a1,a2,a3,c00, & 
       &                 c01,c02,c03,c04,c05,c06,c10,c11,c12,c13,c14,c15, & 
       &                 c16,c20,c21,c22,c23,c24,c25,c26,dd01,dd02,dd03, & 
       &                 dd04,dd05,dd11,dd12,dd13,dd14,dd15,dd21,dd22, & 
       &                 dd23,dd24,dd25,u,gamma,gm1,gm13,gm23,v,w,fb, & 
       &                 gb,gt,fliq,gliq,nu,nu2,nu3,b,c,d,f1, & 
       &                 f2,f3,z,bigj,cos1,cos2,cos3,cos4,cos5,sin1,sin2, & 
       &                 sin3,sin4,sin5,last,deni

  real, parameter   ::    pi      = 3.1415926535897932384e0, & 
       &                  fac1    = 5.0e0 * pi / 3.0e0, & 
       &                  fac2    = 10.0e0 * pi, & 
       &                  fac3    = pi / 5.0e0, & 
       &                  third   = 1.0e0/3.0e0, & 
       &                  twoth   = 2.0e0/3.0e0, & 
       &                  con1    = 1.0e0/5.9302e0, & 
       &                  sixth   = 1.0e0/6.0e0


!!    theta is sin**2(theta_weinberg) = 0.2319 plus/minus 0.00005 (1996)
!!    xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
!!    change theta and xnufam if need be, and the changes will automatically
!!    propagate through the routine. cv and ca are the vector and axial currents.

  real, parameter  ::     theta  = 0.2319e0, & 
       &                  xnufam = 3.0e0, & 
       &                  cv     = 0.5e0 + 2.0e0 * theta, & 
       &                  cvp    = 1.0e0 - cv, & 
       &                  ca     = 0.5e0, & 
       &                  cap    = 1.0e0 - ca, & 
       &                  tfac1  = cv*cv + ca*ca +  & 
       &                           (xnufam-1.0e0) * (cvp*cvp+cap*cap), & 
       &                  tfac2  = cv*cv - ca*ca +  & 
       &                           (xnufam-1.0e0) * (cvp*cvp - cap-cap), & 
       &                  tfac3  = tfac2/tfac1, & 
       &                  tfac4  = 0.5e0 * tfac1, & 
       &                  tfac5  = 0.5e0 * tfac2, & 
       &                  tfac6  = cv*cv + 1.5e0*ca*ca + (xnufam - 1.0e0)* & 
       &                           (cvp*cvp + 1.5e0*cap*cap)


!!    initialize and bail if its too cold
  spair    = 0.0e0
  splas    = 0.0e0
  sphot    = 0.0e0
  sbrem    = 0.0e0
  srecomb  = 0.0e0
  sneut    = 0.0e0


!!    bail from the rest if its too cold
  if (btemp .lt. 1.0e7) return


!!    some frequent factors
  xmue  = abar/zbar
  ye    = zbar/abar
  t9    = btemp * 1.0e-9

  xl    = t9 * con1
  xlp5  = sqrt(xl)
  xl2   = xl*xl
  xl3   = xl2*xl
  xl4   = xl2*xl2
  xl5   = xl*xl4
  xl6   = xl2*xl4
  xl8   = xl4*xl4
  xl9   = xl8 * xl
  xlm1  = 1.0e0/xl
  xlm2  = xlm1*xlm1
  xlm3  = xlm1*xlm2

  rm    = den/xmue
  zeta  = (rm * 1.0e-9)**third * xlm1
  zeta2 = zeta * zeta
  zeta3 = zeta2 * zeta


!!    pair neutrino section
!!    for reactions like e+ + e- => nu_e + nubar_e 
!!    equation 2.8 
  gl = 1.0e0 - 13.04e0*xl2 + 133.5e0*xl4 +1534.0e0*xl6 + 918.6e0*xl8

!!    equation 2.7
  if (t9 .lt. 10.0) then
     xnum = (6.002e19 + 2.084e20*zeta + 1.872e21*zeta2) & 
          &         * exp(-5.5924e0*zeta)
     xden = zeta3 + 9.383e-1*xlm1 - 4.141e-1*xlm2 + 5.829e-2*xlm3
  else
     xnum = (6.002e19 + 2.084e20*zeta + 1.872e21*zeta2) & 
          &         * exp(-4.9924e0*zeta)
     xden = zeta3 + 1.2383e0*xlm1 - 8.141e-1*xlm2 
  end if
  fpair = xnum/xden

!!    equation 2.6
  xnum   = 1.0e0 / (10.7480e0*xl2 + 0.3967e0*xlp5 + 1.005e0)
  xden   = (1.0e0 + rm/(7.692e7*xl3 + 9.715e6*xlp5) )**(-0.3e0)
  qpair  = xnum*xden

!!    equation 2.5
  spair = tfac4*(1.0e0 + tfac3 * qpair)*gl*exp(-2.0e0*xlm1)*fpair




!!    plasma neutrino section 
!!    for collective reactions like gamma_plasmon => nu_e + nubar_e
!!    equation 4.6
  gl2  = 1.1095e11*rm/(btemp*btemp *  & 
       &       sqrt(1.0e0+(1.019e-6*rm)**twoth))
  gl   = sqrt(gl2)
  gl12 = sqrt(gl)
  gl32 = gl * gl12
  gl72 = gl2 * gl32
  gl6  = gl2 * gl2 * gl2

!!    equation 4.7
  ft   = 2.4e0 + 0.6e0*gl12 + 0.51e0*gl + 1.25e0*gl32

!!    equation 4.8
  fl   = (8.6e0*gl2 + 1.35e0*gl72)/(225.0e0 - 17.0e0*gl + gl2)

!!    equation 4.9 and 4.10
  cc   = log10(2.0e0*rm)
  xlnt = log10(btemp)
  xnum = sixth * ( 17.5e0 + cc - 3.0e0 * xlnt)
  xden = sixth * (-24.5e0 + cc + 3.0e0 * xlnt)

!!    equation 4.11
  if (abs(xnum) .gt. 0.7e0  .or.  xden .lt. 0.0e0) then
     fxy = 1.0e0
  else 
     fxy = 1.05e0 + (0.39e0 - 1.25e0*xnum - 0.35e0*sin(4.5e0*xnum) & 
          &                 - 0.3e0 * exp(-1.0e0*(4.5e0*xnum + 0.9e0)**2)) & 
          &       * exp(-1.0e0* ( min(0.0e0, xden - 1.6e0 + 1.25e0*xnum) & 
          &                 / (0.57e0 - 0.25e0*xnum) )**2  )
  end if

!!    equation 4.5
  qv = 3.0e21 * xl9 * gl6 * exp(-gl) * (ft + fl) * fxy

!!    equation 4.1
  splas = 0.93153e0 * qv




!!    photoneutrino process section  
!!    for reactions like e- + gamma => e- + nu_e + nubar_e
!!                       e+ + gamma => e+ + nu_e + nubar_e
!!    equation 3.8 for tau, equation 3.6 for cc,
!!    and table 2 written out for speed
  if (btemp .ge. 1.0e7  .and. btemp .lt. 1.0e8) then
     tau  =  log10(btemp * 1.0e-7)
     cc   =  0.5654e0 + tau
     c00  =  1.008e11
     c01  =  0.0e0
     c02  =  0.0e0
     c03  =  0.0e0
     c04  =  0.0e0
     c05  =  0.0e0
     c06  =  0.0e0
     c10  =  8.156e10
     c11  =  9.728e8
     c12  = -3.806e9
     c13  = -4.384e9
     c14  = -5.774e9
     c15  = -5.249e9
     c16  = -5.153e9
     c20  =  1.067e11
     c21  = -9.782e9 
     c22  = -7.193e9
     c23  = -6.936e9
     c24  = -6.893e9
     c25  = -7.041e9
     c26  = -7.193e9
     dd01 =  0.0e0
     dd02 =  0.0e0
     dd03 =  0.0e0
     dd04 =  0.0e0
     dd05 =  0.0e0
     dd11 = -1.879e10
     dd12 = -9.667e9
     dd13 = -5.602e9
     dd14 = -3.370e9
     dd15 = -1.825e9
     dd21 = -2.919e10
     dd22 = -1.185e10
     dd23 = -7.270e9
     dd24 = -4.222e9
     dd25 = -1.560e9

  else if (btemp .ge. 1.0e8  .and. btemp .lt. 1.0e9) then
     tau  =  log10(btemp * 1.0e-8)
     cc   =  1.5654e0
     c00  =  9.889e10 
     c01  = -4.524e8
     c02  = -6.088e6 
     c03  =  4.269e7 
     c04  =  5.172e7 
     c05  =  4.910e7 
     c06  =  4.388e7
     c10  =  1.813e11
     c11  = -7.556e9 
     c12  = -3.304e9  
     c13  = -1.031e9
     c14  = -1.764e9  
     c15  = -1.851e9
     c16  = -1.928e9
     c20  =  9.750e10
     c21  =  3.484e10
     c22  =  5.199e9  
     c23  = -1.695e9  
     c24  = -2.865e9  
     c25  = -3.395e9  
     c26  = -3.418e9
     dd01 = -1.135e8   
     dd02 =  1.256e8   
     dd03 =  5.149e7   
     dd04 =  3.436e7   
     dd05 =  1.005e7
     dd11 =  1.652e9  
     dd12 = -3.119e9  
     dd13 = -1.839e9  
     dd14 = -1.458e9  
     dd15 = -8.956e8
     dd21 = -1.549e10  
     dd22 = -9.338e9  
     dd23 = -5.899e9  
     dd24 = -3.035e9  
     dd25 = -1.598e9

  else if (btemp .ge. 1.0e9) then
     tau  =  log10(t9)
     cc   =  1.5654e0
     c00  =  9.581e10
     c01  =  4.107e8
     c02  =  2.305e8   
     c03  =  2.236e8   
     c04  =  1.580e8   
     c05  =  2.165e8   
     c06  =  1.721e8
     c10  =  1.459e12
     c11  =  1.314e11
     c12  = -1.169e11  
     c13  = -1.765e11  
     c14  = -1.867e11  
     c15  = -1.983e11  
     c16  = -1.896e11
     c20  =  2.424e11
     c21  = -3.669e9
     c22  = -8.691e9  
     c23  = -7.967e9  
     c24  = -7.932e9  
     c25  = -7.987e9  
     c26  = -8.333e9
     dd01 =  4.724e8
     dd02 =  2.976e8   
     dd03 =  2.242e8   
     dd04 =  7.937e7   
     dd05 =  4.859e7
     dd11 = -7.094e11
     dd12 = -3.697e11
     dd13 = -2.189e11  
     dd14 = -1.273e11  
     dd15 = -5.705e10
     dd21 = -2.254e10
     dd22 = -1.551e10
     dd23 = -7.793e9
     dd24 = -4.489e9
     dd25 = -2.185e9
  end if

!!    equation 3.7, compute the expensive trig functions only one time
  cos1 = cos(fac1*tau)
  cos2 = cos(fac1*2.0e0*tau)
  cos3 = cos(fac1*3.0e0*tau)
  cos4 = cos(fac1*4.0e0*tau)
  cos5 = cos(fac1*5.0e0*tau)
  last = cos(fac2*tau)

  sin1 = sin(fac1*tau)
  sin2 = sin(fac1*2.0e0*tau)
  sin3 = sin(fac1*3.0e0*tau)
  sin4 = sin(fac1*4.0e0*tau)
  sin5 = sin(fac1*5.0e0*tau)

  a0 = 0.5e0*c00  & 
       &     + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2 & 
       &     + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4 & 
       &     + c05*cos5 + dd05*sin5 + 0.5e0*c06*last
  a1 = 0.5e0*c10  & 
       &     + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2 & 
       &     + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4 & 
       &     + c15*cos5 + dd15*sin5 + 0.5e0*c16*last
  a2 = 0.5e0*c20  & 
       &     + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2 & 
       &     + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4 & 
       &     + c25*cos5 + dd25*sin5 + 0.5e0*c26*last

!!    equation 3.4
  xnum   = (a0 + a1*zeta + a2*zeta2)*exp(-cc*zeta)
  xden   = zeta3 + 6.290e-3*xlm1 + 7.483e-3*xlm2 + 3.061e-4*xlm3
  fphoto = xnum/xden

!!    equation 3.3
  xnum   = 0.666e0*((1.0e0 + 2.045e0 * xl)**(-2.066e0))
  xden   = 1.0e0 + rm / (1.875e8*xl + 1.653e8*xl2  & 
       &                     + 8.449e8*xl3 - 1.604e8*xl4)
  qphoto = xnum/xden

!!    equation 3.2
  sphot   =  tfac4 * (1.0e0 - tfac3 * qphoto) * rm * xl5 * fphoto





!!    bremsstrahlung neutrino section 
!!    for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
!!                       n  + n     => n + n + nu + nubar
!!                       n  + p     => n + p + nu + nubar
!!    equation 4.3
  den6   = den * 1.0e-6
  tfermi = 5.9302e9*(sqrt(1.0e0+1.018e0*(den6/xmue)**twoth)-1.0e0)

!!    "weak" degenerate electrons only
  if (btemp .gt. 0.3e0 * tfermi) then
     t8     = btemp * 1.0e-8
     t832   = t8 * sqrt(t8)
     t83    = t8 * t8 * t8
     t86    = t83 * t83
     t8m2   = 1.0e0/(t8 * t8)
     t8m5   = t8m2 * t8m2/t8

   !!    equation 5.3
     eta   = rm/(7.05e6 * t832 + 5.12e4 * t83)
     etam1 = 1.0e0/eta
     etam2 = etam1 * etam1

   !!    equation 5.2
     xnum  = 1.0e0/(23.5e0 + 6.83e4*t8m2 + 7.81e8*t8m5)
     xden  = 1.26e0*(1.0e0+etam1)/(1.0e0+1.47e0*etam1+3.29e-2*etam2)
     fbrem = xnum + xden

   !!    equation 5.9
     xnum  = 1.0e0/( (1.0e0 + rm*1e-9) & 
          &                * (230.0e0 + 6.7e5*t8m2 + 7666e9*t8m5) )
     c00   = 7.75e5*t832 + 247.0e0*t8**(3.85e0)
     c01   = 4.07e0 + 0.0240e0 * t8**(1.4e0)
     c02   = 4.59e-5 * t8**(-0.110e0)
     xden  = 1.0e0/(c00/rm + c01 + c02 * den**(0.656e0))
     gbrem = xnum + xden

   !!    equation 5.1
     sbrem  = 0.5738e0*zbar*ye*t86*den * (tfac4*fbrem - tfac5*gbrem)


   !!    liquid metal with c12 parameters (not too different for other elements)
   !!    equation 5.18 and 5.16
  else
     t8    = btemp * 1.0e-8
     t86   = t8 * t8 * t8 * t8 * t8 * t8
     u     = fac3 * (log10(den) - 3.0e0)
     gamma = 2.275e-1 * zbar * zbar/t8 * (den6/abar)**third
     gm1   = 1.0e0/gamma
     gm13  = gm1**third
     gm23  = gm13 * gm13

   !!    equation 5.25 and 5.26
     v = -0.05483e0 - 0.01946e0*gm13 + 1.86310e0*gm23 - 0.78873e0*gm1
     w = -0.06711e0 + 0.06859e0*gm13 + 1.74360e0*gm23 - 0.74498e0*gm1

   !!    compute the expensive trig functions of equation 5.21 only once
     cos1 = cos(u)
     cos2 = cos(2.0e0*u)
     cos3 = cos(3.0e0*u)
     cos4 = cos(4.0e0*u)
     cos5 = cos(5.0e0*u)

     sin1 = sin(u)
     sin2 = sin(2.0e0*u)
     sin3 = sin(3.0e0*u)
     sin4 = sin(4.0e0*u)

   !!    equation 5.21
     fb =  0.5e0 * 0.17946e0  + 0.00945e0*u + 0.34529e0    & 
          &       - 0.05821e0*cos1 - 0.04969e0*sin1 & 
          &       - 0.01089e0*cos2 - 0.01584e0*sin2 & 
          &       - 0.01147e0*cos3 - 0.00504e0*sin3 & 
          &       - 0.00656e0*cos4 - 0.00281e0*sin4 & 
          &       - 0.00519e0*cos5 

   !!    equation 5.22
     ft =  0.5e0 * 0.06781e0 - 0.02342e0*u + 0.24819e0 & 
          &       - 0.00944e0*cos1 - 0.02213e0*sin1 & 
          &       - 0.01289e0*cos2 - 0.01136e0*sin2 & 
          &       - 0.00589e0*cos3 - 0.00467e0*sin3 & 
          &       - 0.00404e0*cos4 - 0.00131e0*sin4 & 
          &       - 0.00330e0*cos5 

   !!    equation 5.23
     gb =  0.5e0 * 0.00766e0 - 0.01259e0*u + 0.07917e0 & 
          &       - 0.00710e0*cos1 + 0.02300e0*sin1 & 
          &       - 0.00028e0*cos2 - 0.01078e0*sin2 & 
          &       + 0.00232e0*cos3 + 0.00118e0*sin3 & 
          &       + 0.00044e0*cos4 - 0.00089e0*sin4 & 
          &       + 0.00158e0*cos5

   !!    equation 5.24
     gt =  -0.5e0 * 0.00769e0  - 0.00829e0*u + 0.05211e0 & 
          &       + 0.00356e0*cos1 + 0.01052e0*sin1 & 
          &       - 0.00184e0*cos2 - 0.00354e0*sin2 & 
          &       + 0.00146e0*cos3 - 0.00014e0*sin3 & 
          &       + 0.00031e0*cos4 - 0.00018e0*sin4 & 
          &       + 0.00069e0*cos5 

   !!    equation 5.19 and 5.20
     fliq = v*fb + (1.0e0 - v)*ft
     gliq = w*gb + (1.0e0 - w)*gt

   !!    equation 5.17
     sbrem = 0.5738e0*zbar*ye*t86*den * (tfac4*fliq - tfac5*gliq)
  end if




!!    recombination neutrino section
!!    for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
!!    equation 6.11 solved for nu
  xnum = 1.10520e8 * den * ye /(btemp*sqrt(btemp))
  nu   = bn_ifermi12(xnum)
  nu2  = nu * nu
  nu3  = nu2 * nu

!!    equation 6.7 and table 12
  zeta = 1.579e5 * zbar * zbar/btemp
  if (nu .ge. -20.0  .and. nu .lt. 0.0) then
     a1 = 1.51e-2
     a2 = 2.42e-1
     a3 = 1.21e0
     b  = 3.71e-2
     c  = 9.06e-1
     d  = 9.28e-1
     f1 = 0.0e0
     f2 = 0.0e0
     f3 = 0.0e0
  else if (nu .ge. 0.0  .and. nu .le. 10.0) then
     a1 = 1.23e-2
     a2 = 2.66e-1
     a3 = 1.30e0
     b  = 1.17e-1
     c  = 8.97e-1
     d  = 1.77e-1
     f1 = -1.20e-2
     f2 = 2.29e-2
     f3 = -1.04e-3
  end if

!!    equation 6.13  and 6.14
  if (nu .ge. -20.0  .and.  nu .le. 10.0) then
     z    = zeta/(1.0e0 + f1*nu + f2*nu2 + f3*nu3)  
     bigj = (a1/z + a2*z**(-2.25) + a3*z**(-4.55)) * exp(nu) & 
          &        / (1.0e0 + b*exp(c*nu)*(1.0e0 + d*z))

   !!    equation 6.5
     srecomb = tfac6 * 2.649e-18 * ye* zbar**13 * den * bigj & 
          &          / (exp(zeta + nu) + 1.0e0) 
  end if


!!    convert from erg/cm^3/s to erg/g/s 
!!    comment these out to duplicate the itoh et al plots
  deni   = 1.0e0/den
  spair   = spair*deni
  splas   = splas*deni
  sphot   = sphot*deni
  sbrem   = sbrem*deni
  srecomb = srecomb*deni


!!    the total neutrino loss rate
  sneut = splas + spair + sphot + sbrem + srecomb

  return

end subroutine bn_sneutx




