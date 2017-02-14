!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/intrfc
!!
!! NAME
!! 
!!  subroutine intrfc
!!
!!
!! SYNOPSIS
!!
!!  call intrfc(integer(IN) :: numIntCells,
!!              integer(IN) :: numCells,
!!              integer(IN) :: guard,
!!              real(IN)    :: rho(numCells),
!!              real(IN)    :: u(numCells),
!!              real(IN)    :: ut(numCells),
!!              real(IN)    :: utt(numCells),
!!              real(IN)    :: p(numCells),
!!              real(INOUT) :: rhol(numCells),
!!              real(INOUT) :: rhor(numCells),
!!              real(INOUT) :: ul(numCells),
!!              real(INOUT) :: ur(numCells),
!!              real(INOUT) :: utl(numCells),
!!              real(INOUT) :: utr(numCells),
!!              real(INOUT) :: uttl(numCells),
!!              real(INOUT) :: uttr(numCells),
!!              real(INOUT) :: pl(numCells),
!!              real(INOUT) :: pr(numCells),
!!              real(OUT)   :: vl(numCells),
!!              real(OUT)   :: vr(numCells),
!!              real(INOUT) :: gamcl(numCells),
!!              real(INOUT) :: gamcr(numCells),
!!              real(IN)    :: game(numCells),
!!              real(INOUT) :: gamel(numCells),
!!              real(INOUT) :: gamer(numCells),
!!              real(IN)    :: gamc(numCells),
!!              real(IN)    :: grav(numCells),
!!              real(IN)    :: eint(numCells),
!!              real(INOUT) :: eintl(numCells),
!!              real(INOUT) :: eintr(numCells),
!!              real(IN)    :: xn(numCells, hy_numXn),
!!              real(INOUT) :: xnl(numCells, hy_numXn),
!!              real(INOUT) :: xnr(numCells, hy_numXn),
!!              real(OUT)   :: v(numCells),
!!              real(IN)    :: dx(numCells),
!!              real(IN)    :: x(numCells))
!!
!! 
!! DESCRIPTION
!!  
!!  Calculate zone interface values of all variables.  Start by using a
!!  high order polynomial to interpolate the zone average values to the
!!  interface (Fry. Eq. 24).  The coefficients of this polynomial are 
!!  constructed (by coeff() ) such that the polynomial reproduces the 
!!  correct average values in each zone.
!!
!!  Once the coefficients are computed, the first guess at the interface 
!!  values are made, and monotonicity must be enforced.  We need to ensure 
!!  that the interface value does not fall outside the range of the 
!!  adjacent zone average values --  this procedure is handled by interp()
!!
!!  intrfc() then looks for contact discontinuities (through detect() ) and
!!  steepens the profiles if found.  This is done for the density and the
!!  mass fractions.
!!
!!  Next, shocks are flattened if they are too thin, via flaten().  Shocks
!!  that are too thin are not treated accurately, and oscillations can result. 
!!  Flattening is applied to any shocks which are only one zone wide.
!!  No flattening should be used for gravitational accelerations.
!!
!!  Finally, we need to ensure that the parabolic reconstruction of the
!!  zone data is monotonic -- each point in the parabolic is required to
!!  fall between the two zone interface values -- this procedure is 
!!  handled by monot().
!!
!! ARGUMENTS
!!
!!  numIntCells :
!!  numCells :
!!  guard :
!!  rho(numCells) :
!!  u(numCells) :
!!  ut(numCells) :
!!  utt(numCells) :
!!  p(numCells) :
!!  rhol(numCells) :
!!  rhor(numCells) :
!!  ul(numCells) :
!!  ur(numCells) :
!!  utl(numCells) :
!!  utr(numCells) :
!!  uttl(numCells) :
!!  uttr(numCells) :
!!  pl(numCells) :
!!  pr(numCells) :
!!  vl(numCells) :
!!  vr(numCells) :
!!  gamcl(numCells) :
!!  gamcr(numCells) :
!!  game(numCells) :
!!  gamel(numCells) :
!!  gamer(numCells) :
!!  gamc(numCells) :
!!  grav(numCells) :
!!  eint(numCells) :
!!  eintl(numCells) :
!!  eintr(numCells) :
!!  xn :
!!  xnl :
!!  xnr :
!!  v :
!!  dx(numCells) :
!!  x(numCells) :
!!
!!
!!***

  subroutine intrfc(sweepDir,numIntCells, numCells, guard, &
       &            rho, u, ut, utt, p, &
       &            rhol, rhor, ul, ur, &
       &            utl, utr, uttl, uttr, &
       &            pl, pr, vl, vr, gamcl, &
       &            gamcr, game, &
       &            gamel,gamer,gamc,grav, eint, eintl, eintr, xn, &
       &            xnl, xnr, v, dx, x,tmp)

  use Hydro_data,   ONLY:  hy_numXn, hy_smlrho, hy_smallx, hy_small, hy_gravl, &
                           hy_drho, hy_rho6, &
                           hy_du, hy_u6, hy_dut, hy_ut6, &
                           hy_dutt, hy_utt6, hy_dp, hy_dgame, hy_game6, &
                           hy_gravr, hy_p6, hy_dgamc, hy_gamc6, hy_dgrav, &
                           hy_grav6, hy_dxn, hy_xn6, &
                           hy_deint, hy_eint6, &
                           hy_pwcubic, hy_pwl, hy_pwr, hy_dpw, hy_pw6l, hy_pw6r, &
                           hy_ppmModifystates, &
                           hy_useSteepening, hy_useCmaFlattening,&
                           hy_epsiln, hy_omg1, hy_omg2, hy_charLimiting



  implicit none
!!------ARGUMENTS-------------------
#include "Flash.h"


  integer, intent(IN) :: sweepDir,numIntCells, numCells, guard
  real, intent(IN),    DIMENSION(numCells,hy_numXn) :: xn
  real, intent(INOUT), DIMENSION(numCells,hy_numXn) :: xnl, xnr
  real, intent(IN), DIMENSION(numCells) :: &
       rho, u, ut, utt, p, gamc, game, grav, eint, dx, x,tmp
  real, intent(INOUT), DIMENSION(numCells) :: &
       rhol, rhor, &
       ul, ur, &
       utl, utr, &
       uttl, uttr, &
       pl, pr,&
       gamcl, gamcr, &
       gamel, gamer, &
       eintl, eintr
  real, intent(OUT), DIMENSION(numCells) :: &
       v, vl, vr
       

!! ----LOCAL ----------------
  real,dimension(numCells) :: rhog, rhogl, rhogr, drhog, rg6,&
                              coeff1, coeff2, coeff3, coeff4, coeff5, &
                              flatn, flatn1

  
  integer :: i, j, k, n, numIntCells5
  
  real :: checkl, checkr, check, &
       dcheckl, dcheckr, dcheck  

  logical :: charLimiting

  numIntCells5 = numIntCells + 5 
  
! get the coefficients of the quartic polynomial through the zones neighboring
! each zone.  The polynomial used is Eq. 1.6 in Colella & Woodward. 

  do i = 1, numCells
     coeff1(i) = 0.e0
     coeff2(i) = 0.e0
     coeff3(i) = 0.e0     
     coeff4(i) = 0.e0
     coeff5(i) = 0.e0
     rhog(i) = 0.e0
     rhogl(i) = 0.e0
     rhogr(i) = 0.e0     
     drhog(i) = 0.e0
     rg6(i) = 0.e0
     flatn(i) = 0.e0
     flatn1(i) = 0.e0
  end do

  call coeff(numIntCells,numCells, dx, coeff1, coeff2, coeff3, coeff4, coeff5)


! interpolate the density to find the interface values, rhol and rhor for
! each zone


  if (hy_charLimiting) then
     ! Apply the limiting using the characteristic variables - this method project
     ! the primitive variables (density, velocity fields, and pressure) onto the
     ! chracteristic variables, apply limitings, and then project them back to the
     ! primitive variables.

     call interp_char(sweepDir, numIntCells, numCells, &
                      rhol, rho, rhor, &
                      ul,   u,   ur,   &
                      utl,  ut,  utr,  &
                      uttl, utt, uttr, &
                      pl,   p,   pr,   &
                      gamc, game,      &
                      coeff1,coeff2,coeff3,coeff4,coeff5)

     if (hy_useSteepening) &
          call detect(numIntCells,numCells,rhol, rho, rhor, &
                      hy_smlrho, rho, p, game, dx, x)


     ! interpolate the mass fractions and look for contacts
     do n = 1, hy_numXn

        call interp(numIntCells, numCells,xnl(1,n), xn(1,n), xnr(1,n),  &
                    coeff1, coeff2, coeff3, coeff4, coeff5)
     
        if (hy_useSteepening .and. n <= NSPECIES) & 
             call detect(numIntCells, numCells,xnl(1,n), xn(1,n), xnr(1,n), hy_smallx,  &
                         rho, p, game, dx, x)
     
     end do

  else

     ! Apply the limiting using the primitive variables - original version of PPM
     call interp(numIntCells,numCells,rhol, rho, rhor, coeff1, coeff2, coeff3, coeff4, coeff5)

     ! search for contact discontinuities, and steepen the profile if we find them
     ! -- this prevents the contact from spreading out over too many zones as it
     ! propagates
     if (hy_useSteepening) &
          call detect(numIntCells,numCells,rhol, rho, rhor, &
                      hy_smlrho, rho, p, game, dx, x)
    

     ! interpolate the mass fractions and look for contacts
     do n = 1, hy_numXn
     
        call interp(numIntCells, numCells,xnl(1,n), xn(1,n), xnr(1,n),  &
                    coeff1, coeff2, coeff3, coeff4, coeff5)
     
        if (hy_useSteepening .and. n <= NSPECIES) & 
             call detect(numIntCells, numCells,xnl(1,n), xn(1,n), xnr(1,n), hy_smallx,  &
                         rho, p, game, dx, x)
     
     end do
      
     ! interpolate the remainder of the variables 
     call interp(numIntCells,numCells,ul,   u,   ur,   coeff1,coeff2,coeff3,coeff4,coeff5)
     call interp(numIntCells,numCells,utl,  ut,  utr,  coeff1,coeff2,coeff3,coeff4,coeff5)
     call interp(numIntCells,numCells,uttl, utt, uttr, coeff1,coeff2,coeff3,coeff4,coeff5)
     call interp(numIntCells,numCells,pl,   p,   pr,   coeff1,coeff2,coeff3,coeff4,coeff5)

  endif

  ! interpolate the adiabatic indices (game and gamc) and gravity
  call interp(numIntCells,numCells,gamel,   game,gamer,   coeff1,coeff2,coeff3,coeff4,coeff5)
  call interp(numIntCells,numCells,gamcl,   gamc,gamcr,   coeff1,coeff2,coeff3,coeff4,coeff5)
  call interp(numIntCells,numCells,hy_gravl,grav,hy_gravr,coeff1,coeff2,coeff3,coeff4,coeff5)
  call interp(numIntCells, numCells,eintl, eint, eintr,&
              coeff1,coeff2,coeff3,coeff4,coeff5)


  ! look for shocks and flatten the structure if they are too thin    
  call flaten(numIntCells, numCells, u, p, flatn, flatn1)

  do i = 4, numIntCells5
     rhol (i) = flatn(i) * rho (i) + flatn1(i) * rhol (i)
     rhor (i) = flatn(i) * rho (i) + flatn1(i) * rhor (i)
     
     ul   (i) = flatn(i) * u   (i) + flatn1(i) * ul   (i)
     ur   (i) = flatn(i) * u   (i) + flatn1(i) * ur   (i)
     
     utl  (i) = flatn(i) * ut  (i) + flatn1(i) * utl  (i)
     utr  (i) = flatn(i) * ut  (i) + flatn1(i) * utr  (i)
     
     uttl (i) = flatn(i) * utt (i) + flatn1(i) * uttl (i)
     uttr (i) = flatn(i) * utt (i) + flatn1(i) * uttr (i)
     
     pl   (i) = flatn(i) * p   (i) + flatn1(i) * pl   (i)
     pr   (i) = flatn(i) * p   (i) + flatn1(i) * pr   (i)
     
     gamel(i) = flatn(i) * game(i) + flatn1(i) * gamel(i)
     gamer(i) = flatn(i) * game(i) + flatn1(i) * gamer(i)
     
     gamcl(i) = flatn(i) * gamc(i) + flatn1(i) * gamcl(i)
     gamcr(i) = flatn(i) * gamc(i) + flatn1(i) * gamcr(i)

     eintl(i) = flatn(i) * eint(i) + flatn1(i) * eintl(i)
     eintr(i) = flatn(i) * eint(i) + flatn1(i) * eintr(i)
  end do
    
  do n = 1, hy_numXn
     do  i = 4, numIntCells5
        xnl(i,n) = flatn(i) * xn(i,n) + flatn1(i) * xnl(i,n)
        xnr(i,n) = flatn(i) * xn(i,n) + flatn1(i) * xnr(i,n)
     end do
  end do


  ! make sure every thing is monotonic -- no new maxima or minima
  ! should have been introduced (see Colella & Woodward Eq. 1.10)

  call monot(numIntCells, numCells,rhol,     rho,  rhor,     hy_drho,  hy_rho6 )
  call monot(numIntCells, numCells,ul,       u,    ur,       hy_du,    hy_u6   )
  call monot(numIntCells, numCells,utl,      ut,   utr,      hy_dut,   hy_ut6  )
  call monot(numIntCells, numCells,uttl,     utt,  uttr,     hy_dutt,  hy_utt6 )
  call monot(numIntCells, numCells,pl,       p,    pr,       hy_dp,    hy_p6   )
  call monot(numIntCells, numCells,gamel,    game, gamer,    hy_dgame, hy_game6)
  call monot(numIntCells, numCells,gamcl,    gamc, gamcr,    hy_dgamc, hy_gamc6)
  call monot(numIntCells, numCells,eintl,    eint, eintr,    hy_deint, hy_eint6)
  call monot(numIntCells, numCells,hy_gravl, grav, hy_gravr, hy_dgrav, hy_grav6)
  
  do n = 1, hy_numXn
     call monot(numIntCells, numCells,xnl(1,n),xn(1,n),xnr(1,n),hy_dxn(1,n), hy_xn6(1,n))
  end do
    
  do i = 4, numIntCells5
     vl(i) = 1.e0 / rhol(i)
     v (i) = 1.e0 / rho (i)
     vr(i) = 1.e0 / rhor(i)
  end do

! if we are doing the modified states version of PPM, subtract the 
! rho*g contribution off the pressure.
!
!  eqns (27)--(31), Zingale et al (2002) ApJ
!

  hy_pwl(4:numIntCells5)  = pl(4:numIntCells5)
  hy_pwr(4:numIntCells5)  = pr(4:numIntCells5)
  
  if (hy_ppmModifystates) then

! for the modified states formalism, interpolate the quantity rho*g, which
! will be removed from the pressure when computing the wave structure

     rhog = rho*grav
     
     call interp(numIntCells, numCells,rhogl, rhog, rhogr, &
                 coeff1,coeff2,coeff3,coeff4,coeff5)
     
     do i = 4, numIntCells5
        rhogl(i) = flatn(i) * rhog(i) + flatn1(i) * rhogl(i)
        rhogr(i) = flatn(i) * rhog(i) + flatn1(i) * rhogr(i)
     end do
     
     call monot(numIntCells, numCells,rhogl, rhog, rhogr, drhog, rg6)
     
     hy_pwcubic = dx/3.e0*rg6
     hy_pw6r    = hy_p6 + 0.5e0*dx*(drhog+rg6)
     hy_pw6l    = hy_p6 + 0.5e0*dx*(drhog-rg6)
     hy_dpw     = hy_dp - dx*(rhogl + 0.5e0*(drhog+rg6))
     
  else
     
     hy_pwcubic = 0.e0
     hy_dpw     = hy_dp
     hy_pw6l    = hy_p6
     hy_pw6r    = hy_p6
     
  endif

! if we are doing the CMA flattening of the abundances, do it now
  if (hy_useCmaFlattening) &
       call cma_flatten(numIntCells, numCells, guard,xn, xnl, xnr, hy_dxn, hy_xn6)


end subroutine intrfc

