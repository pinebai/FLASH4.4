!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/states
!!
!! NAME
!! 
!!  states
!!
!! SYNOPSIS
!!
!!  call states(integer(IN) :: numIntCells,
!!              integer(IN) :: numCells,
!!              integer(IN) :: j,
!!              integer(IN) :: igeom,
!!              real(IN)    :: rho(numCells),
!!              real(IN)    :: u(numCells),
!!              real(INOUT) :: rhol(numCells),
!!              real(INOUT) :: rhor(numCells),
!!              real(INOUT) :: ul(numCells),
!!              real(INOUT) :: ur(numCells),
!!              real(INOUT) :: utl(numCells),
!!              real(INOUT) :: utr(numCells),
!!              real(INOUT) :: uttl(numCells),
!!              real(INOUT) :: uttr(numCells),
!!              real(IN)    :: p(numCells),
!!              real(INOUT) :: pl(numCells),
!!              real(INOUT) :: pr(numCells),
!!              real(IN) ::    gamcl(numCells),
!!              real(IN) ::    gamcr(numCells),
!!              real(IN)    :: ugrid(numCells),
!!              real(IN) ::    ce(numCells),
!!              real(IN) ::    game(numCells),
!!              real(IN) ::    gamer(numCells),
!!              real(IN) ::    gamc(numCells),
!!              real(IN) ::    gamel(numCells),
!!              real(IN)    :: eintl(numCells),
!!              real(IN)    :: eintr(numCells),
!!              real(IN)    :: xnl(numCells, hy_numXn),
!!              real(IN)    :: xnr(numCells, hy_numXn),
!!              real(IN)    :: dtdx(numCells),
!!              real(IN)    :: dt,
!!              real(IN)    :: x(numCells),
!!              real(IN)    :: xl(numCells),
!!              real(IN)    :: radial_coord(numCells),
!!              real(IN) ::    grav(numCells),
!!              real(IN) ::    fict(numCells))
!!
!! 
!! DESCRIPTION
!!  
!!  Computes effective left and right states for input to Riemann
!!  problems in PPM.  Naive guesses for these states based on the
!!  left and right limits of the cubic interpolants are corrected
!!  by taking into account information which reaches each cell
!!  interface via the different characteristics.  This is done by
!!  computing the average of each interpolated variable over the
!!  domain of dependence corresponding to each characteristic
!!  speed.
!!
!! SIDE EFFECTS
!!
!!  Modifies the following arrays exported by Hydro_data:
!!   hy_clft,    hy_crght,
!!   hy_plft,    hy_prght,
!!   hy_ulft,    hy_urght,
!!   hy_vlft,    hy_vrght,
!!   hy_utlft,   hy_utrght,
!!   hy_uttlft,  hy_uttrgt,
!!   hy_gmelft,  hy_gmergt,
!!   hy_gmclft,  hy_gmcrgt,
!!   hy_eiLft,   hy_eiRght
!!
!! ARGUMENTS
!!
!! numIntCells :
!! numCells :
!! j :
!! igeom :
!! rho :
!! u :
!! rhol :
!! rhor :
!! ul :
!! ur :
!! utl :
!! utr :
!! uttl :
!! uttr :
!! p :
!! pl :
!! pr :
!! gamcl :
!! gamcr :
!! ugrid :
!! ce :
!! game :
!! gamer :
!! gamc :
!! gamel :
!! eintl :
!! eintr :
!! xnl :
!! xnr :
!! dtdx :
!! dt :
!! x :
!! xl :
!! radial_coord :
!! grav :
!! fict :
!!
!!
!!***
subroutine states (numIntCells, numCells, &
                   j, igeom,&
                   rho, u, rhol, rhor, ul, ur, &
                   utl, utr, uttl, uttr, p, pl, pr, &
                   gamcl, gamcr, &
                   ugrid, ce, game, gamer, gamc, gamel, &
                   eintl, eintr, &
                   xnl, xnr, &
                   dtdx, dt, &
                   x, xl, radial_coord, grav, fict)



!=====================================================================
!     compute left and right states for input to rieman problem
  use Hydro_data, ONLY : hy_numXn
  use Hydro_data, ONLY : hy_clft, hy_plft, hy_ulft, hy_utlft, hy_gmelft, hy_gmclft, &
                         hy_dp, hy_p6, hy_du, hy_u6, hy_drho, hy_rho6, &
                         hy_gravr, hy_dgrav, hy_grav6, hy_smallp, hy_smlrho, &
                         hy_vlft, hy_uttlft, hy_xnlft, hy_crght, hy_prght, hy_urght, &
                         hy_vrght, hy_utrght, hy_uttrgt, hy_gmergt, hy_gmcrgt, &
                         hy_xnrght, hy_deint, hy_eint6, hy_eiLft, hy_eiRght, &
                         hy_dut, hy_dutt, hy_utt6, hy_dgame, hy_game6, hy_dgamc, hy_gamc6, &
                         hy_dxn, hy_xn6, hy_gravl, hy_ut6, &
                         hy_ppmModifystates, hy_leveque, &
                         hy_pwl, hy_pwr, hy_dpw, hy_pw6l, hy_pw6r, hy_pwcubic 

  implicit none


!!! ------ARGUMENTS -------------------------
  integer, intent(IN):: j, igeom, numIntCells, numCells
  real, intent(IN) :: dt
  real, intent(IN), dimension(numCells) :: rho, u, p, ugrid, x, xl, &
       radial_coord
  real, intent(INOUT), dimension(numCells) :: rhol, rhor, ul, ur
  real, intent(INOUT), dimension(numCells) :: utl, utr, uttl, uttr, pl, pr
  real, intent(IN),    dimension(numCells) :: gamcl, gamc, gamcr
  real, intent(IN),    dimension(numCells) :: ce, game, gamer, gamel
  real, intent(IN),    dimension(numCells) :: eintl, eintr
  real, intent(IN),    dimension(numCells) :: dtdx
  real, intent(IN),    dimension(numCells) :: grav, fict
  real, intent(INOUT), dimension(numCells, hy_numXn) :: xnl, xnr
  


!!! ---------LOCAL -----------------------------
  real, dimension(numCells) :: pallpl, pallml
  real, dimension(numCells) :: ppl, pml, p0l
  real, dimension(numCells) :: upl, uml, u0l
  real, dimension(numCells) :: ut0l, utt0l
  real, dimension(numCells) :: rhopl, rhoml, rho0l
  real, dimension(numCells) :: gravpl, gravml
  real, dimension(numCells) :: game0l, gamc0l
  real, dimension(numCells) :: eint0l
  
  real, dimension(numCells, hy_numXn) :: xn0l
  
  real, dimension(numCells) :: scrch1, scrch2, scrch3, scrch4,&
                               dloga,cflno,urel



!!  local
  integer :: numIntCells5, numIntCells6, numIntCells7, numIntCells8
  integer :: i, n

  real :: eta
  real :: clft_inv, crght_inv
  real :: qdt, hdt, fdt, slamm, slamp, scrch1_,  scrch2_

  real, PARAMETER :: forthd = 4.e00 / 3.e00
    
!! These are for implementing LeVeque's HSE method.
!! ref: LeVeque and Bale, Proc 7th Int'l Conf on Hyperbolic Problems

  real :: capGamma, m, h, ei
  real :: uplus, uminus, udelta
  real :: rplus, rminus, rdelta
  real :: pplus, pminus, pdelta
  real :: eplus, eminus
  real :: dele, delrho


  numIntCells5 = numIntCells+5
  numIntCells6 = numIntCells+6
  numIntCells7 = numIntCells+7
  numIntCells8 = numIntCells+8
  
  qdt  = 0.25e0*dt
  hdt  = 0.50e0*dt

!======================================================================
!  Compute fluid velocities relative to the grid, which may be moving.

  do i = 1,numIntCells8
     urel(i)  = u(i) - ugrid(i)
     cflno(i) = 0.e0
  enddo

!!
!!   LeVeque and Bale, 1988    Proc 7th Intl Conf on Hyperbolic Problems 
!!   but the expressions below are closer to ``Accurate Simulation of
!!   Rayleigh-Taylor-Instabilities'', Ralf Deiterding, Technical University
!!   Cottbus.   Some of L&B's expressions aren't dimensionally correct??!
!!
  if (hy_leveque) then 
     do i = 1,numIntCells8
        capGamma = (game(i)*game(i)-game(i)+2.)/(game(i)-1.)
        m     = rho(i)*u(i)
        ei    = p(i)/(game(i)-1.)+(m*m/(2.*rho(i)))
        if(i<numCells) h = xl(i+1)-xl(i)
        
        delrho = grav(i)*h/((game(i)-1.)*(m/rho(i)*m/rho(i)*capGamma - &
             &      2.*game(i)*ei/rho(i)))
        dele   = -rho(i)*grav(i)*h/(2.*ei*(game(i)-1))* &
             &              ( 1. - m*m*(3.-game(i))/   &
             & ((game(i)-1.)*(capGamma*m*m-2.*ei*game(i)*rho(i))))
        
        eplus  = ei*(1. + dele)
        eminus = ei*(1. - dele)
        
        rplus  = rho(i)*(1. + delrho)
        rminus = rho(i)*(1. - delrho)
        rdelta = rho(i)*delrho
        
        pplus  = (game(i)-1.)*(eplus  - m*m/(2.*rplus))
        pminus = (game(i)-1.)*(eminus - m*m/(2.*rminus))
        pdelta = (pplus-pminus)/2.
        
        uplus  = m/rplus
        uminus = m/rminus
        udelta = (uplus-uminus)/2.
        
        pl(i) = pl(i) + pdelta
        pr(i) = pr(i) - pdelta
        
        rhol(i) = rhol(i) + rdelta
        rhor(i) = rhor(i) - rdelta
        
        ul(i) = ul(i) + udelta
        ur(i) = ur(i) - udelta
     enddo
  endif

! 1. FOR LEFT STATES : FRYXELL SECTION 3.1.5, EQNS 58 - 59
!----------------------------------------------------------------------
!       Compute averages for the left side of each interface.

!       scrch#() here are domain factors used to compute the interpolant
!       averages.  See Colella & Woodward (1984), eq. 1.12.

!       Compute averages over the domain of dependence of the u+c characteri-
!       stic, if it reaches the left sides of the interfaces.

!
!       the pressure is calculated twice.  One will only contain
!       the `wave generating' pressure -- pressure in excess of 
!       hydrostatic -- if ppm_modifystates is .true.   The other
!       contains the full pressure; this is needed for the sound speed.
!
!       p{l,r},hy_dp,hy_p6 are the coefficients for the reconstruction of
!                    the total pressure;
!       pw{l,r},hy_dpw,pw6{l,r},hy_pwcubic are the coefficients for the 
!                    reconstruction of the wave-generating pressure.
!

  ! 1-(a): u+c characteristic
  do i = 5, numIntCells5
     scrch1(i) = dtdx(i-1) * (urel(i-1) + ce(i-1)) ! eqn 59
     cflno(i)  = max (0.e00, scrch1(i))            ! eqn 59
     scrch1(i) = 0.5e00 * min (1.e00, cflno(i))    ! eqn 59: 1/2*beta
     scrch2(i) = 1.e00 - forthd * scrch1(i)        ! eqn 58: 1-2/3*beta

     ! now solve eqn 58
     ppl(i)   = hy_pwr(i-1) - &
          & scrch1(i) * (hy_dpw(i-1)   - scrch2(i) * hy_pw6l(i-1))    &
          & - (scrch1(i)**3)*hy_pwcubic(i)
     
     pallpl(i)   = pr(i-1) - &
          & scrch1(i) * (hy_dp(i-1)   - scrch2(i) * hy_p6(i-1))   
     
     upl(i)   = ur(i-1) - &
          & scrch1(i) * (hy_du(i-1)   - scrch2(i) * hy_u6(i-1))

     rhopl(i) = rhor(i-1) - &
          & scrch1(i) * (hy_drho(i-1) - scrch2(i) * hy_rho6(i-1))
     
     gravpl(i) = hy_gravr(i-1) - &
          & scrch1(i) * (hy_dgrav(i-1) - scrch2(i) * hy_grav6(i-1))
     
     rhopl(i)    = max (hy_smlrho, rhopl(i))
     pallpl(i)   = max (hy_smallp, pallpl(i))
  end do
  
!       Compute averages over the domain of dependence of the u-c characteri-
!       stic, if it reaches the left sides of the interfaces.

  ! 1-(b): u-c characteristic
  do i = 5, numIntCells5
     scrch3(i) = dtdx(i-1) * (urel(i-1) - ce(i-1))            ! eqn 59
     scrch1(i) = 0.5e00 * min (1.e00, max (0.e00, scrch3(i))) ! eqn 59: 1/2*beta
     scrch2(i) = 1.e00 - forthd * scrch1(i)                   ! eqn 58: 1-2/3*beta
     
     pml(i)   = hy_pwr(i-1) - &
          & scrch1(i) * (hy_dpw(i-1)   - scrch2(i) * hy_pw6l(i-1))    &
          & - (scrch1(i)**3)*hy_pwcubic(i)

     uml(i)   = ur(i-1) - &
          & scrch1(i) * (hy_du(i-1)   - scrch2(i) * hy_u6(i-1))

     rhoml(i) = rhor(i-1) - &
          & scrch1(i) * (hy_drho(i-1) - scrch2(i) * hy_rho6(i-1))
     
     gravml(i) = hy_gravr(i-1) - &
          & scrch1(i) * (hy_dgrav(i-1) - scrch2(i) * hy_grav6(i-1))
  end do

  ! 1-(c): u characteristic
  do i = 5, numIntCells5
     scrch4(i) = dtdx(i-1) * urel(i-1)
     scrch1(i) = 0.5e00 * min (1.e00, max (0.e00, scrch4(i)))
     scrch2(i) = 1.e00 - forthd * scrch1(i)
     
     p0l(i)    = hy_pwr(i-1) &
          & - scrch1(i) * (hy_dpw(i-1)    - scrch2(i) * hy_pw6l(i-1)) &
          & - (scrch1(i)**3)*hy_pwcubic(i)
     
     rho0l(i)  = rhor(i-1) &
          & - scrch1(i) * (hy_drho(i-1)  - scrch2(i) * hy_rho6(i-1))
     u0l(i)    = ur (i-1) &
          & - scrch1(i) * (hy_du(i-1)    - scrch2(i) * hy_u6(i-1))
     ut0l(i)   = utr(i-1) &
          & - scrch1(i) * (hy_dut(i-1)   - scrch2(i) * hy_ut6(i-1))
     utt0l(i)  = uttr(i-1) &
          & - scrch1(i) * (hy_dutt(i-1)  - scrch2(i) * hy_utt6(i-1))
     game0l(i) = gamer(i-1) &
          &- scrch1(i) * (hy_dgame(i-1) - scrch2(i) * hy_game6(i-1))
     gamc0l(i) = gamcr(i-1) &
          &- scrch1(i) * (hy_dgamc(i-1) - scrch2(i) * hy_gamc6(i-1))
  end do

  do i = 5, numIntCells5
     eint0l(i) = eintr(i-1) &
          & - scrch1(i) * (hy_deint(i-1) - scrch2(i) * hy_eint6(i-1))
  end do
  
  do n = 1, hy_numXn 
     do i = 5, numIntCells5
        xn0l(i,n) = xnr(i-1,n) &
             & - scrch1(i) * (hy_dxn(i-1,n) - scrch2(i) * hy_xn6(i-1,n))
     end do
  end do

! 2. FOR LEFT STATES : FRYXELL SECTION 3.1.5, EQNS 62 - 67
!-------------------------------------------------------------------------------
!       Now compute the modified states for the left sides of cell interfaces,
!       excluding the effects of source terms and geometry, which are taken
!       into account later.  The computation of the Lagrangian sound speed
!       (hy_clft) takes into account a bug correction relayed by E. Muller (was
!       using gamc(i), should use gamc(i-1)).

!       scrch1() and scrch2() here are, respectively, C*(beta+ +/- beta-) and
!       beta0.  See Colella & Woodward (1984), eq. 3.6 and 3.7.  The arrays
!       scrch3() and scrch4() are domain factors (the ratio of the size of the
!       domain of dependence to the size of the zone) for the u-c and u charac-
!       teristics, respectively.

  do i = 5, numIntCells5
     hy_clft(i)  = sqrt (gamc(i-1) * pallpl(i) * rhopl(i))       ! (rho*c)
     clft_inv = 1.e0 / hy_clft(i)                                ! 1./(rho*c)
     
     scrch1_   = 0.5e00 * (upl(i) - uml(i) - (ppl(i) - pml(i)) & ! C*beta^{-}_l -- a part of eqn Fry 62
          * clft_inv)                                            ! NOTE: beta^{+}_l = 0.
     
     scrch2_   = (ppl(i) - p0l(i)) * clft_inv * clft_inv &       ! eqn Fry 63
          + 1.e00/rhopl(i)
     
     scrch2_   = scrch2_ - 1.e00 / rho0l(i)                      ! beta^{0}, eqn Fry 63
     
     if (-scrch3(i) .GE. 0.e0) scrch1_   = 0.e0
     if (-scrch4(i) .GE. 0.e0) scrch2_   = 0.e0
     
     hy_plft(i)   = ppl(i) + hy_clft(i) * scrch1_
     hy_plft(i)   = max (hy_plft(i), hy_smallp)
     
     hy_ulft(i)   = upl(i) - scrch1_ + hdt*fict(i-1)
     
     hy_vlft(i)   = 1.e00 / rhopl(i) - scrch2_ - scrch1_ * clft_inv !v=1/rho
     hy_vlft(i)   = min (hy_vlft(i), 1.e0/hy_smlrho)
     
     hy_utlft(i)  = ut0l(i)
     hy_uttlft(i) = utt0l(i)
     hy_gmelft(i) = game0l(i)
     hy_gmclft(i) = gamc0l(i)
  end do
  do i = 5, numIntCells5
     hy_eiLft(i) = eint0l(i)
  end do
  do n = 1, hy_numXn
     do i = 5, numIntCells5
        hy_xnlft(i,n) = xn0l(i,n)
     end do
  end do
  
!       Now correct the left states to include geometry terms.
  
  
  if (igeom .eq. 1 .or. igeom .eq. 2 .or. igeom .eq. 4)   then
     if (igeom .eq. 1)   then
        do i = 4, numIntCells5
           dloga(i) = 1.0 / x(i)
        end do
     else if (igeom .eq. 2) then
        do i = 4, numIntCells5
           dloga(i) = 2.0 / x(i)
        end do
     else
        do i = 4, numIntCells5
           dloga(i) = cos(x(i))/(sin(x(i)) * radial_coord(j))
        enddo
     end if
     do i = 4, numIntCells5
        scrch1_  = (abs (u(i) - ugrid(i)) + ce(i)) * dtdx(i)
        eta      = (1.e00 - scrch1_) / (ce(i) * dt * abs (dloga(i)))
        eta      = min (eta, 1.e00)
        dloga(i) = eta * dloga(i)
     end do
  
     do i = 5, numIntCells5
        scrch1_   = 0.5e00 * rho(i-1) * u(i-1) * dt * dloga(i-1)
        hy_vlft(i)   = 1.e00 / hy_vlft(i) - scrch1_
        hy_vlft(i)   = 1.e00 / hy_vlft(i)
        hy_plft(i)   = hy_plft(i) - scrch1_ * ce(i-1)**2
        hy_plft(i)   = max (hy_plft(i), hy_smallp)
        hy_vlft(i)   = min (hy_vlft(i), 1.e0/hy_smlrho)
     end do
  end if
  
  if (.not. (hy_ppmModifystates .or. hy_leveque)) then
     
     do i = 5, numIntCells5
        slamm = urel(i-1) - ce(i-1)
        slamp = urel(i-1) + ce(i-1)
        
        if ( slamm.le.0.e0 ) gravml(i) = 0.e0
        if ( slamp.le.0.e0 ) gravpl(i) = 0.e0
        
        if ( slamm.gt.0.e0 .and. slamp.gt.0.e0 ) then
           fdt = qdt
        else
           fdt = hdt
        end if
        
        hy_ulft(i) = hy_ulft(i) + fdt*(gravpl(i)+gravml(i))
     end do
     
  end if

! 3. FOR RIGHT STATES : FRYXELL SECTION 3.1.5, EQNS 60 - 61
!-------------------------------------------------------------------------------
!       Compute averages for the right side of each interface.

!       scrch#() here are geometry factors used to compute the interpolant
!       averages.  See Colella & Woodward (1984), eq. 1.12.

!       Compute averages over the domain of dependence of the u+c characteri-
!       stic, if it reaches the right sides of the interfaces.

!
!       the pressure is calculated twice.  One will only contain
!       the `wave generating' pressure -- pressure in excess of 
!       hydrostatic -- if ppm_modifystates is .true.   The other
!       contains the full pressure; this is needed for the sound speed.
!
!       p{l,r},hy_dp,hy_p6 are the coefficients for the reconstruction of
!                    the total pressure;
!       p{l,r},hy_dpw,pw6{l,r},hy_pwcubic are the coefficients for the 
!                    reconstruction of the wave-generating pressure.
!

  ! 3-(a): u+c characteristic
  do i= 5, numIntCells5
     scrch3(i) = -dtdx(i) * (urel(i) + ce(i))
     scrch1(i) = 0.5e00 * min (1.e00, max (0.e00, scrch3(i)))
     scrch2(i) = 1.e00 - forthd * scrch1(i)
     
     ppl(i)   = hy_pwl(i)   + scrch1(i) * (hy_dpw(i)   + scrch2(i) * hy_pw6r(i)) &
          & + (scrch1(i)**3)*hy_pwcubic(i)
     
     upl(i)   = ul(i)   + scrch1(i) * (hy_du(i)   + scrch2(i) * hy_u6(i))
     rhopl(i) = rhol(i) + scrch1(i) * (hy_drho(i) + scrch2(i) * hy_rho6(i))
     
     gravpl(i) = hy_gravl(i) +  &
          & scrch1(i) * (hy_dgrav(i) + scrch2(i) * hy_grav6(i))
  end do

!       Compute averages over the domain of dependence of the u-c characteri-
!       stic, if it reaches the right sides of the interfaces.
  ! 3-(b): u-c characteristic
  do i = 5, numIntCells5
     scrch1(i) = -dtdx(i) * (urel(i) - ce(i))
     scrch1(i) = max (0.e00, scrch1(i))
     cflno(i)  = max (cflno(i), scrch1(i))
     scrch1(i) = 0.5e00 * min (1.e00, scrch1(i))
     scrch2(i) = 1.e00 - forthd * scrch1(i)
     
     pml(i)   = hy_pwl(i)   + scrch1(i) * (hy_dpw(i)   + scrch2(i) * hy_pw6r(i)) &
          & + (scrch1(i)**3)*hy_pwcubic(i)
     
     pallml(i)   = pl(i)   + scrch1(i) * (hy_dp(i)   + scrch2(i) * hy_p6(i))
     
     uml(i)   = ul(i)   + scrch1(i) * (hy_du(i)   + scrch2(i) * hy_u6(i))
     rhoml(i) = rhol(i) + scrch1(i) * (hy_drho(i) + scrch2(i) * hy_rho6(i))
     
     gravml(i) = hy_gravl(i) + &
          & scrch1(i) * (hy_dgrav(i) + scrch2(i) * hy_grav6(i))
     
     pallml(i) = max (pallml(i), hy_smallp)
     rhoml(i)  = max (rhoml(i), hy_smlrho)
  end do

  ! 3-(c): u characteristic
  do i = 5, numIntCells5
     scrch4(i) = -dtdx(i) * urel(i)
     scrch1(i) = 0.5e00 * min (1.e00, max (0.e00, scrch4(i)))
     scrch2(i) = 1.0e00 - forthd * scrch1(i)
     
     p0l(i)    = hy_pwl(i)   &
          &        + scrch1(i) * (hy_dpw(i)    + scrch2(i) * hy_pw6r(i)) &
          & + (scrch1(i)**3)*hy_pwcubic(i)
     
     rho0l(i)  = rhol(i) &
          &        + scrch1(i) * (hy_drho(i)  + scrch2(i) * hy_rho6(i))
     u0l(i)    = ul(i)   &
          &        + scrch1(i) * (hy_du(i)    + scrch2(i) * hy_u6(i))
     ut0l(i)   = utl(i)  &
          &        + scrch1(i) * (hy_dut(i)   + scrch2(i) * hy_ut6(i))
     utt0l(i)  = uttl(i) &
          &        + scrch1(i) * (hy_dutt(i)  + scrch2(i) * hy_utt6(i))
     game0l(i) = gamel(i) &
          &        + scrch1(i) * (hy_dgame(i) + scrch2(i) * hy_game6(i))
     gamc0l(i) = gamcl(i) &
          &        + scrch1(i) * (hy_dgamc(i) + scrch2(i) * hy_gamc6(i))       
  end do
  
  do i = 5, numIntCells5
     eint0l(i) = eintl(i) &
          &        + scrch1(i) * (hy_deint(i) + scrch2(i) * hy_eint6(i))
  end do

  do n = 1, hy_numXn
     do i = 5, numIntCells5
        xn0l(i,n) = xnl(i,n) &
             &        + scrch1(i) * (hy_dxn(i,n) + scrch2(i) * hy_xn6(i,n))
     end do
  end do
  
! 4. FOR RIGHT STATES : FRYXELL SECTION 3.1.5, EQNS 62 - 67
!-------------------------------------------------------------------------------
!       Now compute the modified states for the right sides of cell interfaces,
!       excluding the effects of source terms and geometry, which are taken
!       into account later.

!       scrch1() and scrch2() here are, respectively, C*(beta+ +/- beta-) and
!       beta0.  See Colella & Woodward (1984), eq. 3.6 and 3.7.  The arrays
!       scrch3() and scrch4() are domain factors (the ratio of the size of the
!       domain of dependence to the size of the zone) for the u+c and u charac-
!       teristics, respectively.

  do i = 5, numIntCells5
     hy_crght(i)  = sqrt (gamc(i) * pallml(i) * rhoml(i))           ! (rho*c)
     crght_inv = 1.e0 / hy_crght(i)                                 ! 1./(rho*c)
     
     scrch1_   = -0.5e00 * (uml(i) - upl(i) + (pml(i) - ppl(i))  &  ! C*beta^{+}_r -- a part of eqn Fry 62
          * crght_inv)                                              ! NOTE: beta^{-1}_r = 0.
     
     scrch2_   = (pml(i) - p0l(i)) * crght_inv * crght_inv  &       ! eqn Fry 63
          + 1.e00 / rhoml(i)
     
     scrch2_   = scrch2_ - 1.e00 / rho0l(i)                         ! beta^{0}_r, eqn Fry 63
     
     if (-scrch3(i) .GE. 0.e0) scrch1_   = 0.e0
     if (-scrch4(i) .GE. 0.e0) scrch2_   = 0.e0
     
     hy_prght(i)  = pml(i) + hy_crght(i) * scrch1_                  ! eqn Fry 65
     hy_prght(i)  = max (hy_prght(i), hy_smallp)

     hy_urght(i)  = uml(i) + scrch1_ + hdt*fict(i)                  ! eqn Fry 66
     
     hy_vrght(i)  = 1.e00 / rhoml(i) - scrch2_ - scrch1_ *crght_inv ! v=1/rho
     hy_vrght(i)  = min (hy_vrght(i), 1.e0/hy_smlrho)
     
     hy_utrght(i) = ut0l(i)
     hy_uttrgt(i) = utt0l(i)
     hy_gmergt(i) = game0l(i)
     hy_gmcrgt(i) = gamc0l(i)

  end do

  do i = 5, numIntCells5
     hy_eiRght(i) = eint0l(i)
  end do

  do n = 1, hy_numXn
     do i = 5, numIntCells5
        hy_xnrght(i,n) = xn0l(i,n)
     end do
  end do
  
  !       Now correct the right states to include geometry terms.


  if (igeom .eq. 1 .or. igeom .eq. 2 .or. igeom .eq. 4)   then
     do i = 5, numIntCells5
        scrch1_   = 0.5e00 * rho(i) * u(i) * dt * dloga(i)
        hy_vrght(i)  = 1.e00 / hy_vrght(i) - scrch1_
        hy_vrght(i)  = 1.e00 / hy_vrght(i)
        hy_prght(i)  = hy_prght(i) - scrch1_ * ce(i)**2
        hy_prght(i)  = max (hy_prght(i), hy_smallp)
        hy_vrght(i)  = min (hy_vrght(i), 1.e0/hy_smlrho)
     end do
  end if

  if (.not. (hy_ppmModifystates .or. hy_leveque)) then
     
     do i = 5, numIntCells5
        slamm = urel(i) - ce(i)
        slamp = urel(i) + ce(i)
        
        if ( slamm.ge.0.e0 ) gravml(i) = 0.e0
        if ( slamp.ge.0.e0 ) gravpl(i) = 0.e0
        
        if ( slamm.lt.0.e0 .and. slamp.lt.0.e0 ) then
           fdt = qdt
        else
           fdt = hdt
        end if
        
        hy_urght(i) = hy_urght(i) + fdt*(gravpl(i)+gravml(i))
     end do
     
  end if

  if(xl(5) .ne. 0.e00) return
  
  return
end subroutine states

