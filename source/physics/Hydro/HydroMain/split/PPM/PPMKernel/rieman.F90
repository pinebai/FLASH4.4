!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/rieman
!!
!! NAME
!! 
!!  rieman
!!
!! SYNOPSIS
!!
!!       call rieman ( integer(IN) :: numIntCells, 
!!                     integer(IN) :: numCells, 
!!                     real(OUT)   :: rhoav(numCells), 
!!                     real(OUT)   :: uav(numCells), 
!!                     real(OUT)   :: utav(numCells), 
!!                     real(OUT)   :: uttav(numCells), 
!!                     real(OUT)   :: pav(numCells), 
!!                     real(OUT)   :: urell(numCells), 
!!                     real(IN)    :: ugrdl(numCells), 
!!                     real(IN)    :: game(numCells), 
!!                     real(OUT)   :: gameav(numCells), 
!!                     real(OUT)   :: eintAv(numCells),
!!                     real(OUT)   :: xnav(numCells, hy_numXn), 
!!                     real(IN) :: x(numCells))
!!                 
!!
!!
!! DESCRIPTION
!!  
!!  Solve riemann shock tube problem for a general equation of state using 
!!  the method of Colella and Glaz.  Use a two shock approximation, and
!!  linearly interpolation between the head and tail of a rarefaction to
!!  treat rarefactions.
!!
!!  Take as input the effective left and right states, obtained by 
!!  integrating the parabolic reconstructions of the data over the domain
!!  of dependence of each of the characteristics on the respective side
!!  of the interface.  This is accomplished by states().  Return the
!!  solution to Riemann's problem on axis -- this is used in computing
!!  the fluxes.
!!
!!  The Riemann problem for the Euler's equation produces 4 states, 
!!  separated by the three characteristics (u - cs, u, u + cs):
!!
!!
!!        l_1      t    l_2       l_3
!!         \       ^   .       /
!!          \  *L  |   . *R   /
!!           \     |  .     /
!!            \    |  .    /
!!        L    \   | .   /    R
!!              \  | .  /
!!               \ |. / 
!!                \|./
!!       ----------+----------------> x
!!
!!       l_1 = u - cs   eigenvalue
!!       l_2 = u        eigenvalue (contact)
!!       l_3 = u + cs   eigenvalue
!!
!!       only density jumps across l_2
!!
!!  References:
!!
!!   CG:   Colella & Glaz 1985, JCP, 59, 264.
!!
!!   CW:   Colella & Woodward 1984, JCP, 54, 174.
!!
!!   Fry:  Fryxell et al. 2000, ApJS, 131, 273.
!!
!!   Toro: Toro 1999, ``Riemann Solvers and Numerical Methods for Fluid
!!         Dynamcs: A Practical Introduction, 2nd Ed.'', Springer-Verlag
!!
!! ARGUMENTS
!!
!! numIntCells :
!! numCells :
!! rhoav :
!! uav :
!! utav :
!! uttav :
!! pav :
!! urell :
!! ugrdl :
!! game :
!! gameav :
!! eintAv :
!! xnav :
!! x :
!!
!!
!!
!!***
subroutine rieman (numIntCells, numCells, &
                   rhoav, uav, utav, uttav, pav, &
                   urell, ugrdl, game, gameav, eintAv, xnav, x)
                   

  use Hydro_data, ONLY: hy_numXn, &
                        hy_gmelft, hy_gmergt, &
                        hy_plft,   hy_prght,  &
                        hy_clft,   hy_crght,  &
                        hy_ulft,   hy_urght,  &
                        hy_vlft,   hy_vrght,  &
                        hy_utlft,  hy_utrght, &
                        hy_uttlft, hy_uttrgt, &
                        hy_eiLft,  hy_eiRght, &
                        hy_xnlft,  hy_xnrght, &
                        hy_smallp, hy_smallu, &
                        hy_smlrho, hy_nriem,  &
                        hy_gmclft, hy_gmcrgt,hy_pstor, &
                        hy_riemanTol
  use Driver_interface, ONLY : Driver_abortFlash


  implicit none
!! Arguments ---------------------------------- 

  integer, intent (IN) :: numIntCells,numCells
  real, intent(IN), DIMENSION(numCells) :: x
  real, intent(IN), DIMENSION(numCells) :: ugrdl, game
  real, intent(OUT), DIMENSION(numCells) :: uav, rhoav, utav, uttav, pav, &
                               urell, gameav, eintAv
  real, intent(OUT), DIMENSION(numCells,hy_numXn) :: xnav


!! Local variable ---------------------------------
  real, DIMENSION(numCells) :: wlft, wrght, pstar, ustar, vstar, cestar, &
       rhostr, westar, ps, us, uts, utts, vs, rhos, ces, ws, wes, &
       gmstar, games, gamcs

  
  real, DIMENSION(numCells) :: pstar1, pstar2, gmstrl, gmstrr, &
       &   wlft1, wrght1, gmin, gmax, &
       &   gamfac, aux
  real, dimension(numCells) :: scrch1, scrch2, scrch3, scrch4

  real  ::  ge, gc, ustrl1, ustrr1, ustrl2, ustrr2, &
       & delu1, delu2, pres_err
  
  integer :: i, j, k, n, numIntCells5, ierr
  
  
  character(len=1), save :: dirs(3) = (/ 'x', 'y', 'z' /)
  
  real, parameter :: small_dp = 1.e2 * epsilon(1.e0)

  
    !--------------------------------------------------------------------------

    ! We carry around two adiabatic indices (gammas).  gamc is usually
    ! referred to as gamma_1, and is d(log P)/d(log rho) (CG Eq. 7).  game
    ! is CG Eq. 8.      

    ! calculate limits on gamma based on the values in the neighboring zones
    ! gamfac is the max gamma factor in CG Eq. 31, and is used over and over,
    ! so store it for efficiency.

  numIntCells5 = numIntCells + 5

  do i = 5, numIntCells5
     aux(i)    = sqrt (0.5e0 * (game(i) - 1.0e0) / game(i))
     ge        = 0.5e0 * (hy_gmelft(i) + hy_gmergt(i))
     gc        = 0.5e0 * (hy_gmclft(i) + hy_gmcrgt(i))
     gamfac(i) = (1.e0 - ge / gc) * (ge - 1.e0)
     gmin(i)   = min (game(i-1), game(i), game(i+1))
     gmax(i)   = max (game(i-1), game(i), game(i+1))
  enddo
  
    ! construct first guess for secant iteration by assuming that the nonlinear 
    ! wave speed is equal to the sound speed -- the resulting expression is the
    ! same as Toro, Eq. 9.28 in the Primitive Variable Riemann Solver (PVRS).
    ! See also Fry Eq. 72.
    
  do i = 5, numIntCells5
     pstar1(i) = hy_prght(i) - hy_plft(i) - hy_crght(i) * (hy_urght(i) - hy_ulft(i))
     pstar1(i) = hy_plft(i) + pstar1(i) * (hy_clft(i) / (hy_clft(i) + hy_crght(i)))
     pstar1(i) = max (hy_smallp, pstar1(i))
  enddo

    ! calculate approximation jump in gamma acrosss the interface based on the 
    ! first guess for the pressure jump.  There is a left and right 'star' region,
    ! so we need gamma add both places.  Use CG Eq. 31 and 32, with definitions
    ! as in CG Eq. 33.
    
  do i = 5, numIntCells5
     gmstrl(i) = gamfac(i) * (pstar1(i) - hy_plft(i))
     gmstrl(i) = hy_gmelft(i) + 2.e0 * gmstrl(i) / (pstar1(i) + hy_plft(i))
     
     gmstrr(i) = gamfac(i) * (pstar1(i) - hy_prght(i))
     gmstrr(i) = hy_gmergt(i) + 2.e0 * gmstrr(i) / (pstar1(i) + hy_prght(i))
     
     gmstrl(i) = max (gmin(i), min( gmstrl(i), gmax(i)))
     gmstrr(i) = max (gmin(i), min( gmstrr(i), gmax(i)))
  enddo

    ! calculate nonlinear wave speeds for the left and right moving waves based
    ! on the first guess for the pressure jump.  Again, there is a left and a 
    ! right wave speed.  Compute this using CG Eq. 34.
    
  do i = 5, numIntCells5
     scrch1(i) = pstar1(i) - (gmstrl(i) - 1.e0) * hy_plft(i) &
          & / (hy_gmelft(i) - 1.e0)
     if (scrch1(i) .EQ. 0.e0) scrch1(i) = hy_smallp
     
     wlft1(i)  = pstar1(i) + 0.5e0 * (gmstrl(i) - 1.e0) &
          & * (pstar1(i) + hy_plft(i))
     wlft1(i)  = (pstar1(i) - hy_plft(i)) * wlft1(i) / (hy_vlft(i) * scrch1(i))
     wlft1(i)  = sqrt(abs(wlft1(i)))
     

     scrch2(i) = pstar1(i) - (gmstrr(i) - 1.e0) * hy_prght(i) /(hy_gmergt(i) - 1.e0)
     
     if (scrch2(i) .EQ. 0.e0) scrch2(i) = hy_smallp
     
     wrght1(i) = pstar1(i) + 0.5e0 * (gmstrr(i) - 1.e0) &
          & * (pstar1(i) + hy_prght(i))
     wrght1(i) = (pstar1(i) - hy_prght(i)) * wrght1(i) / (hy_vrght(i) * scrch2(i))
     wrght1(i) = sqrt(abs(wrght1(i)))
     
       ! if the pressure jump is small, the wave speed is just the sound speed

     if (abs (pstar1(i) - hy_plft(i)) < small_dp*(pstar1(i) + hy_plft(i))) wlft1(i) = hy_clft(i)
     wlft1(i)  = max (wlft1(i),  aux(i) * hy_clft(i))
     
     if (abs (pstar1(i) - hy_prght(i)) < small_dp*((pstar1(i) + hy_prght(i)))) wrght1(i) = hy_crght(i)
     wrght1(i) = max (wrght1(i), aux(i) * hy_crght(i))
  enddo

    ! construct second guess for the pressure using the nonlinear wave speeds
    ! from the first guess.  This is basically the same thing we did to get
    ! pstar1, except now we are using the better wave speeds instead of the 
    ! sound speed.

  do i = 5, numIntCells5
     pstar2(i) = hy_prght(i) - hy_plft(i) - wrght1(i) * (hy_urght(i) - hy_ulft(i))
     pstar2(i) = hy_plft(i) + pstar2(i) * wlft1(i) / (wlft1(i) + wrght1(i))
     pstar2(i) = max (hy_smallp, pstar2(i))
  enddo

    ! begin the secant iteration -- see CG Eq. 17 for details.  We will continue to
    ! interate for convergence until the error falls below tol (in which case, 
    ! things are good), or we hit hy_nriem iterations (in which case we have a 
    ! problem, and we spit out an error).

  do i = 5, numIntCells5
     
     hy_pstor(1) = pstar1(i)
     hy_pstor(2) = pstar2(i)
     
     do n = 1, hy_nriem
        
        ! new values for the gamma at the "star" state -- again, using CG Eq. 31
          
        gmstrl(i) = gamfac(i) * (pstar2(i) - hy_plft(i))
        gmstrl(i) = hy_gmelft(i) + 2.e0 * gmstrl(i) / (pstar2(i) + hy_plft(i))
        
        gmstrr(i) = gamfac(i) * (pstar2(i) - hy_prght(i))
        gmstrr(i) = hy_gmergt(i) + 2.e0 * gmstrr(i) / (pstar2(i) + hy_prght(i))
        
        gmstrl(i) = max (gmin(i), min (gmax(i), gmstrl(i)))
        gmstrr(i) = max (gmin(i), min (gmax(i), gmstrr(i)))
        
        ! new nonlinear wave speeds, using CG Eq. 34 and the updated gammas
          
        scrch1(i) = pstar2(i) - (gmstrl(i) - 1.e0) * hy_plft(i) &
             & / (hy_gmelft(i) - 1.e0)
        if (scrch1(i) .EQ. 0.e0) scrch1(i) = hy_smallp
        
        wlft(i)   = pstar2(i) + 0.5e0 * (gmstrl(i) - 1.e0) &
             & * (pstar2(i) + hy_plft(i))
        wlft(i)   = (pstar2(i) - hy_plft(i)) * wlft(i) / (hy_vlft(i) * scrch1(i))
        wlft(i)   = sqrt(abs(wlft(i)))

        scrch2(i) = pstar2(i) - (gmstrr(i) - 1.e0) * hy_prght(i) /(hy_gmergt(i) - 1.e0)

        if (scrch2(i) .EQ. 0.e0) scrch2(i) = hy_smallp
        
        wrght(i)  = pstar2(i) + 0.5e0 * (gmstrr(i) - 1.e0) &
             & * (pstar2(i) + hy_prght(i))
        wrght(i)  = (pstar2(i) - hy_prght(i)) * wrght(i) / (hy_vrght(i) * scrch2(i))
        wrght(i)  = sqrt(abs(wrght(i)))
        
        ! if the pressure jump is small, the wave speed is just the sound speed

        if (abs (pstar2(i) - hy_plft(i)) < small_dp*(pstar2(i) + hy_plft(i))) wlft(i) = hy_clft(i)
        wlft(i)  = max (wlft(i), aux(i) * hy_clft(i))
        
        if (abs (pstar2(i) - hy_prght(i)) < small_dp*(pstar2(i) + hy_prght(i))) wrght(i) = hy_crght(i)
        wrght(i) = max (wrght(i), aux(i) * hy_crght(i))

        ! compute the velocities in the "star" state -- using CG Eq. 18 -- ustrl2 and
        ! ustrr2 are the velocities they define there.  ustrl1 and ustrl2 seem to be
        ! the velocities at the last time, since pstar1 is the old 'star' pressure, and
        ! wlft1 is the old wave speed.
        
        ustrl1    =  hy_ulft(i) - (pstar1(i) -  hy_plft(i)) /  wlft1(i)
        ustrr1    = hy_urght(i) + (pstar1(i) - hy_prght(i)) / wrght1(i)
        ustrl2    =  hy_ulft(i) - (pstar2(i) -  hy_plft(i)) /   wlft(i)
        ustrr2    = hy_urght(i) + (pstar2(i) - hy_prght(i)) /  wrght(i)
        
        delu1     = ustrl1 - ustrr1
        delu2     = ustrl2 - ustrr2
        scrch1(i) = delu2  - delu1
        
        if (abs(pstar2(i)-pstar1(i)) .le. hy_smallp) scrch1(i) = 0.e0
        
        if (abs(scrch1(i)) .lt. hy_smallu) then
           delu2 = 0.e0
           scrch1(i) = 1.e0
        endif

        ! pressure at the "star" state -- using CG Eq. 18

        pstar(i)  = pstar2(i) - delu2 * (pstar2(i) - pstar1(i)) / scrch1(i)
        pstar(i)  = max (hy_smallp, pstar(i))
        
        ! check for convergence of iteration, hy_riemanTol is a run-time parameter
        
        pres_err = abs(pstar(i)-pstar2(i)) / pstar(i)
        if (pres_err .lt. hy_riemanTol) goto 10
        
        ! reset variables for next iteration
          
        pstar1(i) = pstar2(i)
        pstar2(i) = pstar(i)
        hy_pstor(n+2) = pstar(i)
        
        wlft1(i)  = wlft(i)
        wrght1(i) = wrght(i)
        
     enddo
     
     n = n - 1
     
     ! print error message and stop code if iteration fails to converge
     
     print *, ' '
     print *, 'Nonconvergence in subroutine rieman'
     print *, ' '
     print *, 'Zone index       = ', i
     print *, 'Zone center      = ', x(i)
     print *, 'Iterations tried = ', n+2
     print *, 'Pressure error   = ', pres_err
     print *, 'rieman_tol       = ', hy_riemanTol
     print *, ' '
     print *, 'pL       = ', hy_plft(i),   ' pR       =', hy_prght(i)
     print *, 'uL       = ', hy_ulft(i),   ' uR       =', hy_urght(i)
     print *, 'cL       = ', hy_clft(i),   ' cR       =', hy_crght(i)
     print *, 'gamma_eL = ', hy_gmelft(i), ' gamma_eR =', hy_gmergt(i)
     print *, 'gamma_cL = ', hy_gmclft(i), ' gamma_cR =', hy_gmcrgt(i)
     print *, ' '
     print *, 'Iteration history:'
     print *, ' '
     print '(A4, 2X, A20)', 'n', 'p*'
     do j = 1, n+2
        print '(I4, 2X, E20.12)', j, hy_pstor(j)
     enddo
     print *, ' '
     print *, 'Terminating execution.'
     call Driver_abortFlash('Nonconvergence in subroutine rieman')
       
       ! land here if the iterations have converged
       
10     continue
     
  enddo

! end of secant iteration

! calculate fluid velocity for the "star" state -- this comes from the shock
! jump equations, Fry Eq. 68 and 69.  The ustar velocity can be computed
! using either the jump eq. for a left moving or right moving shock -- we use
! the average of the two.
! NOTE: Also look at Fry Eqn. 75 and 76.

  do i = 5, numIntCells5
     scrch3(i) = hy_ulft (i) - (pstar(i) -  hy_plft(i)) /  wlft(i)
     scrch4(i) = hy_urght(i) + (pstar(i) - hy_prght(i)) / wrght(i)
     ustar(i)  = 0.5e0 * (scrch3(i) + scrch4(i))
  enddo

! account for grid velocity

  do i = 5, numIntCells5
     urell(i)   = ustar(i) - ugrdl(i)
!!$     scrch1(i)  = sign (1.e0, urell(i))
     scrch1(i)  = sign (0.5, urell(i)) - sign(0.5,-urell(i))     
  enddo

! decide which state is located at the zone iterface based on the values 
! of the wave speeds.  This is just saying that if ustar > 0, then the state
! is U_L.  if ustar < 0, then the state on the axis is U_R.

  do i = 5, numIntCells5
     
     scrch2(i) = 0.5e0 * ( 1.e0 + scrch1(i))
     scrch3(i) = 0.5e0 * ( 1.e0 - scrch1(i))
     
     ps(i)    = hy_plft(i)   * scrch2(i) + hy_prght(i)  * scrch3(i)
     us(i)    = hy_ulft(i)   * scrch2(i) + hy_urght(i)  * scrch3(i)
     uts(i)   = hy_utlft(i)  * scrch2(i) + hy_utrght(i) * scrch3(i)
     utts(i)  = hy_uttlft(i) * scrch2(i) + hy_uttrgt(i) * scrch3(i)
     vs(i)    = hy_vlft(i)   * scrch2(i) + hy_vrght(i)  * scrch3(i) !v for v=1/rho
     games(i) = hy_gmelft(i) * scrch2(i) + hy_gmergt(i) * scrch3(i)
     gamcs(i) = hy_gmclft(i) * scrch2(i) + hy_gmcrgt(i) * scrch3(i)
     
     rhos(i)  = 1.e0 / vs(i)
     rhos(i)  = max (hy_smlrho, rhos(i))
     
     vs(i)    = 1.e0 / rhos(i)
     ws(i)    = wlft(i) * scrch2(i) + wrght(i) * scrch3(i)
     ces(i)   = sqrt (gamcs(i) * ps(i) * vs(i))
     
     ! compute rhostar, using the shock jump condition (Fry Eq. 80)
     
     vstar(i)  = vs(i) - (pstar(i) - ps(i)) / ws(i) / ws(i)
     rhostr(i) = 1.e0 / vstar(i)
     cestar(i) = sqrt (gamcs(i) * pstar(i) * vstar(i))
     
! compute some factors, Fry Eq. 81 and 82       

     wes(i)    = ces(i)    - scrch1(i) * us(i)
     westar(i) = cestar(i) - scrch1(i) * ustar(i)
     
     scrch4(i) = ws(i) * vs(i) - scrch1(i) * us(i)
     
     
     if (pstar(i) - ps(i) .ge. 0.e0) then
        wes(i)    = scrch4(i)
        westar(i) = scrch4(i)
     endif
     
     wes(i)    = wes(i)    + scrch1(i) * ugrdl(i)
     westar(i) = westar(i) + scrch1(i) * ugrdl(i)
  enddo


  ! compute Fry Eq. 86
  do i = 5, numIntCells5
     gamfac(i) = (1.e0 - games(i) / gamcs(i)) * (games(i) - 1.e0)
     gmstar(i) = gamfac(i) * (pstar(i) - ps(i))
     gmstar(i) = games(i) + 2.e0 * gmstar(i) / (pstar(i) + ps(i))
     gmstar(i) = max (gmin(i), min (gmax(i), gmstar(i)))
  enddo

  do i = 5, numIntCells5
     eintAv(i) = hy_eiLft(i) * scrch2(i) + hy_eiRght(i) * scrch3(i)
  enddo
  
  do n = 1, hy_numXn 
     do i = 5, numIntCells5
        xnav(i,n) = hy_xnlft(i,n) * scrch2(i) + hy_xnrght(i,n) * scrch3(i)
     enddo
  enddo
  
! compute correct state for rarefaction fan by linear interpolation

  do i = 5, numIntCells5
     scrch1(i) = max (wes(i) - westar(i), wes(i) + westar(i), hy_smallu)
     scrch1(i) =     (wes(i) + westar(i)) / scrch1(i)
     
     scrch1(i) = 0.5e0 * (1.e0 + scrch1(i))
     scrch2(i) =          1.e0 - scrch1(i)
     
     rhoav(i)  = scrch1(i) * rhostr(i) + scrch2(i) * rhos (i)
     uav  (i)  = scrch1(i) * ustar(i)  + scrch2(i) * us(i)
     utav (i)  = uts(i)
     uttav(i)  = utts(i)
     pav   (i) = scrch1(i) * pstar(i)  + scrch2(i) * ps(i)
     gameav(i) = scrch1(i) * gmstar(i) + scrch2(i) * games(i)
  enddo
  
  do i = 5, numIntCells5
     if (westar(i) .ge. 0.e0) then
        rhoav(i)  = rhostr(i)
        uav(i)    = ustar(i)
        pav(i)    = pstar(i)
        gameav(i) = gmstar(i)
     endif
     
     if (wes(i) .lt. 0.e0) then
        rhoav(i)  = rhos(i)
        uav(i)    = us(i)
        pav(i)    = ps(i)
        gameav(i) = games(i)
     endif
     
     urell(i) = uav(i) - ugrdl(i)
  enddo

  return
end subroutine rieman
  
