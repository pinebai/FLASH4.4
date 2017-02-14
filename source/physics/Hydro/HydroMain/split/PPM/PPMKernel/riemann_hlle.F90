!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/riemann_hlle
!!
!! NAME
!!
!!  riemann_hlle
!!
!!
!! SYNOPSIS
!!
!!  call riemann_hlle ( real(in) :: rho_left,  
!!                real(in) :: rho_right, 
!!                real(in) :: u_left,
!!                real(in) :: u_right, 
!!                real(in) :: ut_left, 
!!                real(in) :: ut_right, 
!!                real(in) :: utt_left,
!!                real(in) :: utt_right, 
!!               real(in) :: p_left,
!!               real(in) :: p_right, 
!!               real(in) :: e_left,
!!               real(in) :: e_right, 
!!               real(in) :: ei_left,
!!               real(in) :: ei_right,
!!               real(in) :: game_left,
!!               real(in) :: game_right,
!!               real(in) :: gamc_left,
!!               real(in) :: gamc_right,
!!               real(in) :: grav_left, 
!!               real(in) :: grav_right,
!!               real(in) :: dt,
!!               real(inout) :: rho_flux, 
!!               real(inout) :: u_flux, 
!!               real(inout) ::ut_flux, 
!!               real(inout) :: utt_flux, 
!!               real(inout) :: e_flux, 
!!               real(inout) :: ei_flux)
!!
!! DESCRIPTION
!!
!!  An approximate-flux Riemann solver that implements the HLLE 
!!  (Harten-Lax-Van Leer-Einfeldt) algorithm.  This is used in shocks
!!  to avoid carbuncle pheonomenon or odd-even decoupling.  In contrast
!!  to other Riemann solvers, which find an approximation to the pressure
!!  and density at the contact discontinuity, HLLE approximates the 
!!  largest and smallest signal velocities (Einfeldt 1998).
!!
!!  Given the states on the left and right of the interface, compute the
!!  fluxes through the interface.  
!!
!!  see Einfeldt, B., Munz, C.D., Roe, P.L., Sjoegren, B., 1991, JCP, 92, 273
!!
!!      Einfeldt, B., 1988, SIAM JNA, 25, 294
!!
!!      Roe, P.L., 1981, JCP, 43, 357
!!
!!      Einfeldt, B., 1988, in "Shock Tubes and Waves", ed. H. Groenig, 
!!         Proceedings of the Sixteenth Int.Symp. on Shock Tubes and Waves, 
!!         Aachen July 26-31, 1987, p. 671
!!
!!   (equation numbers refer to the first reference, unless otherwise noted)
!!
!!
!! NOTES
!! 
!!  -- This operates on the zone average values to the left and the
!!     right of the interface, *NOT* the reconstructed interface values
!!     as returned from states.  Since we are only applying this inside 
!!     shocks, this is not really a problem.
!!
!!  -- The fluxes from the previous Riemann solver should be passed in, 
!!     and if the right conditions are met, they will be replaced with the 
!!     appropriate HLLE fluxes.  There is no pressure 'flux', so pav from 
!!     the regular Riemann solver will be used.
!!
!!  -- We ignore the grid velocity in here.
!!
!!  -- This code is adapted from AMRA, http://www.camk.edu.pl/~tomek/AMRA/
!!
!! ARGUMENTS
!!
!! rho_left:  
!! rho_right: 
!! u_left:
!! u_right: 
!! ut_left: 
!! ut_right: 
!! utt_left:
!! utt_right: 
!! p_left:
!! p_right: 
!! e_left:
!! e_right: 
!! ei_left:
!! ei_right:
!! game_left:
!! game_right:
!! gamc_left:
!! gamc_right:
!! grav_left: 
!! grav_right:
!! dt:
!! rho_flux: 
!! u_flux: 
!! ut_flux: 
!! utt_flux: 
!! e_flux: 
!! ei_flux:
!!
!!
!!  
!!***

subroutine riemann_hlle(rho_left,  rho_right, &
                        u_left,    u_right, &
                        ut_left,   ut_right, &
                        utt_left,  utt_right, &
                        p_left,    p_right, &
                        e_left,    e_right, &
                        ei_left,   ei_right, &
                        game_left, game_right, &
                        gamc_left, gamc_right, &
                        grav_left, grav_right, &
                        dt, &
                        rho_flux, u_flux, ut_flux, utt_flux, e_flux, ei_flux)

  use Hydro_data, ONLY : hy_small

  implicit none

  real, INTENT(in) :: rho_left, rho_right
  real, INTENT(in) :: u_left, u_right
  real, INTENT(in) :: ut_left, ut_right
  real, INTENT(in) :: utt_left, utt_right
  real, INTENT(in) :: p_left, p_right
  real, INTENT(in) :: e_left, e_right
  real, INTENT(in) :: ei_left, ei_right
  real, INTENT(in) :: game_left, game_right
  real, INTENT(in) :: gamc_left, gamc_right
  real, INTENT(in) :: grav_left, grav_right

  real, INTENT(in) :: dt
  real, INTENT(inout) :: rho_flux, u_flux, ut_flux, utt_flux, e_flux, ei_flux

  real :: rls, rrs, rssi
  real :: rvm

  real :: as_left, as_right, am

  real :: a1, a4
  real :: chalf, uhalf

  real :: urell_l, urell_r

  real :: bl, br

  real :: betag

  real :: bm, bp, bdif, bdifi

  real :: f1l, f2l, f3l, f4l, f5l, f6l
  real :: f1r, f2r, f3r, f4r, f5r, f6r

  real :: u1g_l, u1g_r

! which approximation to the speeds do we use -- eventually, we'll make this a 
! parameter
  integer, parameter :: metcav = 2


! compute the sqrt of the density on the left and right of the interface
! and a common factor (see Eq. 5.3a)

  rls  = sqrt(rho_left)
  rrs  = sqrt(rho_right)
  rssi = 1.0/(rls+rrs)


! compute the average velocity (see Eq. 5.3b)

  rvm = (rls*u_left + rrs*u_right)*rssi


! compute the sound speeds for the left and right states
  as_left = sqrt(gamc_left*p_left/rho_left)
  as_right = sqrt(gamc_right*p_right/rho_right)


! compute the average sound speed.  (see Einfeldt, SIAM JNA, Eq. 5.7)

  am = sqrt((rls*as_left**2 + rrs*as_right**2)*rssi &
       + 0.5*rls*rrs*rssi**2 *(u_right - u_left)**2 )


! Roe's eigenvalues
  
  a1 = rvm - am
  a4 = rvm + am


! numerical signal velocities -- there are three choices here, outlined
! in Einfeldt's paper.  They have differing dissipation properties

  if ( metcav == 0 ) then

! (4.7a), (4.7b) [most dissipative]

     chalf = 0.5*(a4 - a1)
     uhalf = 0.5*(a4 + a1)

     urell_l = u_left
     urell_r = u_right

     br = max(abs(uhalf)    + chalf, &
              abs(urell_l) + as_left, &
              abs(urell_r) + as_right)

     bl = -br


  elseif ( metcav == 1 ) then

! (4.5a), (4.5b) [standard]

     bl = min( a1, u_left - as_left )
     br = max( a4, u_right + as_right )


  elseif ( metcav == 2 ) then

! (4.9b), (4.10a), (4.10b) [sharpest]

     betag = 0.5*( game_left + game_right )
     betag = sqrt( 0.5*(betag-1.e0)/betag )

     bl = min(a1, u_left - betag*as_left)
     br = max(a4, u_right + betag*as_right)

  end if


! see comment to (4.4b)

  bm = min( 0.e0, bl )
  bp = max( 0.e0, br )

  bdif = bp - bm

! only apply the HLLE fluxes if the difference in the signal speeds is
! large enough to notice

  if ( abs(bdif) > hy_small*max(abs(bp),abs(bm)) ) then

     bdifi = 1.e0/bdif

     u1g_l = u_left  + 0.5*dt*grav_left
     u1g_r = u_right + 0.5*dt*grav_right

     
! we are ignoring the grid velocity, so the relative velocity is just the
! velocity

     urell_l = u1g_l
     urell_r = u1g_r


!----------------------------------------------------------------------------
! left base fluxes
!
!   these are just the fluxes (Eq. 1.2) from the Euler equation written as 
!   U_t + F_x = 0.  We evaluate the fluxes, F(U), using the values in the
!   zone to the left of the interface.  The HLLEing will come later.
!
!   notice that pressure does not appear in momentum flux as it is 
!   incorporated in momentum equation
!----------------------------------------------------------------------------

! density
     f1l = rho_left*urell_l

! momentum
     f2l = f1l*u1g_l

#if N_DIM >= 2
     f3l = f1l*ut_left
#else
     f3l = 0.e0
#endif

#if N_DIM == 3
     f4l = f1l*utt_left
#else
     f4l = 0.e0
#endif
     
! total energy
     f5l = f1l*e_left + u1g_l*p_left

! internal energy -- note, our internal energy flux has a pressure term
     f6l = f1l*ei_left + u1g_l*p_left


!----------------------------------------------------------------------------
! right fluxes
!
!   these are just the fluxes (Eq. 1.2) from the Euler equation written as 
!   U_t + F_x = 0.  We evaluate the fluxes, F(U), using the values in the
!   zone to the right of the interface.  The HLLEing will come later.
!----------------------------------------------------------------------------

! density
     f1r = rho_right*urell_r

! momentum
     f2r = f1r*u1g_r

#if N_DIM >= 2
     f3r = f1r*ut_right
#else
     f3r = 0.e0
#endif

#if N_DIM == 3
     f4r = f1r*utt_right
#else
     f4r = 0.e0
#endif

! total energy
     f5r = f1r*e_right + u1g_r*p_right

! internal energy -- note, our internal energy flux has a pressure term
     f6r = f1r*ei_right + u1g_r*p_right


!----------------------------------------------------------------------------
! HLLE fluxes (4.4b)
!----------------------------------------------------------------------------

     rho_flux = (bp*f1l - bm*f1r + bp*bm*(rho_right - rho_left))*bdifi

     u_flux = (bp*f2l - bm*f2r + &
          bp*bm*(rho_right*u_right - rho_left*u_left))*bdifi


#if N_DIM >= 2
     ut_flux = (bp*f3l - bm*f3r + &
          bp*bm*(rho_right*ut_right - rho_left*ut_left))*bdifi
#else
     ut_flux = 0.e0
#endif

#if N_DIM == 3
     utt_flux = (bp*f4l - bm*f4r + &
          bp*bm*(rho_right*utt_right - rho_left*utt_left))*bdifi
#else
     utt_flux = 0.e0
#endif

     e_flux = (bp*f5l - bm*f5r + &
          bp*bm*(rho_right*e_right - rho_left*e_left))*bdifi

     ei_flux = (bp*f6l - bm*f6r + &
          bp*bm*(rho_right*ei_right - rho_left*ei_left))*bdifi

  end if

  return
end subroutine riemann_hlle
