!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_states
!!
!! NAME
!!
!!   hy_rhd_states
!!
!! SYNOPSIS
!!   hy_rhd_states( real   (IN) :: Vc(:),
!!                  real  (OUT) :: Vm(:),
!!                  real  (OUT) :: Vp(:),
!!                  real   (IN) :: grav(:),
!!                  real   (IN) :: dt,
!!                  real   (IN) :: r(:),
!!                  real   (IN) :: dr(:),
!!                  real   (IN) :: dvol(:),
!!                  integer(IN) :: n,
!!                  integer(IN) :: dir)
!!
!!
!! DESCRIPTION
!!
!!  Advance rhd equations in primitive form using characteristic
!!  tracing with upwind limiting.
!!
!!  Node arrangements:
!!
!!                     i-1/2          i          i+1/2
!!                          |m--------c--------p|
!!
!!                Vm(i) = lim     V(x)
!!                     x->x(i-1/2)+
!!
!!                Vp(i) = lim     V(x)
!!                     x->x(i+1/2)-
!!
!!
!!  ARGUMENTS
!!
!!    Vc   -  1-D array of primitive values, center value
!!    Vm   -  1-D array of primitive values, left interface limit
!!    Vp   -  1-D array of primitive values, right interface limit,
!!    grav -  gravity (not used for now)
!!    dt   -  time step
!!    r    -  1-D array with cell centers
!!    dr   -  1-D array with mesh spacing
!!    dvol -  1-D array with volume spacing
!!    n    -  number of points
!!    dir  -  sweep direction
!!      
!!
!!***

subroutine hy_rhd_states(Vc, Vm, Vp, grav, dt, r, dr, dvol, n, dir)
  
  use Hydro_data, ONLY : hy_meshGeom

  use Grid_data,  ONLY : gr_iguard,gr_jguard,gr_kguard

  use hy_rhd_interface, ONLY:  hy_rhd_reconstruct, &
                               hy_rhd_sources,     &
                               hy_rhd_soundSpeed2, &
                               hy_rhd_checkBoundaryValues

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "RHD.h"

  !! Argument list ---------------------------------------------
  integer, INTENT(IN)                  :: n, dir
  real, DIMENSION(NUNK_VARS,n), INTENT(IN)  :: Vc
  real, DIMENSION(NUNK_VARS,n), INTENT(OUT) :: Vm, Vp
  real, DIMENSION(n), INTENT(IN)       :: grav, r, dr, dvol
  real, INTENT(IN)                     :: dt
  !! -----------------------------------------------------------
  
  integer                        :: i, k, nv, ifail
  integer, PARAMETER             :: rhd_char = 5

  real, DIMENSION(n)             :: h, cs2
  real, DIMENSION(rhd_char)      :: alpha_p, alpha_m, scrh_p, scrh_m, lambda
  real, DIMENSION(rhd_char)      :: drp, drm, xip, xim
  real, DIMENSION(NUNK_VARS)          :: dV, V6, vp_ref, vm_ref, dvp, dvm
  real, DIMENSION(NUNK_VARS,rhd_char) :: vp_av, vm_av
  real, DIMENSION(NUNK_VARS,n)        :: src
  real     ::  cs, lorentz, eta, lorentz2_1d, scrh, D
  real     ::  a, b, c, r2, r5, delta_2
  real     ::  vx, vy, vz, volm, volp

  integer  :: VELN_VAR,VELT1_VAR,VELT2_VAR
  integer  :: nGuard,sweepBegin,sweepEnd


  select case(dir)
  case (SWEEP_X)
     VELN_VAR = VELX_VAR
     VELT1_VAR= VELY_VAR
     VELT2_VAR= VELZ_VAR
     nGuard   = gr_iguard
  case (SWEEP_Y)
     VELN_VAR = VELY_VAR
     VELT1_VAR= VELZ_VAR
     VELT2_VAR= VELX_VAR
     nGuard   = gr_jguard
  case (SWEEP_Z)
     VELN_VAR = VELZ_VAR
     VELT1_VAR= VELX_VAR
     VELT2_VAR= VELY_VAR
     nGuard   = gr_kguard
  end select

  ! Set the required array length
  sweepBegin = nGuard
  sweepEnd   = n-nGuard+1

  ! Check boundary values  
  call hy_rhd_checkBoundaryValues(Vc,n,dir)
  
  ! ----------------------------------
  !  Uncomment the next three lines 
  !  if you want flat reconstruction
  ! ----------------------------------                 
  
  ! Vp = Vc
  ! Vm = Vc
  ! return
  
  call hy_rhd_reconstruct(Vc, Vm, Vp, r, dr, dvol, n, dir)  ! get left and right states
  call hy_rhd_enthalpy(Vc, h, 1, n, n)                      ! get enthalpy
  call hy_rhd_soundSpeed2(Vc, h, cs2, 1, n , n)             ! get sound speed
  call hy_rhd_sources(Vc, h, cs2, src, r, dr, sweepBegin, sweepEnd, n, dir, HY_PRIMITIVE)

  do i = sweepBegin, sweepEnd
     dV(:) = Vp(:,i) - Vm(:,i)
     V6(:) = 2.e0*Vc(:,i) - (Vp(:,i) + Vm(:,i))
     
     ! ------------------------------------------------
     !  Compute eigenvalues (lambda), and define
     ! some useful quantities
     ! ------------------------------------------------
     
     vx = Vc(VELN_VAR,i)        ! normal velocity
     vy = Vc(VELT1_VAR,i)       ! 1st tangential component
     vz = Vc(VELT2_VAR,i)       ! 2nd tangential component

     eta     = vy*vy + vz*vz
     a       = eta + vx*vx
     eta     = sqrt(1.e0 - vx*vx - cs2(i)*eta)
     lorentz = 1.e0/sqrt(1.e0 - a)
     delta_2 = 1.e0/(1.e0 - a*cs2(i))
     cs      = sqrt(cs2(i))

     lambda(1)   = (vx*(1.e0 - cs2(i)) - cs/lorentz*eta)*delta_2
     lambda(2:4) =  vx
     lambda(5)   = (vx*(1.e0 - cs2(i)) + cs/lorentz*eta)*delta_2
     
     lorentz2_1d = 1.e0/(1.e0 - vx*vx)
     D  = Vc(DENS_VAR,i)*lorentz; !  Lab - Density 
     r2 = D/(cs*eta)
     r5 = 1.e0/(h(i)*cs2(i))
     c  = lorentz*eta*vx
     a  = c - cs
     b  = c + cs
     c  = lorentz2_1d/(lorentz*D)
     
     ! --------------------------------------------------
     !  For each characteristic "k", find the domain of 
     !  dependence for each interface (scrh_p for i+1/2,
     !  scrh_m for i-1/2);
     ! --------------------------------------------------
     
     do k = 1, 5 
        drp(k)  = dt*max(0.e0,  lambda(k))
        drm(k)  = dt*max(0.e0, -lambda(k))
     end do

     if ((hy_meshGeom == SPHERICAL) .and.( dir == SWEEP_X)) then
        
        volm = (r(i) - 0.5*dr(i))**3.0/3.e0
        volp = volm + dvol(i)
        do k = 1, 5
           xip(k) = (r(i) + drp(k))**3.e0/3.e0
           xim(k) = (r(i) - drm(k))**3.e0/3.e0
           
           scrh_p(k) =  (volp - xip(k))/dvol(i)
           scrh_m(k) = -(volm - xim(k))/dvol(i)
        end do
        
     else if ((hy_meshGeom == CYLINDRICAL) .and. (dir == SWEEP_X)) then
        
        volm = 0.5*(r(i) - 0.5*dr(i))*(r(i) - 0.5*dr(i))
        volp = volm + dvol(i)
        do k = 1, 5
           xip(k) = 0.5*(r(i) + drp(k))*(r(i) + drp(k))
           xim(k) = 0.5*(r(i) - drm(k))*(r(i) - drm(k))
           
           scrh_p(k) =  (volp - xip(k))/dvol(i)
           scrh_m(k) = -(volm - xim(k))/dvol(i)
        end do
        
     else   ! -- CARTESIAN is always nicer... --
        
        scrh_p(:) = drp(:)/dr(i)
        scrh_m(:) = drm(:)/dr(i)
        
     end if
     
     ! --------------------------------------------------
     !  Calculate the average of V over the
     !  domain of dependence lying to the left (for
     !  Vp) or to the right (for Vm) of the interface.
     ! --------------------------------------------------
     
     do k = 1, 5 
        if (lambda(k) > 0.e0) then
           vp_av(:,k) = Vp(:,i) - 0.5d0*scrh_p(k)*   &
                (dV(:) - (3.e0 - 2.e0*scrh_p(k))*V6(:))
           vm_av(:,k) = Vm(:,i)
           
        else
           vp_av(:,k) = Vp(:,i) 
           vm_av(:,k) = Vm(:,i) + 0.5d0*scrh_m(k)*   &
                (dV(:) + (3.e0 - 2.e0*scrh_m(k))*V6(:))
        end if
     end do
     
     !  ----  define the reference states  ---- 
     
     vp_ref(:) = vp_av(:, rhd_char)
     vm_ref(:) = vm_av(:, 1)
     
     ! ------------------------------------------------
     !  Use upwind limiting to select only those
     !  characteristics which contribute to the 
     !  effective left and right states. 
     !  alpha_p(k) and alpha_m(k) are projection
     !  coefficients on the left eigenvectors
     !  of the primitive relativistic equations, i.e.
     !
     !   alpha+(k) = l(k) . dV+(k)
     !   alpha-(k) = l(k) . dV-(k)
     ! ------------------------------------------------
     
     if (vx .gt. 0.e0) then       ! contribution to i+1/2  (p)
        
        if (lambda(1) .gt. 0.e0) then 
           dV(VELN_VAR)    = vp_ref(VELN_VAR) - vp_av(VELN_VAR, 1)
           dV(PRES_VAR)    = vp_ref(PRES_VAR) - vp_av(PRES_VAR, 1)
           alpha_p(1) = 0.5d0*(-r2*dV(VELN_VAR) + r5*dV(PRES_VAR))                         
           alpha_m(1) = 0.e0
        else
           dV(VELN_VAR)    = vm_ref(VELN_VAR) - vm_av(VELN_VAR, 1)
           dV(PRES_VAR)    = vm_ref(PRES_VAR) - vm_av(PRES_VAR, 1)
           alpha_p(1) = 0.e0                         
           alpha_m(1) = 0.5d0*(-r2*dV(VELN_VAR) + r5*dV(PRES_VAR))
        end if
        
        dV(DENS_VAR)    = vp_ref(DENS_VAR) - vp_av(DENS_VAR, 2)
        dV(PRES_VAR)    = vp_ref(PRES_VAR) - vp_av(PRES_VAR, 2)
        alpha_p(2) = dV(DENS_VAR) - r5*dV(PRES_VAR)                                    
        alpha_m(2) = 0.e0
        
        dV(VELN_VAR)    = vp_ref(VELN_VAR)  - vp_av(VELN_VAR, 3)
        dV(VELT1_VAR)   = vp_ref(VELT1_VAR) - vp_av(VELT1_VAR, 3)
        dV(PRES_VAR)    = vp_ref(PRES_VAR) - vp_av(PRES_VAR, 3)
        alpha_p(3) = vx*vy*lorentz2_1d*dV(VELN_VAR) + dV(VELT1_VAR) + vy*c/h(i)*dV(PRES_VAR)   
        alpha_m(3) = 0.e0
        
        dV(VELN_VAR)    = vp_ref(VELN_VAR)  - vp_av(VELN_VAR, 4)
        dV(VELT2_VAR)   = vp_ref(VELT2_VAR) - vp_av(VELT2_VAR, 4)
        dV(PRES_VAR)    = vp_ref(PRES_VAR) - vp_av(PRES_VAR, 4)
        alpha_p(4) = vx*vz*lorentz2_1d*dV(VELN_VAR) + dV(VELT2_VAR) + vz*c/h(i)*dV(PRES_VAR)        
        alpha_m(4) = 0.e0
        
        dV(VELN_VAR)    = vp_ref(VELN_VAR) - vp_av(VELN_VAR, 5)
        dV(PRES_VAR)    = vp_ref(PRES_VAR) - vp_av(PRES_VAR, 5)
        alpha_p(5) = 0.5d0*(r2*dV(VELN_VAR) + r5*dV(PRES_VAR))                                  
        alpha_m(5) = 0.e0
        
     else                           ! contribution to i-1/2  (m)
        
        dV(VELN_VAR)    = vm_ref(VELN_VAR) - vm_av(VELN_VAR, 1)
        dV(PRES_VAR)    = vm_ref(PRES_VAR) - vm_av(PRES_VAR, 1)
        alpha_p(1) = 0.e0                         
        alpha_m(1) = 0.5d0*(-r2*dV(VELN_VAR) + r5*dV(PRES_VAR))
        
        dV(DENS_VAR)    = vm_ref(DENS_VAR) - vm_av(DENS_VAR, 2)
        dV(PRES_VAR)    = vm_ref(PRES_VAR) - vm_av(PRES_VAR, 2)
        alpha_p(2) = 0.e0
        alpha_m(2) = dV(DENS_VAR) - r5*dV(PRES_VAR)                                    
        
        dV(VELN_VAR)    = vm_ref(VELN_VAR)  - vm_av(VELN_VAR, 3)
        dV(VELT1_VAR)   = vm_ref(VELT1_VAR) - vm_av(VELT1_VAR, 3)
        dV(PRES_VAR)    = vm_ref(PRES_VAR)  - vm_av(PRES_VAR, 3)
        alpha_p(3) = 0.e0
        alpha_m(3) = vx*vy*lorentz2_1d*dV(VELN_VAR) + dV(VELT1_VAR) + vy*c/h(i)*dV(PRES_VAR)   
        
        dV(VELN_VAR)    = vm_ref(VELN_VAR)  - vm_av(VELN_VAR, 4)
        dV(VELT2_VAR)   = vm_ref(VELT2_VAR) - vm_av(VELT2_VAR, 4)
        dV(PRES_VAR)    = vm_ref(PRES_VAR)  - vm_av(PRES_VAR, 4)
        alpha_p(4) = 0.e0
        alpha_m(4) = vx*vz*lorentz2_1d*dV(VELN_VAR) + dV(VELT2_VAR) + vz*c/h(i)*dV(PRES_VAR)        
        
        if (lambda(5) .lt. 0.e0) then 
           dV(VELN_VAR)    = vm_ref(VELN_VAR) - vm_av(VELN_VAR, 5)
           dV(PRES_VAR)    = vm_ref(PRES_VAR) - vm_av(PRES_VAR, 5)
           alpha_p(5) = 0.e0
           alpha_m(5) = 0.5d0*(r2*dV(VELN_VAR) + r5*dV(PRES_VAR))                                  
        else
           dV(VELN_VAR)    = vp_ref(VELN_VAR) - vp_av(VELN_VAR, 5)
           dV(PRES_VAR)    = vp_ref(PRES_VAR) - vp_av(PRES_VAR, 5)
           alpha_p(5) = 0.5d0*(r2*dV(VELN_VAR) + r5*dV(PRES_VAR))                                  
           alpha_m(5) = 0.e0
        end if
        
     end if
     
     r2 = 1.e0/r2;
     r5 = 1.e0/r5;
     
     ! ---------------------------------------------------
     !  sum contributions of alpha in the right
     !  eigenvector expansion, i.e.
     !
     !   dvp = \sum_(k, \lambda_k > 0)  alpha_p(k) * r_k
     !   dvm = \sum_(k, \lambda_k < 0)  alpha_m(k) * r_k
     !
     ! ---------------------------------------------------
     
     dvp(DENS_VAR)  = alpha_p(1) + alpha_p(2) + alpha_p(5)
     dvm(DENS_VAR)  = alpha_m(1) + alpha_m(2) + alpha_m(5)
     dvp(VELN_VAR)  = (-alpha_p(1) + alpha_p(5))*r2
     dvm(VELN_VAR)  = (-alpha_m(1) + alpha_m(5))*r2
     dvp(VELT1_VAR) = cs*c*vy*(alpha_p(1)*a - alpha_p(5)*b) + alpha_p(3)
     dvm(VELT1_VAR) = cs*c*vy*(alpha_m(1)*a - alpha_m(5)*b) + alpha_m(3)
     dvp(VELT2_VAR) = cs*c*vz*(alpha_p(1)*a - alpha_p(5)*b) + alpha_p(4)
     dvm(VELT2_VAR) = cs*c*vz*(alpha_m(1)*a - alpha_m(5)*b) + alpha_m(4)
     dvp(PRES_VAR)  = (alpha_p(1) + alpha_p(5))*r5
     dvm(PRES_VAR)  = (alpha_m(1) + alpha_m(5))*r5
     
     ! ----------------------------------------
     !   compute time-centered states for 
     !   Riemann solver
     ! ----------------------------------------
     Vp(DENS_VAR,i)  = vp_ref(DENS_VAR)  - dvp(DENS_VAR)  + 0.5d0*dt*src(DENS_VAR,i)
     Vp(VELN_VAR,i)  = vp_ref(VELN_VAR)  - dvp(VELN_VAR)  + 0.5d0*dt*src(VELN_VAR,i)
     Vp(VELT1_VAR,i) = vp_ref(VELT1_VAR) - dvp(VELT1_VAR) + 0.5d0*dt*src(VELT1_VAR,i)
     Vp(VELT2_VAR,i) = vp_ref(VELT2_VAR) - dvp(VELT2_VAR) + 0.5d0*dt*src(VELT2_VAR,i)
     Vp(PRES_VAR,i)  = vp_ref(PRES_VAR)  - dvp(PRES_VAR)  + 0.5d0*dt*src(PRES_VAR,i)

     Vm(DENS_VAR,i)  = vm_ref(DENS_VAR)  - dvm(DENS_VAR)  + 0.5d0*dt*src(DENS_VAR,i)
     Vm(VELN_VAR,i)  = vm_ref(VELN_VAR)  - dvm(VELN_VAR)  + 0.5d0*dt*src(VELN_VAR,i)
     Vm(VELT1_VAR,i) = vm_ref(VELT1_VAR) - dvm(VELT1_VAR) + 0.5d0*dt*src(VELT1_VAR,i)
     Vm(VELT2_VAR,i) = vm_ref(VELT2_VAR) - dvm(VELT2_VAR) + 0.5d0*dt*src(VELT2_VAR,i)
     Vm(PRES_VAR,i)  = vm_ref(PRES_VAR)  - dvm(PRES_VAR)  + 0.5d0*dt*src(PRES_VAR,i)

   
     ! ----------------------------------------
     !           Check states:
     ! 
     !  revert to 1st order when at least one 
     !  of the following circumstances is 
     !  encountered:
     !
     !   - the total velocity exceeds 1
     !   - the density is negative
     !   - the pressure is negative 
     ! ----------------------------------------
     
     a = Vp(VELX_VAR,i)*Vp(VELX_VAR,i) + Vp(VELY_VAR,i)*Vp(VELY_VAR,i) +  Vp(VELZ_VAR,i)*Vp(VELZ_VAR,i)
     b = Vm(VELX_VAR,i)*Vm(VELX_VAR,i) + Vm(VELY_VAR,i)*Vm(VELY_VAR,i) +  Vm(VELZ_VAR,i)*Vm(VELZ_VAR,i)
     
     ifail = 0
     if (a .ge. 1.e0) then
        !print *," ! vyp > 1 in states "
        ifail = 1
     end if
     if (b .ge. 1.e0) then
        !print *," ! vym > 1 in states "
        ifail = 1
     end if
     
     if (Vp(PRES_VAR,i) .le. 0.e0) then
        !print *," Negative pressure + in states"
        ifail = 1
     end if
     if (Vm(PRES_VAR,i) .le. 0.e0) then
        !print *," Negative pressure - in states"
        ifail = 1
     end if
     if (Vp(DENS_VAR,i) .le. 0.e0) then
        !print *," Negative pressure + in states"
        ifail = 1
     end if
     if (Vm(DENS_VAR,i) .le. 0.e0) then
        !print *," Negative pressure - in states"
        ifail = 1
     end if
     if (ifail == 1) then 
        Vp = Vc
        Vm = Vc
     end if
  end do
end subroutine hy_rhd_states
