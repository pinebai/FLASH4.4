!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_reconstruct
!!
!! NAME
!!
!!  hy_rhd_reconstruct
!!
!! SYNOPSIS
!!  hy_rhd_reconstruct( real    (IN) :: a(:),
!!                   real   (OUT) :: am(:),
!!                   real   (OUT) :: ap(:),
!!                   real    (IN) :: r(:),
!!                   real    (IN) :: dr(:),
!!                   real    (IN) :: dvol(:),
!!                   integer (IN) :: n,
!!                   integer (IN) :: dir)
!!
!!
!! DESCRIPTION
!!
!!   Provide left and right interpolated values using
!!   either  linear or PPM reconstruction
!!
!!   Node arrangements:
!!
!!                         i-1/2      i        i+1/2
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
!!    a   -    volume averages of primitive quantities (?)
!!    am  -    1-D array of primitive values, left interface limit
!!    ap  -    1-D array of primitive values, right interface limit,
!!    r    -   1-D array with cell centers
!!    dr   -   1-D array with mesh spacing
!!    dvol -   1-D array with volume spacing
!!    n    -   number of points
!!    dir  -   sweep direction
!!
!!
!!
!!***


subroutine hy_rhd_reconstruct (a, am, ap, r, dr, dvol, n, dir)

  use Hydro_data, ONLY : hy_meshGeom, hy_reconType
  implicit none

#include "constants.h"
#include "Flash.h"
#include "RHD.h"

  !! Argument list ---------------------------------------
  integer, INTENT(in)                  :: n, dir
  real, DIMENSION(NUNK_VARS,n), INTENT(in)  :: a
  real, DIMENSION(NUNK_VARS,n), INTENT(out) :: am, ap
  real, DIMENSION(n),      INTENT(in)  :: r, dr, dvol
  !! -----------------------------------------------------

  real                    :: ald, ard, eta_t, etaj
  real, DIMENSION(NUNK_VARS)   :: scrh1, scrh2, scrh3, scrh4
  real, DIMENSION(NUNK_VARS,n) :: da, da_lim
  real                    :: dam, dap, dac

  real, DIMENSION(n)    :: delta2, f_t, rho2,  fj
  real, DIMENSION(n)    :: dp, d2p, min_p
  real, DIMENSION(n)    :: c1, c2, c3
  real, DIMENSION(NUNK_VARS) :: steepness

  real, PARAMETER       :: EPS2 = 1.e0
  real, PARAMETER       :: OME1 = 0.52d0
  real, PARAMETER       :: OME2 = 10.e0

  real                   :: beta, dvn, fltn, v2l, v2r, v2_max, v2_min
  logical, DIMENSION(n)  :: l_shock
  integer,save           :: j, nv,  sj, ivn

  integer, PARAMETER :: RECON_LINEAR=1,RECON_PARABOLIC=2

  ! -------------------------------------------------------
  !            Use Piecewise Linear Interpolants
  !            when use_plm = .TRUE.
  ! ------------------------------------------------------- 

  ivn = VELX_VAR + dir - IAXIS

  if(hy_reconType == RECON_LINEAR) then
     do j = 2,n-1
        do nv = 1, NUNK_VARS

           dap = a(nv, j + 1) - a(nv, j)
           dam = a(nv, j)     - a(nv, j - 1)
           dac = a(nv, j + 1) - a(nv, j - 1)

           if (dam*dap .gt. 0.e0) then 
              !  -- MC limiter --
              !        scrh1(nv) = min(2.e0*min(abs(dam), abs(dap)), 0.5d0*abs(dac))
              !        scrh1(nv) = scrh1(nv)*sign(1.e0, dac)

              ! -- Van Leer Limiter --
              !        scrh1(nv) = 2.e0*dam*dap/(dam + dap + 1.e-12)

              ! -- Minmod Limiter
              scrh1(nv) = min(abs(dam), abs(dap))
              scrh1(nv) = scrh1(nv)*sign(1.e0, dap)
              
           else
              scrh1(nv) = 0.e0
           endif
        end do

        am(:,j) = a(:,j) - 0.5d0*scrh1(:)
        ap(:,j) = a(:,j) + 0.5d0*scrh1(:)
     end do

  end if

  if(hy_reconType == RECON_PARABOLIC) then
     steepness = 2.e0

     ! ---- PPM interpolation coefficients ----
     if (hy_meshGeom .eq. SPHERICAL .and. dir .eq. IAXIS) then

        do j = 2, n - 2
           beta  = 1.0/(dvol(j-1) + dvol(j) + dvol(j+1) + dvol(j+2))
           c2(j) = beta*dvol(j+1)*(dvol(j+1) + dvol(j+2))/   &
                (dvol(j) + 2.e0*dvol(j+1))
           c3(j) = beta*dvol(j)*(dvol(j-1) + dvol(j))/   &
                (2.0*dvol(j) + dvol(j+1))
           c1(j) = dvol(j) + 2.0*(dvol(j+1)*c3(j) - dvol(j)*c2(j))
           c1(j) = c1(j)/(dvol(j) + dvol(j+1))
        end do

     else if (hy_meshGeom .eq. CYLINDRICAL .and. dir .eq. IAXIS) then
        
        do j = 2, n - 2
           beta  = 1.0/(dvol(j-1) + dvol(j) + dvol(j+1) + dvol(j+2))
           c2(j) = beta*dvol(j+1)*(dvol(j+1) + dvol(j+2))/   &
                   (dvol(j) + 2.e0*dvol(j+1))
           c3(j) = beta*dvol(j)*(dvol(j-1) + dvol(j))/   &
                   (2.0*dvol(j) + dvol(j+1))
           c1(j) = dvol(j) + 2.0*(dvol(j+1)*c3(j) - dvol(j)*c2(j))
           c1(j) = c1(j)/(dvol(j) + dvol(j+1))
        end do
        
     else   ! in CARTESIAN geometry coefficients reduce to this
        
        c1(:) = 0.5d0
        c2(:) = 1.e0/6.e0
        c3(:) = 1.e0/6.e0
        
     end if
     
! ---------------------------------------------------------
!           MAIN LOOP ON VARIABLES
!  ---------------------------------------------------------  

     do j = 2, n
        da(:,j) = a(:,j) - a(:,j - 1)
     end do

! ---------------------------------------------------------
!                           step #1
!  --------------------------------------------------------- 

     do j = 2, n - 1
        scrh1(:) = da(:,j + 1)*da(:,j)
        do nv = 1,NUNK_VARS 
           if (scrh1(nv) .gt. 0.e0) then 
              scrh2(nv)    = steepness(nv)*min(abs(da(nv, j)), abs(da(nv, j + 1)))
              scrh2(nv)    = min(0.5e0*abs(da(nv, j + 1) + da(nv, j)), scrh2(nv))
              da_lim(nv,j) = scrh2(nv)*sign(1.e0, da(nv,j))
           else
              da_lim(nv,j) = 0.e0
           end if
        end do
     end do
     
     do j = 2, n - 2 
        scrh1(:) = a(:,j) + c1(j)*da(:,j + 1) &
             + c2(j)*da_lim(:,j) - c3(j)*da_lim(:,j + 1)
        ap(:, j)     = scrh1(:) 
        am(:, j + 1) = scrh1(:)
     end do

! ---------------------------------------------------
!          step # 2 :    CONTACT STEEPENING
! --------------------------------------------------- 

!         no contact steepening for now

! -------------------------------------------------------
!                 S T E P     # 3: Parabolic Limiter
!  ------------------------------------------------------- 

     do j = 3, n - 3
        scrh1(:) = (ap(:,j) - a(:,j))*(a(:,j) - am(:,j))
        do nv = 1,NUNK_VARS 
           if (scrh1(nv) .le. 0.e0) then 
              am(nv, j) = a(nv, j)
              ap(nv, j) = a(nv, j)
           end if
        end do
     end do
     
     do j = 3, n - 3
        scrh2(:) = ap(:,j) - am(:,j)
        scrh3(:) = scrh2*(a(:,j) - 0.5d0*(am(:,j) + ap(:,j)))
        scrh4(:) = scrh2(:)*scrh2(:)/6.e0;
        do nv = 1,NUNK_VARS 
           if ( scrh3(nv) .gt. scrh4(nv)) am(nv, j) = 3.e0*a(nv,j) - 2.e0*ap(nv,j)
           if (-scrh4(nv) .gt. scrh3(nv)) ap(nv, j) = 3.e0*a(nv,j) - 2.e0*am(nv,j)
        end do
     end do
  end if
  ! End of parabolic reconstruction

! --------------------------------------------------------
!      S T E P   # 4: Flattening
!
!  Flattening is applied to either parabolic or linear 
!  interpolations. 
! -------------------------------------------------------- 

  do j = 2, n - 1
    dp(j)    = a(PRES_VAR, j + 1) - a(PRES_VAR, j - 1)
    min_p(j) = min(a(PRES_VAR, j + 1), a(PRES_VAR, j - 1))
  end do  
  do j = 3, n - 2
    d2p(j) = a(PRES_VAR, j + 2) - a(PRES_VAR, j - 2)
  end do  

! ---------------------------------------------------------
!        CHECK TO SEE IF ZONE J IS INSIDE A SHOCK
!  --------------------------------------------------------- 


  do j = 3, n - 2
    beta = abs(dp(j)) / min_p(j)
    dvn = a(ivn, j + 1) - a(ivn, j - 1)
    if (beta < EPS2 .or. dvn > 0.e0) then
      f_t(j) = 0.e0
    else
      fltn   = OME2*(dp(j)/d2p(j) - OME1)
      fltn   = min(1.e0, fltn)
      f_t(j) = max(0.e0, fltn)
    end if
  end do

  do j = 4, n - 3
    if (dp(j) .gt. 0.e0) then 
      sj = -1
    else
      sj =  1
    endif
    fj(j) = max (f_t(j), f_t(j + sj))
  end do

  do j = 4, n - 3
    scrh1(:) = a(:,j)*fj(j)
    scrh2(:) = 1.e0 - fj(j)
    am(:,j) = scrh1(:) + am(:,j)*scrh2(:)
    ap(:,j) = scrh1(:) + ap(:,j)*scrh2(:)
  end do

  ! ------------------------------------------
  !        Relativistic limiter :
  !
  !   revert to 1st order whenever the total
  !   reconstructed velocity exceeds 1
  ! ------------------------------------------

  do j = 4, n - 3

    v2l = am(VELX_VAR,j)*am(VELX_VAR,j) + &
          am(VELY_VAR,j)*am(VELY_VAR,j) + &
          am(VELZ_VAR,j)*am(VELZ_VAR,j)
    v2r = ap(VELX_VAR,j)*ap(VELX_VAR,j) + &
          ap(VELY_VAR,j)*ap(VELY_VAR,j) + &
          ap(VELZ_VAR,j)*ap(VELZ_VAR,j)

    if (v2l .ge. 1.e0 .or. v2r .ge. 1.e0) then
      am(:,j) = a(:,j)
      ap(:,j) = a(:,j)
    end if
  end do

end subroutine hy_rhd_reconstruct
