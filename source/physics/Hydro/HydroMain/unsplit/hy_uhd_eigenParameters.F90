!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_eigenParameters
!!
!! NAME
!!
!!  hy_uhd_eigenParameters
!!
!! SYNOPSIS
!!
!!  hy_uhd_eigenParameters( real (IN)           :: V(HY_VARINUM2),
!!                          integer(IN)         :: dir,
!!                          real (OUT)          :: U_normal,
!!                          real (OUT)          :: C_fast,
!!                          real (OUT),optional :: C_alfn,
!!                          real (OUT),optional :: C_slow,
!!                          real (OUT),optional :: A_f,
!!                          real (OUT),optional :: A_s,
!!                          real (OUT),optional :: B_beta(MDIM),
!!                          real (OUT),optional :: C_hyp)
!!
!! DESCRIPTION
!!
!!  This routine calculates several parameters that are used in calculations of
!!  MHD/Hydro eigenvalues and eigenvectors.
!!
!!
!! ARGUMENTS
!!
!!  V        - Primitive variables + gammas (dens,velx,vely,velz,pres,(magx,magy,magz),gamc,game)
!!  dir      - x,y,z direction
!!  U_normal - Fluid velocity in normal direction
!!  C_fast   - Fast magnetoacoustic speed for MHD/Sound speed for Hydro
!!  C_alfn   - Alfven speed (needed for MHD only)
!!  C_slow   - Slow magnetoacoustic speed (needed for MHD only)
!!  A_f      - Normalization coefficient (needed for MHD only)
!!  A_s      - Normalization coefficient (needed for MHD only)
!!  B_beta   - Alfven velcoities in transversal direction (needed for MHD only)
!!  C_hyp    - advection wave speed for GLM-MHD
!!
!!***


Subroutine hy_uhd_eigenParameters(V,dir,U_normal,C_fast,C_alfn,C_slow,A_f,A_s,B_beta,C_hyp)

  use Hydro_data!,        ONLY : hy_meshMe,hy_forceHydroLimit
  use Driver_data,       ONLY : dr_nStep 
  use Driver_interface,  ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage,Logfile_open,Logfile_close

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! Arguments type declaration ---------------------------
  real, dimension(HY_VARINUM2), intent(IN)  :: V
  integer, intent(IN)  :: dir
  real, intent(OUT) :: U_normal,C_fast
  real, intent(OUT), optional :: C_alfn,C_slow,A_f,A_s
  real, dimension(MDIM), intent(OUT),optional :: B_beta
  real, intent(OUT), optional :: C_hyp
  !! ------------------------------------------------------
  integer :: logUnit
  logical :: logUnitLocal=.true.
  real    :: a2
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  real    :: u2,a,C_slow2,C_fast2,A_s2,A_f2
  real    :: bbx,bby,bbz,bb2,bbT1,bbT2,bbT,sqrtd
#endif

  if (present(C_alfn)) C_alfn = 0.
  if (present(C_slow)) C_slow = 0.
  if (present(A_f)) A_f = 0.
  if (present(A_s)) A_s = 0.
  if (present(B_beta)) B_beta= 0.
  if (present(C_hyp)) C_hyp= 0.


  U_normal = V(HY_VELX+dir-1)
  a2 = V(HY_GAMC)*V(HY_PRES)/V(HY_DENS)
  !C_fast = sqrt(a2)

  if (a2 .le. 0.) then
     call Driver_abortFlash&
          ("[hy_uhd_eigenParameters-A]: Zero or imaginary sound speed has obtained! "//&
           "Please try other (more diffusive) slope limiter, flux, order, cfl, etc.")
  else
     C_fast = sqrt(a2)
  endif

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  !! Alfven velocities
  sqrtd=1./sqrt(V(HY_DENS))
  bbx  = V(HY_MAGX)*sqrtd
  bby  = V(HY_MAGY)*sqrtd
  bbz  = V(HY_MAGZ)*sqrtd
  bb2  = bbx*bbx+bby*bby+bbz*bbz
  !bb2 = dot_product(V(HY_MAGX:HY_MAGZ),V(HY_MAGX:HY_MAGZ))*sqrtd

  !! Normal and tangent components
  if (dir==DIR_X) then
     ! x-sweep
     bbT1     = bby
     bbT2     = bbz
  elseif (dir==DIR_Y) then
     ! y-sweep
     bbT1     = bbz
     bbT2     = bbx
  elseif (dir==DIR_Z) then
     ! z-sweep
     bbT1     = bbx
     bbT2     = bby
  endif

  C_alfn = V(HY_MAGX+dir-1)*sqrtd
  bbT  = sqrt(bbT1*bbT1+bbT2*bbT2)
  u2   = dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))
  C_fast2 = 0.5*(a2 + bb2 + sqrt( (a2-bb2)*(a2-bb2)+4.*a2*bbT**2  ))
  C_slow2 = a2*C_alfn*C_alfn/C_fast2


  if ((a2 .le. 0.) .or. (C_fast2 < C_slow2)) then
     call Driver_abortFlash&
          ("[hy_uhd_eigenParameters-B]: Zero or imaginary sound speed has obtained! "//&
           "Please try other (more diffusive) slope limiter, flux, order, cfl, etc.")
  else
     !! Renormalization coefficients
     if (.not. hy_forceHydroLimit .and. bb2 > 0.) then
        if (abs(C_fast2-C_slow2) > 1.e-16*a2) then
           ! CASE III & CASE IV of Roe & Balara
           ! bbT = 0 & two cases happen depending on the magnitudes of C_alfven
           ! (1) C_alfven > C_sound
           ! (2) C_alfven < C_sound
           A_f2 = min(1., max(0.,(a2-C_slow2)/(C_fast2-C_slow2)))
           A_s2 = 1.-A_f2
        else
           ! CASE V of Roe & Balsara, SIAM, 1996.
           ! C_fast = C_slow case: 
           ! this happens when C_alfven = C_sound & bbT = 0,
           ! then C_fast = C_slow = C_alfven = C_sound.
           ! Also see Stone et al., ApJS, 2008
           A_f2 = 1.
           A_s2 = 0.
        endif
     else
        ! CASE II of Roe & Balsara
        !hydro limit
        A_f2     = 1.
        A_s2     = 0.
        C_fast2  = a2
        C_slow2  = 0.
        C_alfn   = 0.
     endif

     if (bb2 .ne. 0.) then
        !! Avoid indeterminate cases: see Roe & Balsara
        if (bbT > 0.) then
           B_beta = (/bbx, bby, bbz/)/bbT
        else
           B_beta = sqrt(0.5)
        endif
        B_beta(dir)  = 0.
     else
        !! This is the hydro limit with all magnetic components are exactly zero.
        !! In this case, we make sure that the MHD eigensystem reduces to that of
        !! the pure hydro case. See also MHD_StaggeredMesh/hy_uhd_eigenVector.F90.
        B_beta = 0.
     endif


     !! Updates
     C_fast = sqrt(C_fast2)
     C_slow = sqrt(C_slow2)
     A_f    = sqrt(A_f2)
     A_s    = sqrt(A_s2)
  endif

#endif



#ifdef BROKENCODE
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  !! Alfven velocities
  sqrtd=1./sqrt(V(HY_DENS))
  bbx  = V(HY_MAGX)*sqrtd
  bby  = V(HY_MAGY)*sqrtd
  bbz  = V(HY_MAGZ)*sqrtd
  ! bb2  = bbx*bbx+bby*bby+bbz*bbz
  bb2 = dot_product(V(HY_MAGX:HY_MAGZ),V(HY_MAGX:HY_MAGZ))*sqrtd

  !! Normal and tangent components
  if (dir==DIR_X) then
     ! x-sweep
     bbT1     = bby
     bbT2     = bbz
  elseif (dir==DIR_Y) then
     ! y-sweep
     bbT1     = bbz
     bbT2     = bbx
  elseif (dir==DIR_Z) then
     ! z-sweep
     bbT1     = bbx
     bbT2     = bby
  endif

  C_alfn = V(HY_MAGX+dir-1)*sqrtd
  bbT  = sqrt(bbT1*bbT1+bbT2*bbT2)
  u2   = dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))

!!$  if (V(HY_PRES)<0. ) then
!!$     call Logfile_open(logUnit,logUnitLocal)
!!$     write(logUnit,*)'[hy_uhd_eigenParameters] ERROR: negative pressure at nstep=',V(HY_PRES),dr_nStep
!!$     call Logfile_close(logUnitLocal)
!!$     call Driver_abortFlash&
!!$          ("[hy_uhd_eigenParameters] Negative pressure: Please lower CFL or try different limiter.")
!!$  endif
!!$  if (V(HY_DENS)<0. ) then
!!$     call Logfile_open(logUnit,logUnitLocal)
!!$     write(logUnit,*)'[hy_uhd_eigenParameters] ERROR: negative density at nstep=',V(HY_DENS),dr_nStep
!!$     call Logfile_close(logUnitLocal)
!!$     call Driver_abortFlash&
!!$          ("[hy_uhd_eigenParameters] Negative density: Please lower CFL or try different limiter.")
!!$  endif
!!$  if (V(HY_GAMC)<0. ) then
!!$     call Logfile_open(logUnit,logUnitLocal)
!!$     write(logUnit,*)'[hy_uhd_eigenParameters] ERROR: negative gamc at nstep=',V(HY_GAMC),dr_nStep
!!$     call Logfile_close(logUnitLocal)
!!$     call Driver_abortFlash&
!!$          ("[hy_uhd_eigenParameters] Negative gamc: Please lower CFL or try different limiter.")
!!$  endif


!!$  if (a2 .le. 0.) then
!!$     call Driver_abortFlash&
!!$          ("[hy_uhd_eigenParameters]: Zero or imaginary sound speed has obtained! "//&
!!$           "Please try other (more diffusive) slope limiter, flux, order, cfl, etc.")
!!$  else
!     C_fast = sqrt(a2)
!!$  endif

  !! Sound speed and magneto-acoustic speeds
  !a2      = V(HY_GAMC)*V(HY_PRES)/V(HY_DENS)
!!$  if (a2 .le. 0.) then
!!$     call Driver_abortFlash&
!!$          ("[hy_uhd_eigenParameters]: Zero or imaginary sound speed has obtained! "//&
!!$           "Please try other (more diffusive) slope limiter, flux, order, cfl, etc.")
!!$  else

     C_fast2 = 0.5*(a2 + bb2 + sqrt( (a2-bb2)*(a2-bb2)+4.*a2*bbT**2  ))
     C_slow2 = a2*C_alfn*C_alfn/C_fast2


     !! Renormalization coefficients
     if (.not. hy_forceHydroLimit) then
        if (C_fast2-C_slow2 > 0.) then
           A_f2 = min(1., max(0.,(a2-C_slow2)/(C_fast2-C_slow2)))
           A_s2 = 1.-A_f2
        else
           A_f2 = .5
           A_s2 = .5
        endif
     else !hydro limit
        A_f2     = 1.
        A_s2     = 0.
        C_fast2  = a2
        C_slow2  = 0.
        C_alfn   = 0.
     endif

     !! Avoid indeterminate cases: see Roe & Balsara
     if (bbT > 0.)then
        B_beta = (/bbx, bby, bbz/)/bbT
     else
        B_beta = sqrt(0.5)
     endif
     B_beta(dir)  = 0.
 
     !! Updates
     C_fast = sqrt(C_fast2)
     C_slow = sqrt(C_slow2)
     A_f    = sqrt(A_f2)
     A_s    = sqrt(A_s2)
 ! endif
#endif /* for MHD */
#endif /*ifdef DONGWOOK*/

End Subroutine hy_uhd_eigenParameters
