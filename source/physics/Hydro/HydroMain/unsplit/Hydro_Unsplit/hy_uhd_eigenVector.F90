!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_Unsplit/hy_uhd_eigenVector
!!
!! NAME
!!
!!  hy_uhd_eigenVector
!!
!! SYNOPSIS
!!
!!  hy_uhd_eigenVector( real (OUT)         :: LeftEigvec (HY_VARINUM,HY_WAVENUM),
!!                      real (OUT)         :: RightEigvec(HY_VARINUM,HY_WAVENUM)
!!                      real (IN)          :: V(HY_VARINUM2),
!!                      integer(IN)        :: dir,
!!                      logical(IN)        :: cons,
!!                      real (IN)          :: C_fast,
!!                      real (IN),optional :: C_alfn,
!!                      real (IN),optional :: C_slow,
!!                      real (IN),optional :: A_f,
!!                      real (IN),optional :: A_s,
!!                      real (IN),optional :: B_beta(MDIM) )
!!
!! DESCRIPTION
!!
!!  This routine calculates MHD/Hydro eigenvectors in either primitive form (used in
!!  Riemann solver) or conservative form (used in conservative updates).
!!
!!
!! ARGUMENTS
!!
!!  LeftEigvec  - Left eigenvectors
!!  RightEigvec - Right eigenvectors
!!  V           - Primitive variables + gammas:
!!                (dens,velx,vely,velz,pres,(magx,magy,magz),gamc,game)
!!  dir         - x,y,z direction
!!  cons        - A logical switch to choose either primitive or conservative eigenvector
!!  C_fast      - Fast magnetoacoustic speed for MHD/Sound speed for Hydro
!!  C_alfn      - Alfven speed (needed for MHD only)
!!  C_slow      - Slow magnetoacoustic speed (needed for MHD only)
!!  A_f         - Normalization coefficient (needed for MHD only)
!!  A_s         - Normalization coefficient (needed for MHD only)
!!  B_beta      - Alfven velcoities in transversal direction (needed for MHD only)
!!
!!***

!#define DEBUG_HY_EIGEN

Subroutine hy_uhd_eigenVector&
     (LeftEigvec,RightEigvec,V,dir,cons,C_fast,C_alfn,C_slow,A_f,A_s,B_beta)

  use Hydro_data,           ONLY : hy_meshMe
  use Logfile_interface,    ONLY : Logfile_open,Logfile_close

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! Arguments type declaration ---------------------------------------
  real, dimension(HY_VARINUM,HY_WAVENUM), intent(OUT) :: LeftEigvec
  real, dimension(HY_VARINUM,HY_WAVENUM), intent(OUT) :: RightEigvec
  real, dimension(HY_VARINUM2), intent(IN) :: V
  integer, intent(IN) :: dir
  logical, intent(IN) :: cons
  real,    intent(IN) :: C_fast
  real,    intent(IN), optional :: C_alfn,C_slow,A_f,A_s
  real, dimension(MDIM), intent(IN), optional  :: B_beta
  !! ------------------------------------------------------------------

  integer :: ii,jj,logUnit
  integer :: HY_VEL1,HY_VEL2,HY_VEL3
  real    :: k,dinv,a2inv,a,a2,cf2,u2,ahinv,test
  real    :: H,Na
  logical :: logUnitLocal=.true.



  if (dir == DIR_X) then
     HY_VEL1=HY_VELX
     HY_VEL2=HY_VELY
     HY_VEL3=HY_VELZ
  elseif (dir == DIR_Y) then
     HY_VEL1=HY_VELY
     HY_VEL2=HY_VELZ
     HY_VEL3=HY_VELX
  elseif (dir == DIR_Z) then
     HY_VEL1=HY_VELZ
     HY_VEL2=HY_VELX
     HY_VEL3=HY_VELY
  endif


  ! initialize with zeros:
  !! eigen vectors should be initialized with zeros because 
  !! we only compute non-zero entries here.
  LeftEigVec  = 0.
  RightEigVec = 0.


  ! parameters
  dinv  =1./V(HY_DENS)
  k=1.-V(HY_GAME) !! k=1.-game

  a2    = C_fast*C_fast
  a2inv = 1./a2
  ahinv = 0.5/C_fast
  Na    = 0.5*a2inv

  if (.not. cons) then
     !! -----------------------------------------------!
     !! Construct left eigenvectors                    !
     !! -----------------------------------------------!

     ! left going fast wave
     LeftEigvec(HY_VEL1,HY_FASTLEFT) =-ahinv
     LeftEigvec(HY_PRES,HY_FASTLEFT) = dinv*Na

     ! left going slow wave
     LeftEigvec(HY_VEL2,HY_SLOWLEFT) = 1. ! -V(HY_DENS)

     ! entropy wave
     LeftEigvec(HY_DENS,HY_ENTROPY ) = 1.
     LeftEigvec(HY_PRES,HY_ENTROPY ) =-a2inv

     ! right going slow wave
     LeftEigvec(HY_VEL3,HY_SLOWRGHT) = 1. ! V(HY_DENS)

     ! right going fast wave
     LeftEigvec(HY_VEL1,HY_FASTRGHT) = ahinv
     LeftEigvec(HY_PRES,HY_FASTRGHT) = LeftEigvec(HY_PRES,HY_FASTLEFT)

     !! -----------------------------------------------!
     !! Construct right eigenvectors                   !
     !! -----------------------------------------------!
     ! left going fast wave
     RightEigvec(HY_DENS,HY_FASTLEFT) = V(HY_DENS)
     RightEigvec(HY_VEL1,HY_FASTLEFT) =-C_fast
     RightEigvec(HY_PRES,HY_FASTLEFT) = V(HY_DENS)*a2

     ! left going slow wave
     RightEigvec(HY_VEL2,HY_SLOWLEFT) = 1. !-dinv

     ! entropy wave
     RightEigvec(HY_DENS,HY_ENTROPY ) = 1.

     ! right going slow wave
     RightEigvec(HY_VEL3,HY_SLOWRGHT) = 1. !dinv

     ! right going fast wave
     RightEigvec(HY_DENS,HY_FASTRGHT) = V(HY_DENS)
     RightEigvec(HY_VEL1,HY_FASTRGHT) = C_fast
     RightEigvec(HY_PRES,HY_FASTRGHT) = RightEigvec(HY_PRES,HY_FASTLEFT)


#ifdef DEBUG_HY_EIGEN
     do ii = 1,HY_WAVENUM
        test=dot_product(LeftEigvec(:,ii), RightEigvec(:,ii))
        if (abs(test-1.) > 1.e-4) then
           call Logfile_open(logUnit,logUnitLocal)
           write(logUnit,*)'dot product is not unity: test=',test,ii,dir   !DEBUG
           call Logfile_close(logUnitLocal)
        endif
     enddo

     do ii = 1,HY_WAVENUM
        if (ii < HY_WAVENUM) then
           test=dot_product(LeftEigvec(:,ii), RightEigvec(:,ii+1))
           if (abs(test) > 1.e-4) then
              call Logfile_open(logUnit,logUnitLocal)
              write(logUnit,*)'dot product is not zero: test=',test,ii,dir   !DEBUG
              call Logfile_close(logUnitLocal)
           endif
        else
           test=dot_product(LeftEigvec(:,ii), RightEigvec(:,ii-1))
           if (abs(test) > 1.e-4) then
              call Logfile_open(logUnit,logUnitLocal)
              write(logUnit,*)'dot product is not zero: test=',test,ii,dir   !DEBUG
              call Logfile_close(logUnitLocal)
           endif
        endif
     enddo
#endif


  else ! if conserved

     ! enthalpy H = (E+p)/rho
     u2=0.5*dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))
     H = (V(HY_DENS)*u2-V(HY_PRES)/k + V(HY_PRES))/V(HY_DENS)

     !! -----------------------------------------------!
     !! Construct left eigenvectors                    !
     !! -----------------------------------------------!
     ! left going fast wave
     LeftEigvec(HY_DENS,HY_FASTLEFT) = Na*(-k*u2+V(HY_VEL1)*C_fast)
     LeftEigvec(HY_VEL1,HY_FASTLEFT) =-Na*(-k*V(HY_VEL1)+C_fast)
     LeftEigvec(HY_VEL2,HY_FASTLEFT) = Na*k*V(HY_VEL2)
     LeftEigvec(HY_VEL3,HY_FASTLEFT) = Na*k*V(HY_VEL3)
     LeftEigvec(HY_PRES,HY_FASTLEFT) =-Na*k

     ! left going slow wave
     LeftEigvec(HY_DENS,HY_SLOWLEFT) =-V(HY_VEL2)
     !LeftEigvec(HY_VEL1,HY_SLOWLEFT) = 0.
     LeftEigvec(HY_VEL2,HY_SLOWLEFT) = 1.
     !LeftEigvec(HY_VEL3,HY_SLOWLEFT) = 0.
     !LeftEigvec(HY_PRES,HY_SLOWLEFT) = 0.

     ! entropy wave
     LeftEigvec(HY_DENS,HY_ENTROPY ) = 1.+Na*k*u2*2.
     LeftEigvec(HY_VEL1,HY_ENTROPY ) =-k*V(HY_VEL1)*a2inv
     LeftEigvec(HY_VEL2,HY_ENTROPY ) =-k*V(HY_VEL2)*a2inv
     LeftEigvec(HY_VEL3,HY_ENTROPY ) =-k*V(HY_VEL3)*a2inv
     LeftEigvec(HY_PRES,HY_ENTROPY ) = k*a2inv

     ! right going slow wave
     LeftEigvec(HY_DENS,HY_SLOWRGHT) =-V(HY_VEL3) 
     !LeftEigvec(HY_VEL1,HY_SLOWRGHT) = 0.
     !LeftEigvec(HY_VEL2,HY_SLOWRGHT) = 0.
     LeftEigvec(HY_VEL3,HY_SLOWRGHT) = 1.
     !LeftEigvec(HY_PRES,HY_SLOWRGHT) = 0.

     ! right going fast wave
     LeftEigvec(HY_DENS,HY_FASTRGHT) = Na*(-k*u2-V(HY_VEL1)*C_fast)
     LeftEigvec(HY_VEL1,HY_FASTRGHT) =-Na*(-k*V(HY_VEL1)-C_fast)
     LeftEigvec(HY_VEL2,HY_FASTRGHT) = Na*k*V(HY_VEL2)
     LeftEigvec(HY_VEL3,HY_FASTRGHT) = Na*k*V(HY_VEL3)
     LeftEigvec(HY_PRES,HY_FASTRGHT) =-Na*k


     !! -----------------------------------------------!
     !! Construct right eigenvectors                   !
     !! -----------------------------------------------!
     ! left going fast wave
     RightEigvec(HY_DENS,HY_FASTLEFT) = 1.
     RightEigvec(HY_VEL1,HY_FASTLEFT) = V(HY_VEL1) - C_fast
     RightEigvec(HY_VEL2,HY_FASTLEFT) = V(HY_VEL2)
     RightEigvec(HY_VEL3,HY_FASTLEFT) = V(HY_VEL3)
     RightEigvec(HY_PRES,HY_FASTLEFT) = H - V(HY_VEL1)*C_fast

     ! left going slow wave
     !RightEigvec(HY_DENS,HY_SLOWLEFT) = 0.
     !RightEigvec(HY_VEL1,HY_SLOWLEFT) = 0.
     RightEigvec(HY_VEL2,HY_SLOWLEFT) = 1.
     !RightEigvec(HY_VEL3,HY_SLOWLEFT) = 0.
     RightEigvec(HY_PRES,HY_SLOWLEFT) = V(HY_VEL2)

     ! entropy wave
     RightEigvec(HY_DENS,HY_ENTROPY ) = 1.
     RightEigvec(HY_VEL1,HY_ENTROPY ) = V(HY_VEL1)
     RightEigvec(HY_VEL2,HY_ENTROPY ) = V(HY_VEL2)
     RightEigvec(HY_VEL3,HY_ENTROPY ) = V(HY_VEL3)
     RightEigvec(HY_PRES,HY_ENTROPY ) = u2

     ! right going slow wave
     !RightEigvec(HY_DENS,HY_SLOWRGHT) = 0.
     !RightEigvec(HY_VEL1,HY_SLOWRGHT) = 0.
     !RightEigvec(HY_VEL2,HY_SLOWRGHT) = 0.
     RightEigvec(HY_VEL3,HY_SLOWRGHT) = 1.
     RightEigvec(HY_PRES,HY_SLOWRGHT) = V(HY_VEL3)

     ! right going fast wave
     RightEigvec(HY_DENS,HY_FASTRGHT) = 1.
     RightEigvec(HY_VEL1,HY_FASTRGHT) = V(HY_VEL1) + C_fast
     RightEigvec(HY_VEL2,HY_FASTRGHT) = V(HY_VEL2)
     RightEigvec(HY_VEL3,HY_FASTRGHT) = V(HY_VEL3)
     RightEigvec(HY_PRES,HY_FASTRGHT) = H + V(HY_VEL1)*C_fast


  endif

End Subroutine hy_uhd_eigenVector
