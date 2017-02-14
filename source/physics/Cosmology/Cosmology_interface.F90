!!****h* source/physics/Cosmology/Cosmology_interface
!! 
!! NAME
!!  Cosmology_interface
!! 
!! SYNOPSIS
!!   use Cosmology_interface : ONLY
!!
!! DESCRIPTION
!!   This is the interface definition file for the Cosmology unit.
!!
!!***

Module Cosmology_interface

  implicit none

# include "Flash.h"
# include "constants.h"

    interface
       subroutine Cosmology_cdmPowerSpectrum(k, Anorm, znorm, n, Omega0, h, Lambda0, powerSpectrum)
         implicit none
         real, intent(IN) :: k, Anorm, znorm, n, Omega0, h, Lambda0
         real, intent(OUT) :: powerSpectrum
       end subroutine Cosmology_cdmPowerSpectrum
    end interface

  interface
     subroutine Cosmology_computeVariance (lambda, Mass, Delta0, dDelta0dM, N, f,&
          Anorm, znorm, npspc, Omega0, h, Lambda0)
       implicit none
       integer, intent(in)  :: N
       real,intent(inout)     :: lambda(N), Delta0(N), f
       real,intent(inout)     :: dDelta0dM(N), Mass(N)
       real,intent(inout)     :: Anorm, znorm, npspc, Omega0, h, Lambda0
     end subroutine Cosmology_computeVariance
  end interface

  interface
     subroutine Cosmology_computeDeltaCrit (z, dcrit, dcritdz, D, N, Omega0, H0, Lambda0, Omegatot)
       implicit none
       integer, intent(in) :: N
       real, intent(in)   :: z(N), Omega0, H0, Lambda0, Omegatot 
       real, intent(inout) :: D(N), dcrit(N), dcritdz(N)
     end subroutine Cosmology_computeDeltaCrit
  end interface

  interface
     subroutine Cosmology_computeDt (dt_cosmo)
       implicit none
       real, INTENT(inout)    :: dt_cosmo
     end subroutine Cosmology_computeDt
  end interface

  interface
     subroutine Cosmology_init(restart)
       implicit none
       logical, intent(IN) :: restart
     end subroutine Cosmology_init
  end interface

  interface
     subroutine Cosmology_finalize()
       implicit none
     end subroutine Cosmology_finalize
  end interface

  interface 
     subroutine Cosmology_getParams(hubble,omega,baryon,lambda)
       implicit none
       real, intent(OUT) :: hubble
       real, intent(OUT) :: omega
       real, intent(OUT) :: baryon
       real, intent(OUT) :: lambda
     end subroutine Cosmology_getParams
  end interface

  interface
     subroutine Cosmology_getOldRedshift(zOld)
       implicit none
       real, intent(OUT) :: zOld
     end subroutine Cosmology_getOldRedshift
  end interface

  interface 
     subroutine Cosmology_getRedshift(z)
       implicit none
       real, intent(OUT) :: z
     end subroutine Cosmology_getRedshift
  end interface

  interface
     subroutine Cosmology_massToLength (M, lambda)
       implicit none
       real, intent(IN) :: M
       real, intent(OUT) :: lambda
     end subroutine Cosmology_massToLength
  end interface

  interface
     subroutine Cosmology_redshiftHydro( blockCount, blockList)
       implicit none

       integer, intent(IN) :: blockCount
       integer, dimension(blockCount), intent(IN) :: blockList 

     end subroutine Cosmology_redshiftHydro
  end interface

  interface 
     subroutine Cosmology_redshiftToTime(z,t)
       implicit none
       real, intent(IN) :: z
       real, intent(OUT) :: t
     end subroutine Cosmology_redshiftToTime
  end interface

  interface
     subroutine Cosmology_sendOutputData()
       implicit none
     end subroutine Cosmology_sendOutputData
  end interface


  interface
     subroutine Cosmology_solveFriedmannEqn(timeEndAdvance, dt)
       implicit none
       real, intent(IN) :: timeEndAdvance
       real, intent(IN) :: dt

     end subroutine Cosmology_solveFriedmannEqn
  end interface

  interface
     subroutine Cosmology_unitTest(  fileUnit, perfect )
       implicit none
       integer, intent(IN) :: fileUnit
       logical, intent(INOUT) :: perfect
     end subroutine Cosmology_unitTest
  end interface


end Module Cosmology_interface
