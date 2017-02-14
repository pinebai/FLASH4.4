!!****ih* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_interface
!!
!! NAME
!!  hy_rhd_interface
!!
!! SYNOPSIS
!!  use hy_rhd_interface
!!
!! DESCRIPTION
!!  This is an interface specific for the 8Wave 
!!  MHD module that defines its public interfaces.
!!
!!***
Module hy_rhd_interface

#include "constants.h"
#include "Flash.h"
#include "RHD.h"

  implicit none

  interface
     subroutine hy_rhd_checkBoundaryValues (Vc, n, dir)
       integer, INTENT(IN) :: n, dir
       real, DIMENSION(NUNK_VARS,n), INTENT(IN) :: Vc
     end subroutine hy_rhd_checkBoundaryValues
  end interface


  interface
     subroutine hy_rhd_conserveToPrimitive (U, ibeg, iend, n)
       integer, INTENT(in) :: ibeg, iend, n
       real, DIMENSION(NUNK_VARS,n), INTENT(inout)  :: U
     end subroutine hy_rhd_conserveToPrimitive
  end interface


  interface
     subroutine hy_rhd_enthalpy(u, h, ibeg, iend, n)
       integer, INTENT(in) :: ibeg,iend,n
       real, DIMENSION(NUNK_VARS,n), INTENT(in) :: u
       real, DIMENSION(n), INTENT(out) :: h
     end subroutine hy_rhd_enthalpy
  end interface


  interface
     subroutine hy_rhd_flux (u, flux, ucns, ibeg, iend, n, dir)
       integer, INTENT(in) :: ibeg, iend, n, dir
       real, DIMENSION(NUNK_VARS,n), INTENT(in)  :: u
       real, DIMENSION(NUNK_VARS,n), INTENT(out) :: ucns
       real, DIMENSION(NFLUXES,n), INTENT(out) :: flux
     end subroutine hy_rhd_flux
  end interface


  interface
     subroutine hy_rhd_hlle&
          (Vc, Vm, Vp, Flux, x, dx, dt, speed, vint, ibeg, iend, n, dir)
       implicit none
       integer, INTENT(in) :: n, dir, ibeg, iend
       real, DIMENSION(NUNK_VARS,n), INTENT(in)  :: Vm, Vp, Vc
       real, DIMENSION(NFLUXES,n), INTENT(out) :: Flux
       real, DIMENSION(n), INTENT(in)  :: x, dx
       real, DIMENSION(n), INTENT(out) :: speed, vint
       real, INTENT(in)                :: dt
     end subroutine hy_rhd_hlle
  end interface


  interface
     subroutine hy_rhd_maxChSpeed (u, h, speed, ibeg, iend, n , dir)
       integer, INTENT(in) :: ibeg, iend, n, dir
       real, DIMENSION(NUNK_VARS,n), INTENT(in) :: u
       real, DIMENSION(n), INTENT(in)  :: h
       real, DIMENSION(n), INTENT(out) :: speed
     end subroutine hy_rhd_maxChSpeed
  end interface


  interface
     subroutine hy_rhd_pressureFunc(u, x, fp, dfdp, iflag)
       integer, INTENT(out) :: iflag
       real, DIMENSION(NUNK_VARS),INTENT(in) :: u
       real, INTENT(in)  :: x
       real, INTENT(out) :: fp, dfdp
     end subroutine hy_rhd_pressureFunc
  end interface


  interface
     subroutine hy_rhd_primitiveToConserve (U, ibeg, iend, n)
       integer, INTENT(in) :: ibeg, iend, n
       real, DIMENSION(NUNK_VARS,n), INTENT(inout)  :: U
     end subroutine hy_rhd_primitiveToConserve
  end interface


  interface
     subroutine hy_rhd_reconstruct (a, am, ap, r, dr, dvol, n, dir)
       integer, INTENT(in) :: n, dir
       real, DIMENSION(NUNK_VARS,n), INTENT(in)  :: a
       real, DIMENSION(NUNK_VARS,n), INTENT(out) :: am, ap
       real, DIMENSION(n), INTENT(in)  :: r, dr, dvol
     end subroutine hy_rhd_reconstruct
  end interface


  interface
     subroutine hy_rhd_riemann&
          (Vc, Vm, Vp, Flux, x, dx, dt, speed, vint, n, dir)
       implicit none
       integer, INTENT(in) :: n, dir
       real, DIMENSION(NUNK_VARS,n), INTENT(in) :: Vm, Vp, Vc
       real, DIMENSION(NFLUXES,n), INTENT(out)  :: Flux
       real, DIMENSION(n), INTENT(in)  :: x, dx
       real, DIMENSION(n), INTENT(out) :: speed, vint
       real, INTENT(in)                :: dt
     end subroutine hy_rhd_riemann
  end interface


  interface
     subroutine hy_rhd_setTstep( dtime,i,j,k,blockID)
       real, INTENT(in) :: dtime
       integer, INTENT(in) ::  i,j,k,blockID
     end subroutine hy_rhd_setTstep
  end interface


  interface
     subroutine hy_rhd_shock&
          (tau0, u0, p0, g0, V0, h0, p1, u1,  dudp, zeta, istate)
       integer, INTENT(in) :: istate
       real, INTENT(in)    :: tau0, u0,p0, g0, V0, h0, p1
       real, INTENT(out)   :: u1, dudp, zeta
     end subroutine hy_rhd_shock
  end interface


  interface
     subroutine hy_rhd_soundSpeed2(u, h, a2, ibeg, iend, n)
       integer, INTENT(in) :: ibeg, iend, n
       real, DIMENSION(NUNK_VARS,n), INTENT(in) :: u
       real, DIMENSION(n), INTENT(in)  :: h
       real, DIMENSION(n), INTENT(out) :: a2
     end subroutine hy_rhd_soundSpeed2
  end interface


  interface
     subroutine hy_rhd_sources&
          (Vc, h, cs2, src, r, dr, ibeg, iend, n, dir, eqnType)
       integer, INTENT(IN) :: ibeg, iend, n, dir,eqnType
       real, DIMENSION(NUNK_VARS,n), INTENT(OUT) :: src
       real, DIMENSION(NUNK_VARS,n), INTENT(IN)  :: Vc
       real, DIMENSION(n), INTENT(IN) :: cs2, h, r, dr
     end subroutine hy_rhd_sources
  end interface


  interface
     subroutine hy_rhd_states(Vc, Vm, Vp, grav, dt, r, dr, dvol, n, dir)
       integer, INTENT(IN) :: n, dir
       real, DIMENSION(NUNK_VARS,n), INTENT(IN) :: Vc
       real, DIMENSION(NUNK_VARS,n), INTENT(OUT) :: Vm, Vp
       real, DIMENSION(n), INTENT(IN) :: grav, r, dr, dvol
       real, INTENT(IN)  :: dt
     end subroutine hy_rhd_states
  end interface


  interface
     subroutine hy_rhd_sweep&
          (blockCount,blockList,timeEndAdv, dt, dtOld, sweepDir)
       implicit none
       integer, intent(IN) ::  blockCount
       integer, intent(IN), dimension(blockCount) :: blockList
       real, intent (IN) :: timeEndAdv, dt, dtOld
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_rhd_sweep
  end interface

end Module hy_rhd_interface
