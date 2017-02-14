!!****ih* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/multiTemp/hy_ppm_kernelInterface
!!
!! NAME
!!
!!  hy_ppm_kernelInterface
!!
!! SYNOPSIS
!!
!!  use hy_ppm_kernelInterface
!!
!! DESCRIPTION
!!
!!  Interface module for internal use within the split PPM Hydro implementation -
!!  PPM kernel routines.
!!
!!***

Module hy_ppm_kernelInterface
  implicit none
#include "Hydro_components.h"
#include "constants.h"

  interface
     subroutine hydro_1d (blockID,numIntCells,numCells, guard,bcs,        &
          xyzswp, hy_meshMe, dt, dt_old,                 &
          jCell, kCell,                             &
          igeom, useGravity,                             &
          xbot, xtop,                               &
          ybot, ytop, ylft, yrgt,                   &
          zlft, zrgt, ugrid,                        &
          primaryCoord ,     &
          primaryLeftCoord , &
          primaryRghtCoord , &
          primaryDx        , &
          secondCoord      , &
          thirdCoord       , &
          radialCoord     , &
          u, ut, utt, rho, p, e, eintIn, tmp, game, gamc,   &
          xn, utbt, uttp, utlt, utrt,               &
          shock_multid,                             &
          dtdx, areaLeft, area, cvol, grav, ngrav, fict, &
          rhoflx, uflx, pav, utflx, uttflx,         &
          eflx, eintflx, eiaflx, oneflx, voFlx, xnflx)
       use Hydro_data, ONLY : hy_dirGeom, hy_movingGrid, hy_useCellAreasForFluxes, &
            hy_epsiln, hy_omg1, hy_omg2
       use Hydro_data, ONLY : hy_cvisc, hy_useCmaAdvection, hy_smallx, hy_smallp
       use Hydro_data, ONLY : hy_numXn, hy_hybridRiemann
       use Hydro_data, ONLY : hy_ppmEnerFluxConstructionMeth
       use Hydro_data, ONLY : hy_eMass, hy_pMass
       implicit none
       !--arguments-------------------------
       integer, intent(IN) ::  blockID,jCell, kCell, numIntCells,numCells,&
            xyzswp, hy_meshMe,igeom, guard
       real, intent(IN) :: dt, dt_old

       logical, intent(IN) :: useGravity
       integer,intent(IN),dimension(2,MDIM):: bcs

       real, DIMENSION(numCells), intent(IN) :: rho, u, ut, utt, tmp
       real, DIMENSION(numCells), intent(IN) ::  primaryCoord ,     &
            primaryLeftCoord , &
            primaryRghtCoord , &
            primaryDx        , &
            secondCoord      , &
            thirdCoord       , &
            radialCoord     , &
            shock_multid
       real, DIMENSION(numCells,0:HYDRO_NUM_E_COMPONENTS),intent(IN) :: p, e
       real, DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS),intent(IN) :: eintIn
       real, DIMENSION(numCells), intent(IN)    :: cvol
       real, DIMENSION(numCells), intent(INOUT) :: grav, areaLeft
       real, DIMENSION(numCells), intent(OUT)   :: area
       real, DIMENSION(numCells, hy_numXn), intent(OUT) :: xnflx
       real, DIMENSION(numCells), intent(OUT) :: dtdx, ngrav, fict
       real, DIMENSION(numCells), intent(OUT) :: rhoflx, uflx, &
            utflx, uttflx
       real, DIMENSION(numCells), intent(OUT) :: oneflx, voFlx
       real, DIMENSION(numCells,0:HYDRO_NUM_E_COMPONENTS), intent(OUT) :: eflx, pav
       real, DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS), intent(OUT) :: eintflx
       real, DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS), intent(OUT) :: eiaflx



       real, DIMENSION(numCells, hy_numXn), intent(INOUT) :: xn
       real, intent(IN) :: xbot, xtop, ybot 
       real, intent(IN) :: ytop, ylft, yrgt, zlft, zrgt
       real, intent(IN), DIMENSION(numCells) :: ugrid
       real, DIMENSION(numCells,0:HYDRO_NUM_GAME_COMPONENTS), intent(INout) :: game
       real, intent(INOUT), DIMENSION(numCells) :: gamc, utbt, &
            uttp, utlt, utrt
     end subroutine hydro_1d
  end interface


  interface
     subroutine intrfc(sweepDir,numIntCells, numCells, guard, &
          &            rho, u, ut, utt, p, &
          &            rhol, rhor, ul, ur, &
          &            utl, utr, uttl, uttr, &
          &            pl, pr, vl, vr, gamcl, &
          &            gamcr, game, &
          &            gamel,gamer,gamc,grav,&
                       eint, eintl, eintr,etot,etotl,etotr,pFrac,pFracl,pFracr, xn, &
          &            xnl, xnr, v, dx, x, tmp)

       use Hydro_data,   ONLY:  hy_numXn, hy_smlrho, hy_smallx, hy_small, hy_gravl, &
            hy_drho, hy_rho6, &
            hy_du, hy_u6, hy_dut, hy_ut6, &
            hy_dutt, hy_utt6, hy_dp, hy_dgame, hy_game6, &
            hy_gravr, hy_p6, hy_dgamc, hy_gamc6, hy_dgrav, &
            hy_grav6, hy_dxn, hy_xn6, &
            hy_deint, hy_eint6, hy_detot,hy_etot6,hy_dpFrac,hy_pFrac6, &
            hy_pwcubic, hy_pwl, hy_pwr, hy_dpw, hy_pw6l, hy_pw6r, &
            hy_ppmModifystates, &
            hy_useSteepening, hy_useCmaFlattening,&
            hy_epsiln, hy_omg1, hy_omg2
       implicit none
       integer, intent(IN) :: sweepDir,numIntCells, numCells, guard
       real, intent(IN),    DIMENSION(numCells,hy_numXn) :: xn
       real, intent(INOUT), DIMENSION(numCells,hy_numXn) :: xnl, xnr
       real, intent(IN), DIMENSION(numCells) :: &
            rho, u, ut, utt, p, gamc, grav, dx, x, tmp
       real, intent(INOUT), DIMENSION(numCells) :: &
            rhol, rhor, &
            ul, ur, &
            utl, utr, &
            uttl, uttr, &
            pl, pr, &
            gamcl, gamcr
       real, intent(OUT), DIMENSION(numCells) :: &
            v, vl, vr
       real,intent(IN),DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS) :: eint
       real,intent(IN),DIMENSION(numCells,0:HYDRO_NUM_E_COMPONENTS) :: etot
       real,intent(IN),DIMENSION(numCells,1:HYDRO_NUM_E_COMPONENTS) :: pFrac
       real,intent(IN),DIMENSION(numCells,0:HYDRO_NUM_GAME_COMPONENTS) :: game
       real,intent(INOUT),DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS) :: eintl, eintr
       real,intent(INOUT),DIMENSION(numCells,0:HYDRO_NUM_E_COMPONENTS) :: etotl, etotr
       real,intent(INOUT),DIMENSION(numCells,1:HYDRO_NUM_E_COMPONENTS) :: pFracl, pFracr
       real,intent(INOUT),DIMENSION(numCells,0:HYDRO_NUM_GAME_COMPONENTS) :: gamel, gamer
     end subroutine intrfc
  end interface

  interface
     subroutine states (numIntCells, numCells, &
          j, igeom,&
          rho, u, rhol, rhor, ul, ur, &
          utl, utr, uttl, uttr, p, pl, pr, &
          gamcl, gamcr, &
          ugrid, ce, game, gamer, gamc, gamel, &
          eintl, eintr, &
          etotl,etotr, &
          pFracl,pFracr, &
          xnl, xnr, &
          dtdx, dt, &
          x, xl, radial_coord, grav, fict)
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
       integer, intent(IN):: j, igeom, numIntCells, numCells
       real, intent(IN) :: dt
       real, intent(IN), dimension(numCells) :: rho, u, p, ugrid, x, xl, &
            radial_coord
       real, intent(INOUT), dimension(numCells) :: rhol, rhor, ul, ur
       real, intent(INOUT), dimension(numCells) :: utl, utr, uttl, uttr, pl, pr
       real, intent(IN),    dimension(numCells) :: gamcl, gamc, gamcr
       real, intent(IN),    dimension(numCells) :: ce, game
       real, intent(IN),    DIMENSION(numCells,0:HYDRO_NUM_GAME_COMPONENTS) :: gamel, gamer
       real, intent(IN),    DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS) :: eintl, eintr
       real, intent(IN),    DIMENSION(numCells,0:HYDRO_NUM_E_COMPONENTS) :: etotl, etotr
       real, intent(IN),    DIMENSION(numCells,1:HYDRO_NUM_E_COMPONENTS) :: pFracl, pFracr
       real, intent(IN),    dimension(numCells) :: dtdx
       real, intent(IN),    dimension(numCells) :: grav, fict
       real, intent(INOUT), dimension(numCells, hy_numXn) :: xnl, xnr
     end subroutine states
  end interface

  interface
     subroutine rieman (numIntCells, numCells, &
          rhoav, uav, utav, uttav, pav, &
          urell, ugrdl, kappa, game, gameav, pFracAv,eintAv,etotAv, xnav, x)

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
       implicit none
       integer, intent (IN) :: numIntCells,numCells
       real, intent(IN), DIMENSION(numCells) :: x
       real, intent(IN), DIMENSION(numCells) :: ugrdl, game
       real, intent(OUT), DIMENSION(numCells) :: uav, rhoav, utav, uttav, pav, &
            urell
       real, intent(OUT), DIMENSION(numCells) :: kappa
       real, intent(OUT), DIMENSION(numCells,hy_numXn) :: xnav
       real, intent(OUT), DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS) :: eintAv
       real, intent(OUT), DIMENSION(numCells,0:HYDRO_NUM_E_COMPONENTS) :: etotAv
       real, intent(OUT), DIMENSION(numCells,1:HYDRO_NUM_E_COMPONENTS) :: pFracAv
       real, intent(OUT), DIMENSION(numCells,0:HYDRO_NUM_EINT_COMPONENTS) :: gameav
     end subroutine rieman
  end interface

  interface
     subroutine hy_dbgWriteStates(vname,var,varl,varr, &
          varav, varflx,&
          rho,rhol,rhor, rhoav,&
          urell,dt,&
          numIntCells, numCells, nguard,&
          xCoords,xL,xR,blockID)
       implicit none
       character(len=*),intent(IN) :: vname
       integer,intent(IN) :: numIntCells, numCells, nguard
       real,intent(IN),dimension(numCells) :: var,varl,varr
       real,intent(IN),dimension(numCells) :: varav,varflx
       real,intent(IN),dimension(numCells) :: rho,rhol,rhor,rhoav
       real,intent(IN),dimension(numCells) :: xCoords,xL,xR,urell
  real,intent(IN) :: dt
       integer,intent(IN) :: blockID
     end subroutine hy_dbgWriteStates
  end interface

end Module hy_ppm_kernelInterface
