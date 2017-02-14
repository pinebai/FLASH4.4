!!****if* source/Simulation/SimulationMain/FlameChannel/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!  call Simulation_initBlock(integer(IN) :: blockID)
!!
!! DESCRIPTION
!!  This routine applies initial conditions to the specified block.
!!
!!  Initialization of flame in pressure equilibrium.
!!  The flame front can be planar or spherical in multiple dimensions
!!  depending on the value of the pseudo_1d parameter.
!!  See description of parameters in Config file for more info.
!! 
!! ARGUMENTS
!!
!!  blockID -         the number of the block to update
!!
!! PARAMETERS
!!
!!  eosModeInit -     after this routine sets up initial conditions,
!!                    the grid package calls Eos to insure the
!!                    values are thermodynamically consistent.  This
!!                    parameter controls the mode of application of
!!                    the Eos.
!!  others      -   See description of in Config file for more info.
!!
!! SEE ALSO
!!
!!  Eos_wrapped
!!
!! HISTORY
!!  Dean Townsley 2008
!!***
!

subroutine Simulation_initBlock(blockID)
  
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkBoundBox, Grid_getDeltas, Grid_putPointData, &
    Grid_getCellCoords, Grid_getDeltas
  use Flame_interface, ONLY : Flame_getProfile, Flame_rhJump
  use fl_effData, ONLY : fl_effDeltae
  use Eos_interface, ONLY : Eos

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  
  integer, intent(in) :: blockID

  integer :: i, j, k

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: cell
  integer :: isizeGC, jsizeGC, ksizeGC
  real, allocatable, dimension(:) :: iCenter, jCenter, kCenter
  real, allocatable, dimension(:) :: iLeft, jLeft, kLeft
  real, allocatable, dimension(:) :: iRight, jRight, kRight
  real, dimension(MDIM) :: deltas

  real, dimension(EOS_NUM) :: state
  real :: dyi_qn, dqbar_qn, velx, vely, velz
  real :: flam, flam_l, flam_r, flam_c, p_l, p_r, p_c
  real :: rho_l, rho_r, rho_c, e_l, e_r, e_c, vx_l, vx_r, vx_c
  real :: vy_l, vy_c, vy_r, vz_l, vz_c, vz_r
  real :: fsurf_x_position
  real :: fsurf_distance_c, fsurf_distance_l, fsurf_distance_r
  real :: ye, yi
  real :: kine, alpha, alpha_inv
  real :: zctrVortex, dp_vortex, vortex_stream1, vortex_stream2
  real :: y_vort1, y_vort2, z_vort

!==============================================================================

  ! get essential info about this block - index limits and cell coordinates
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  isizeGC = blkLimitsGC(HIGH,IAXIS)
  allocate(iCenter(isizeGC))
  allocate(iLeft(isizeGC))
  allocate(iRight(isizeGC))
  jsizeGC = blkLimitsGC(HIGH,JAXIS)
  allocate(jCenter(jsizeGC))
  allocate(jLeft(jsizeGC))
  allocate(jRight(jsizeGC))
  ksizeGC = blkLimitsGC(HIGH,KAXIS)
  allocate(kCenter(ksizeGC))
  allocate(kLeft(ksizeGC))
  allocate(kRight(ksizeGC))
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,iCenter,isizeGC)
  call Grid_getCellCoords(IAXIS,blockID,LEFT_EDGE,.true.,iLeft,isizeGC)
  call Grid_getCellCoords(IAXIS,blockID,RIGHT_EDGE,.true.,iRight,isizeGC)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,jCenter,jsizeGC)
  call Grid_getCellCoords(JAXIS,blockID,LEFT_EDGE,.true.,jLeft,jsizeGC)
  call Grid_getCellCoords(JAXIS,blockID,RIGHT_EDGE,.true.,jRight,jsizeGC)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,kCenter,ksizeGC)
  call Grid_getCellCoords(KAXIS,blockID,LEFT_EDGE,.true.,kLeft,ksizeGC)
  call Grid_getCellCoords(KAXIS,blockID,RIGHT_EDGE,.true.,kRight,ksizeGC)

  call Grid_getDeltas(blockID, deltas)

  !-----------------------------------------------
  ! loop over all zones and init
  !-----------------------------------------------
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           if (.not. sim_ignite) then
              ! no burned material, only unburned
              state(:) = sim_eosData_u(:)
              flam = 0.0
              velx = sim_inflowVx
                 
           else ! ignite
              !-----------------------------------------------
              ! initialize, including a burned region
              !-----------------------------------------------

              ! planar flame surface with normal tilted up by angle theta from x direction
              fsurf_x_position = sim_fracPerturb*(sim_xmax-sim_xmin)
              fsurf_distance_c = iCenter(i) - fsurf_x_position
              fsurf_distance_l = iLeft(i) - fsurf_x_position
              fsurf_distance_r = iRight(i) - fsurf_x_position

              ! determine local state in this zone
              ! assume deltas are equal (cells are cuboid) and don't worry
              ! too much about the corners
              if ( fsurf_distance_l > 1.5*sim_laminarWidth ) then
                 ! whole cell unburned material
                 state(:) = sim_eosData_u(:)
                 flam = 0.0
              else if ( fsurf_distance_r < -1.5*sim_laminarWidth ) then
                 ! fully burned to NSE
                 state(:) = sim_eosData_b(:)
                 flam = 1.0
              else
                 ! partially burned
                 ! at least one cell will fall here (necessary to get initial refinement right)
!                 ! get left density
!                 call Flame_getProfile(fsurf_distance_l, flam_l)
!
!                 ! calculate propertise for partially burned material
!                 ! note, in fact ye_f and ye_a should be equal
!                 yi = 1.0/sim_eosData_u(EOS_ABAR)*(1.0-flam_l) + &
!                    1.0/sim_eosData_b(EOS_ABAR)*flam_l
!                 ye = sim_eosData_u(EOS_ZBAR)/sim_eosData_u(EOS_ABAR)*(1.0-flam_l) + &
!                    sim_eosData_b(EOS_ZBAR)/sim_eosData_b(EOS_ABAR)*flam_l
!                 state(:) = sim_eosData_u(:)
!                 state(EOS_ABAR) = 1.0/yi
!                 state(EOS_ZBAR) = ye/yi
!
!                 ! put this in pressure equilibrium with unburned material
!                 call Flame_rhJump(sim_eosData_u, state, flam_l*fl_effDeltae, sim_flamespeed, MODE_DENS_TEMP)
!
!                 rho_l = state(EOS_DENS)
!                 p_l = state(EOS_PRES)
!                 vx_l = sim_inflowVx - sim_flamespeed * &
!                    ( sim_eosData_u(EOS_DENS) / rho_l - 1.0 )
!                 e_l = state(EOS_EINT) + 0.5*vx_l**2 - flam_l*fl_effDeltae
!
!                 ! get right density
!                 call Flame_getProfile(fsurf_distance_r, flam_r)
!
!                 ! calculate propertise for partially burned material
!                 ! note, in fact ye_f and ye_a should be equal
!                 yi = 1.0/sim_eosData_u(EOS_ABAR)*(1.0-flam_l) + &
!                    1.0/sim_eosData_b(EOS_ABAR)*flam_l
!                 ye = sim_eosData_u(EOS_ZBAR)/sim_eosData_u(EOS_ABAR)*(1.0-flam_l) + &
!                    sim_eosData_b(EOS_ZBAR)/sim_eosData_b(EOS_ABAR)*flam_l
!                 state(:) = sim_eosData_u(:)
!                 state(EOS_ABAR) = 1.0/yi
!                 state(EOS_ZBAR) = ye/yi
!
!                 ! put this in pressure equilibrium with unburned material
!                 call Flame_rhJump(sim_eosData_u, state, flam_r*fl_effDeltae, sim_flamespeed, MODE_DENS_TEMP)
!
!                 rho_r = state(EOS_DENS)
!                 p_r = state(EOS_PRES)
!                 vx_r = sim_inflowVx - sim_flamespeed * &
!                    ( sim_eosData_u(EOS_DENS) / rho_r - 1.0 )
!                 e_r = state(EOS_EINT) + 0.5*vx_r**2 - flam_r*fl_effDeltae
!
!                 ! get cell-centered density
!                 call Flame_getProfile(fsurf_distance_c, flam_c)
!
!                 ! calculate propertise for partially burned material
!                 ! note, in fact ye_f and ye_a should be equal
!                 yi = 1.0/sim_eosData_u(EOS_ABAR)*(1.0-flam_l) + &
!                    1.0/sim_eosData_b(EOS_ABAR)*flam_l
!                 ye = sim_eosData_u(EOS_ZBAR)/sim_eosData_u(EOS_ABAR)*(1.0-flam_l) + &
!                    sim_eosData_b(EOS_ZBAR)/sim_eosData_b(EOS_ABAR)*flam_l
!                 state(:) = sim_eosData_u(:)
!                 state(EOS_ABAR) = 1.0/yi
!                 state(EOS_ZBAR) = ye/yi
!
!                 ! put this in pressure equilibrium with unburned material
!                 call Flame_rhJump(sim_eosData_u, state, flam_c*fl_effDeltae, sim_flamespeed, MODE_DENS_TEMP)
!
!                 rho_c = state(EOS_DENS)
!                 p_c = state(EOS_PRES)
!                 vx_c = sim_inflowVx - sim_flamespeed * &
!                    ( sim_eosData_u(EOS_DENS) / rho_c - 1.0 )
!                 e_c = state(EOS_EINT) + 0.5*vx_c**2 - flam_c*fl_effDeltae
!
!                 ! init velocity field, calculate cell average
!                 state(EOS_DENS) = ( rho_l + 4.0*rho_c + rho_r ) / 6.0
!                 state(EOS_PRES) = ( p_l + 4.0*p_c + p_r ) / 6.0
!                 velx = ( vx_l*rho_l + 4.0*vx_c*rho_c + vx_r*rho_r ) / 6.0 / &
!                    state(EOS_DENS)
!                 flam = ( flam_l*rho_l + 4.0*flam_c*rho_c + flam_r*rho_r ) / 6.0 / &
!                    state(EOS_DENS)
!                 state(EOS_EINT) = ( e_l*rho_l + 4.0*e_c*rho_c + e_r*rho_r ) / &
!                    6.0 / state(EOS_DENS) - 0.5*velx**2 + flam*fl_effDeltae
!
!                 ! calculate propertise for partially burned material
!                 ! note, in fact ye_f and ye_a should be equal
!                 yi = 1.0/sim_eosData_u(EOS_ABAR)*(1.0-flam) + &
!                    1.0/sim_eosData_b(EOS_ABAR)*flam
!                 ye = sim_eosData_u(EOS_ZBAR)/sim_eosData_u(EOS_ABAR)*(1.0-flam) + &
!                    sim_eosData_b(EOS_ZBAR)/sim_eosData_b(EOS_ABAR)*flam
!                 state(EOS_ABAR) = 1.0/yi
!                 state(EOS_ZBAR) = ye/yi
!
!                 !call Eos(MODE_DENS_EI,1, state)
!                 call Eos(MODE_DENS_PRES,1, state)
!
!                 velx = sim_inflowVx - sim_flamespeed * &
!                    ( sim_eosData_u(EOS_DENS) / state(EOS_DENS) - 1.0 )
!
                 ! ****** just set it to the center *******!
                 call Flame_getProfile(fsurf_distance_c, flam)

                 ! calculate propertise for partially burned material
                 ! note, in fact ye_f and ye_a should be equal
                 yi = 1.0/sim_eosData_u(EOS_ABAR)*(1.0-flam) + &
                    1.0/sim_eosData_b(EOS_ABAR)*flam
                 ye = sim_eosData_u(EOS_ZBAR)/sim_eosData_u(EOS_ABAR)*(1.0-flam) + &
                    sim_eosData_b(EOS_ZBAR)/sim_eosData_b(EOS_ABAR)*flam
                 state(:) = sim_eosData_u(:)
                 state(EOS_ABAR) = 1.0/yi
                 state(EOS_ZBAR) = ye/yi

                 ! put this in pressure equilibrium with unburned material
                 call Flame_rhJump(sim_eosData_u, state, flam*fl_effDeltae, sim_flamespeed, MODE_DENS_TEMP)

              endif

              alpha = sim_eosData_u(EOS_DENS) / state(EOS_DENS)
              alpha_inv = 1.0 / alpha
              velx = sim_inflowVx - sim_flamespeed * ( alpha - 1.0 )

           endif ! sim_ignite

           if ( (iCenter(i) < sim_xBeginVortex) .OR. &
                (iCenter(i) > sim_xEndVortex ) ) then

              vely = 0.0
              velz = 0.0

           else

!              ! compute left edge
!              y_vort1 = jLeft(j) - sim_yctrVortex
!              ! mirror about the periodic boundary
!              y_vort2 = jLeft(j) - &
!                 ( sim_ymin - sim_yctrVortex + sim_ymax )
!              ! center of domain
!              z_vort = kLeft(k) - 0.5*(sim_zmax+sim_zmin)
!              vortex_stream1 = sim_vortexStrength * &
!                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!              ! counter-rotating, so make minus
!              vortex_stream2 = -sim_vortexStrength * &
!                 exp( -( y_vort2**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!
!              vy_l = alpha_inv * z_vort / sim_vortexSize**2 * &
!                 ( vortex_stream1 + vortex_stream2 )
!              vz_l = -alpha_inv / sim_vortexSize**2 * &
!                 ( y_vort1*vortex_stream1 + y_vort2*vortex_stream2 )
!
!              ! compute right edge
!              y_vort1 = jRight(j) - sim_yctrVortex
!              ! mirror about the periodic boundary
!              y_vort2 = jRight(j) - &
!                 ( sim_ymin - sim_yctrVortex + sim_ymax )
!              ! center of domain
!              z_vort = kRight(k) - 0.5*(sim_zmax+sim_zmin)
!              vortex_stream1 = sim_vortexStrength * &
!                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!              ! counter-rotating, so make minus
!              vortex_stream2 = -sim_vortexStrength * &
!                 exp( -( y_vort2**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!
!              vy_r = alpha_inv * z_vort / sim_vortexSize**2 * &
!                 ( vortex_stream1 + vortex_stream2 )
!              vz_r = -alpha_inv / sim_vortexSize**2 * &
!                 ( y_vort1*vortex_stream1 + y_vort2*vortex_stream2 )
!
!              ! compute cell-center
!              y_vort1 = jCenter(j) - sim_yctrVortex
!              ! mirror about the periodic boundary
!              y_vort2 = jCenter(j) - &
!                 ( sim_ymin - sim_yctrVortex + sim_ymax )
!              ! center of domain
!              z_vort = kCenter(k) - 0.5*(sim_zmax+sim_zmin)
!              vortex_stream1 = sim_vortexStrength * &
!                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!              ! counter-rotating, so make minus
!              vortex_stream2 = -sim_vortexStrength * &
!                 exp( -( y_vort2**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )
!
!              vy_c = alpha_inv * z_vort / sim_vortexSize**2 * &
!                 ( vortex_stream1 + vortex_stream2 )
!              vz_c = -alpha_inv / sim_vortexSize**2 * &
!                 ( y_vort1*vortex_stream1 + y_vort2*vortex_stream2 )

!**Only one vortex
              y_vort1 = jLeft(j) - sim_yctrVortex
              z_vort = kLeft(k) - 0.5*(sim_zmax+sim_zmin)
              vortex_stream1 = sim_vortexStrength * &
                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )

              ! vy = d/dz (stream), vz = -d/dy (stream)
              vy_l = -z_vort / sim_vortexSize**2 * vortex_stream1
              vz_l = y_vort1 / sim_vortexSize**2 * vortex_stream1

              y_vort1 = jRight(j) - sim_yctrVortex
              z_vort = kRight(k) - 0.5*(sim_zmax+sim_zmin)
              vortex_stream1 = sim_vortexStrength * &
                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )

              ! vy = d/dz (stream), vz = -d/dy (stream)
              vy_r = -z_vort / sim_vortexSize**2 * vortex_stream1
              vz_r = y_vort1 / sim_vortexSize**2 * vortex_stream1

              y_vort1 = jCenter(j) - sim_yctrVortex
              z_vort = kCenter(k) - 0.5*(sim_zmax+sim_zmin)
              vortex_stream1 = sim_vortexStrength * &
                 exp( -( y_vort1**2 + z_vort**2 ) / 2.0 / sim_vortexSize**2 )

              ! vy = d/dz (stream), vz = -d/dy (stream)
              vy_c = -z_vort / sim_vortexSize**2 * vortex_stream1
              vz_c = y_vort1 / sim_vortexSize**2 * vortex_stream1

              ! compute cell-average (simplified since density is constant in
              ! y-z plane)
              vely = (vy_l + 4.0*vy_c + vy_r) / 6.0
              velz = (vz_l + 4.0*vz_c + vz_r) / 6.0

           endif
           
           kine = 0.5 * ( velx**2 + vely**2 + velz**2 )

           !-----------------------------------------------
           !  Now store all this info on the grid
           !-----------------------------------------------
           cell(IAXIS) = i
           cell(JAXIS) = j
           cell(KAXIS) = k
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, cell, state(EOS_DENS))
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, cell, state(EOS_TEMP))

           call Grid_putPointData(blockId, CENTER, FLAM_MSCALAR, EXTERIOR, cell, flam)

           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, cell, velx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, cell, vely)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, cell, velz)

           !  usually I would just call the EOS, but we happen to have all this data
           !  so we'll just put it in.
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, cell, state(EOS_EINT)+kine)
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, cell, state(EOS_EINT))
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, cell, state(EOS_PRES))
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, cell, state(EOS_GAMC))
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, cell, &
                                       state(EOS_PRES)/(state(EOS_DENS)*state(EOS_EINT))+1.0)
        enddo
     enddo
  enddo
  
  deallocate(iCenter)
  deallocate(iLeft)
  deallocate(iRight)
  deallocate(jCenter)
  deallocate(jLeft)
  deallocate(jRight)
  deallocate(kCenter)
  deallocate(kLeft)
  deallocate(kRight)

  return
  
end subroutine Simulation_initBlock
