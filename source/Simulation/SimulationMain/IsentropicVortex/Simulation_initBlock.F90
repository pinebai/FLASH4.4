!!****if* source/Simulation/SimulationMain/IsentropicVortex/Simulation_initBlock
!!
!!  NAME
!!    Simulation_initBlock
!!
!!  SYNOPSIS
!!    call Simulation_initBlock( integer(IN) :: blockID,
!!                               
!!
!!  DESCRIPTION
!!    Initialize fluid data for a specified block.
!!
!!    This version initializes conditions for the isentropic_vortex
!!    problem. Given a constant background state specified by
!!    density rho_o, pressure p_o, and velocity (u_o, v_o).
!!    Define the following perturbations:
!!     du = -y (vortex_strength/ (2 pi)) exp ((1 - r^2)/2)
!!     dv =  x (vortex_strength/ (2 pi)) exp ((1 - r^2)/2)
!!     dT = - ((gamma-1)/2) (vortex_strength/ (2 pi))^2 exp (1 - r^2)
!!    where x and y are the distance from the vortex center, and
!!    r^2 = x^2 + y^2
!!
!!    Then the conserved variables at time=0 are
!!     rho   = (T_o + dT)^(1/(gamma-1))
!!     rho u = rho (u_o + du)
!!     rho v = rho (v_o + dv)
!!     rho e = (rho^gamma)/(gamma-1) + rho (u^2 + v^2)/2
!!    
!!  ARGUMENTS
!!    blockID   --    The number of the block to initialize.
!!
!!  PARAMETERS
!!    rho_ambient       \
!!    p_ambient          | Background density, pressure, x- and y-velocity
!!    u_ambient          |
!!    v_ambient         /
!!
!!    vortex_strength   Vortex strength (beta in the reference);
!!                        see description
!!
!!    xctr              \ The coordinates of the vortex center
!!    yctr              / 
!!
!!    nx_subint         \ For estimating initial cell averages, the
!!    ny_subint         /  number of subintervals in x- and y-directions
!!
!!  USES
!!   logfile
!!   physical_constants: get_constant_from_db
!!   multifluid_database
!!   runtime_parameters
!!   dBase: nxb, nyb, nzb, nguard, ionmax, k2d, k3d, ndim,
!!    iHi_gc, jHi_gc, kHi_gc, dBasePropertyInteger, dBaseKeyNumber,
!!    dBaseSpecies, dBaseGetCoords, dBasePutData
!! 
!!  NOTES
!!    Reference:  Yee, Vinokur & Djomehri, J. Comp. Phys 162 
!!
!!***
subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_gamma, sim_uAmbient, sim_vAmbient,&
       sim_vortexStrength, sim_xctrTrue, sim_yctrTrue, sim_nxSubint, &
       sim_nySubint,  sim_imax, sim_jmax, sim_imin, sim_jmin, &
       sim_imidDomain, sim_jmidDomain, sim_diDomain, sim_djDomain,&
       sim_tStarAmbient, sim_rbar, sim_constAmbient, sim_smlrho, sim_smallx,&
       sim_eosData, sim_eosMassFr

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putRowData, Grid_getRowData
  use Eos_interface, ONLY : Eos_wrapped, Eos

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer, intent(IN) :: blockID

  real :: xctr_closest, yctr_closest
  real :: temp_coeff, temp_exp, gm1i, rbari
  real ::  dx_loc, dy_loc, xp_loc, yp_loc, rp2_loc, deni
  real :: del_u, del_v, t_star, del_t_star
  real :: rho_sum, rhou_sum, rhov_sum, rhow_sum, rhoe_sum
  real :: rhou, rhov, rhow, rhoe
  real :: rho_loc, u_loc, v_loc, w_loc, t_loc, p_loc, e_loc, gamma_loc

! Needed for eos call.
  real :: abar, zbar, dpt, dpd, det, ded, c_v, c_p, pel, xxne, eta

  integer :: i, j, k, n
  integer :: ii, jj
  real :: entropy, dst, dsd
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC,eosRange
  integer,dimension(MDIM) :: startingPos
  integer :: sizeX,sizeY,sizeZ,vecLen=1
  logical :: gcell=.false.

  real,allocatable,dimension(:)::xCenter,xLeft,xRight
  real,allocatable,dimension(:)::yCenter,yLeft,yRight
  real,allocatable,dimension(:)::zCenter,zLeft,zRight
  real,allocatable,dimension(:)::rho, p, t, e, u, v, w, etot, game, gamc
#if NSPECIES > 0
  real,allocatable,dimension(:,:)::xn
#endif

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1  ! These are INTERIOR sizes
  sizeY=blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
  sizeZ=blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

  allocate(xCenter(sizeX))
  allocate(xRight(sizeX))
  allocate(xLeft(sizeX))
  allocate(yCenter(sizeY))
  allocate(yRight(sizeY))
  allocate(yLeft(sizeY))
  allocate(zCenter(sizeZ))
  allocate(zRight(sizeZ))
  allocate(zLeft(sizeZ))
  allocate(rho(sizeX))
  allocate(p(sizeX))
  allocate(t(sizeX))
  allocate(e(sizeX))
  allocate(u(sizeX))
  allocate(v(sizeX))
  allocate(w(sizeX))
  allocate(etot(sizeX))
  allocate(game(sizeX))
  allocate(gamc(sizeX))
#if NSPECIES > 0
  allocate(xn(sizeX,NSPECIES))

!   Assume a single species.
    xn(:,:) = sim_smallx
    xn(:,1) = 1.
#endif

  if (NDIM > 2) then
     call Grid_getCellCoords (KAXIS, blockId, CENTER, gcell, zCenter, sizeZ)
     call Grid_getCellCoords (KAXIS, blockId, LEFT_EDGE, gcell, zLeft, sizeZ)
     call Grid_getCellCoords (KAXIS, blockId, RIGHT_EDGE, gcell, zRight, sizeZ)
  end if
  if(NDIM>1) then
     call Grid_getCellCoords (JAXIS, blockId, CENTER, gcell, yCenter, sizeY)
     call Grid_getCellCoords (JAXIS, blockId, LEFT_EDGE, gcell, yLeft, sizeY)
     call Grid_getCellCoords (JAXIS, blockId, RIGHT_EDGE, gcell, yRight, sizeY)
  end if
  call Grid_getCellCoords (IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords (IAXIS, blockId, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords (IAXIS, blockId, RIGHT_EDGE, gcell, xRight, sizeX)
! Initialize the flowfield.

  temp_coeff = sim_vortexStrength/(2.0*PI)
  gm1i = 1.0/(sim_gamma-1.0)
  rbari = 1.0/sim_rbar

  do k = 1,sizeZ
     do j = 1,sizeY
        do i = 1,sizeX
           
           !     these are cell sizes and cell subinterval sizes.
           dx_loc = (xRight(i) - xLeft(i))/sim_nxSubint
           dy_loc = (yRight(j) - yLeft(j))/sim_nySubint
           
           ! find the 'nearest' vortex - might be at the specified
           ! center (sim_xctrTrue,sim_yctrTrue) or one of its periodic images
           ! the search is based on which quadrant of the (square) domain
           ! the point of interest (xCenter(i), yCenter(j)) is in
           if (xCenter(i) < sim_imidDomain) then    ! image left of 
              !xmin might be closer
              if ( abs(xCenter(i) -  sim_xctrTrue)             <     &
                   abs(xCenter(i) - (sim_xctrTrue - sim_diDomain)) ) then 
                   ! true is closer
                 xctr_closest = sim_xctrTrue
              else                      ! image is closer
                 xctr_closest = sim_xctrTrue - sim_diDomain
              endif
           else                      ! image right of sim_imax might be closer
              if ( abs(xCenter(i) -  sim_xctrTrue)             <     &
                   abs(xCenter(i) - (sim_xctrTrue + sim_diDomain)) ) then  
                 ! true is closer
                 xctr_closest = sim_xctrTrue
              else                                ! image is closer
                 xctr_closest = sim_xctrTrue + sim_diDomain
              endif
           endif
           if (yCenter(j) < sim_jmidDomain) then    ! image below 
              ! ymin might be closer
              if ( abs(yCenter(j) -  sim_yctrTrue)             <     &
                   abs(yCenter(j) - (sim_yctrTrue - sim_djDomain)) ) then  
                 ! true is closer
                 yctr_closest = sim_yctrTrue
              else                    ! image is closer
                 yctr_closest = sim_yctrTrue - sim_djDomain
              endif
           else                     ! image above ymax might be closer
              if ( abs(yCenter(j) -  sim_yctrTrue)             <     &
                   abs(yCenter(j) - (sim_yctrTrue + sim_djDomain)) ) then  
                 ! true is closer
                 yctr_closest = sim_yctrTrue
              else                              ! image is closer
                 yctr_closest = sim_yctrTrue + sim_djDomain
              endif
           endif
           
           !   Obtain cell averages by numerical integration (summation).
           !   Initialize cell sums to zero.
           rho_sum  = 0.0
           rhou_sum = 0.0
           rhov_sum = 0.0
           rhow_sum = 0.0
           rhoe_sum = 0.0
           
           !  Each cell is subdivided into sim_nxSubint * sim_nySubint 
           !   subcells.

           do jj = 1, sim_nySubint
              do ii = 1, sim_nxSubint

                 
                 xp_loc = xLeft(i) + (float(ii-1) + 0.5)*dx_loc - xctr_closest
                 yp_loc = yLeft(j) + (float(jj-1) + 0.5)*dy_loc - yctr_closest
                 rp2_loc = xp_loc**2 + yp_loc**2
                 
                 temp_exp = exp(0.5*(1.0-rp2_loc))
                 
                 del_u = -yp_loc*temp_coeff*temp_exp
                 del_v =  xp_loc*temp_coeff*temp_exp
                 del_t_star = -((sim_gamma - 1.0)/sim_gamma)*0.5*&
                               (temp_coeff*temp_exp)**2
                 t_star = sim_tStarAmbient + del_t_star

                 u_loc = sim_uAmbient + del_u
                 v_loc = sim_vAmbient + del_v
                 w_loc = 0.0
                 t_loc = t_star*rbari
                 rho_loc = sim_constAmbient*t_loc**gm1i
                 sim_eosData(EOS_DENS)=rho_loc
                 sim_eosData(EOS_TEMP)= t_loc
                 !       Get e_loc from rho_loc and t_loc.

                 call Eos(MODE_DENS_TEMP,vecLen,sim_eosData,sim_eosMassFr)


                 e_loc=sim_eosData(EOS_EINT)
                 p_loc=sim_eosData(EOS_PRES)

                 rho_sum  = rho_sum  + rho_loc
                 rhou_sum = rhou_sum + rho_loc*u_loc
                 rhov_sum = rhov_sum + rho_loc*v_loc
                 rhow_sum = rhow_sum + rho_loc*w_loc
                 rhoe_sum = rhoe_sum + rho_loc*(e_loc +  &
                               0.5*(u_loc**2 + v_loc**2 + w_loc**2)) 
              enddo ! ii-loop
           enddo ! jj-loop

           deni = 1.0/float(sim_nxSubint*sim_nySubint)
           rho(i) = rho_sum *deni
           rhou   = rhou_sum*deni
           rhov   = rhov_sum*deni
           rhow   = rhow_sum*deni
           rhoe   = rhoe_sum*deni
           
           u(i)    = rhou/rho(i)
           v(i)    = rhov/rho(i)
           w(i)    = rhow/rho(i)
           etot(i) = rhoe/rho(i)
           e(i)    = etot(i) - 0.5*(u(i)**2 + v(i)**2 + w(i)**2)
           
        enddo  ! i-loop

       !this routine initializes a row at a time...
        startingPos(IAXIS) = 1
        startingPos(JAXIS) = j
        startingPos(KAXIS) = k
#if NSPECIES > 0
        do n=0,NSPECIES-1
           call Grid_putRowData(blockID,CENTER,SPECIES_BEGIN+n,INTERIOR,IAXIS,&
                                startingPos,xn(:,n+1),sizeX)
        end do
#endif
 
        call Grid_putRowData(blockId, CENTER, DENS_VAR, INTERIOR, &
                             IAXIS, startingPos, rho, sizeX)
        call Grid_putRowData(blockId, CENTER, VELX_VAR, INTERIOR, &
                             IAXIS, startingPos, u, sizeX)
        call Grid_putRowData(blockId, CENTER, VELY_VAR, INTERIOR, &
                             IAXIS, startingPos, v, sizeX)
        call Grid_putRowData(blockId, CENTER, VELZ_VAR, INTERIOR, &
                             IAXIS, startingPos, w, sizeX)
        call Grid_putRowData(blockId, CENTER, ENER_VAR, INTERIOR, &
                             IAXIS, startingPos, etot, sizeX)
        call Grid_putRowData(blockID,CENTER,EINT_VAR,INTERIOR,&
                             IAXIS, startingPos, e, sizeX)

!   Get p, t from rho, e
!   Since these are INTERIOR calls, the range should not include the guard cells
!   With INTERIOR usage, the index 1 == the INTERIOR edge.  However, Eos_wrapped
!   can't really work with the INTERIOR mode, it has to collect ALL the data and
!   must be indexed with EXTERIOR mode.
!   So we need to offset j,k by the guard cell index.
        eosRange(LOW,IAXIS)=blkLimits(LOW,IAXIS)
        eosRange(HIGH,IAXIS)=blkLimits(HIGH,IAXIS)
        eosRange(:,JAXIS) = blkLimits(LOW,JAXIS) + j-1
        eosRange(:,KAXIS) = blkLimits(LOW,KAXIS) + k-1
        call Eos_wrapped(MODE_DENS_EI,eosRange,blockID)

        call Grid_getRowData(blockID,CENTER,PRES_VAR,INTERIOR,&
                             IAXIS, startingPos, p, sizeX)
        game = p/(rho*e) + 1.0
        
        call Grid_putRowData(blockId, CENTER, GAME_VAR, INTERIOR, &
                             IAXIS, startingPos, game, sizeX)


     enddo  ! j-loop
  enddo  ! k-loop
  
  deallocate(xCenter)
  deallocate(xRight)
  deallocate(xLeft)
  deallocate(yCenter)
  deallocate(yRight)
  deallocate(yLeft)
  deallocate(zCenter)
  deallocate(zRight)
  deallocate(zLeft)
  deallocate(rho)
  deallocate(p)
  deallocate(t)
  deallocate(e)
  deallocate(u)
  deallocate(v)
  deallocate(w)
  deallocate(etot)
  deallocate(game)
  deallocate(gamc)
#if NSPECIES > 0
  deallocate(xn)
#endif

  return
end
