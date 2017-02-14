!!****if* source/Simulation/SimulationMain/NeiTest/hy_ppm_updateSoln
!! NAME
!!
!!  hy_ppm_updateSoln
!!
!!
!! SYNOPSIS
!!
!!  hy_ppm_updateSoln(integer, intent(IN)  :: rangeSwitch, 
!!                integer, intent(IN)  :: xyzswp, 
!!                real, intent(IN)     :: dt,           
!!                integer, intent(IN)  :: blkLimits(HIGH,MDIM),
!!                integer, intent(IN)  :: blkLimitsGC(HIGH,MDIM),
!!                integer, intent(IN)  :: numCells,  
!!                real, intent(IN)     :: tempArea(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, intent(IN)     :: tempGrav1d_o(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                     blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                     blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, intent(IN)     :: tempGrav1d(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                   blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                   blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, intent(IN)     :: tempDtDx(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, intent(IN)     :: tempFict(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, intent(IN)     :: tempFlx(NFLUXES, &
!!                                                 blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),
!!                real, pointer,       :: solnData(:,:,:,:)
!!
!! DESCRIPTION
!!
!!  Update the cell average quantities in time based on the fluxes through the 
!!  boundaries.  This is a general update routine for conservation laws in a 
!!  directionally split, finite-volume method.  Given the time-averaged, 
!!  cell-edge fluxes (tempFlx,...), update the solution vector to the new
!!  time in the direction specified by xyzswp.
!!
!!  rangeSwitch specifies which cells to update.  When performing flux 
!!  conservation, it is sometimes useful to delay the update at the boundaries
!!  until after the corrected fluxes have been computed.  rangeSwitch lets
!!  you update either all the cells, all the cells away from the boundary,
!!  or only the cells abutting a boundary.  
!!
!! ARGUMENTS
!!
!!  rangeSwitch --   If UPDATE_INTERIOR, update only "interior" cells. That means here,
!!                    all cells except the first and last layer in the sweep direction.
!!                   If UPDATE_BOUND, update only the first and the last layer of cells
!!                    in the sweep direction.
!!                   Otherwise update all cells.
!!                   Guard cells are never updated.
!!
!!  xyzswp --        The direction of the current sweep, one of SWEEP_X, SWEEP_Y, SWEEP_Z
!!
!!  dt --            The timestep to advance the solution by
!!  blkLimits  -- Index limits of the block interior
!!  blkLimitsGC -- Index limits of the block including guardcell
!!  numCells  -- number of cells in one dimension
!!
!!  tempArea --      Temp data from hydro_1d
!!  tempGrav1d_o -
!!  tempGrav1d -
!!  tempDtDx -       The timestep divided by the cell volume (with 
!!                   geometrical factors)
!!  tempFict -       Geometry related forces, e.g., centrifugal; 
!!                   used to update velocities 
!!
!!  tempFlx --       The fluxes through the boundary
!!
!!  solnData --      Pointer to a block of data
!!
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


#define CIP_
#ifdef DEBUG_ALL
#define DEBUG_HYDRO
#endif

subroutine hy_ppm_updateSoln(rangeSwitch,                        &
                         xyzswp, dt,                          &
                         blkLimits,blkLimitsGC,numCells,               &
                         tempArea, tempGrav1d_o, tempGrav1d,  &
                         tempDtDx, tempFict,                  &
                         tempFlx,  solnData )
!==============================================================================

  use Hydro_data, ONLY: hy_numXn
  use Hydro_data, ONLY: hy_smlrho, hy_smallp, hy_eintSwitch, hy_useCmaAdvection

  implicit none  
#include "constants.h"
#include "PPM.h"
#include "Flash.h"

  integer, intent(IN) :: rangeSwitch
  integer, intent(IN) :: xyzswp
  real,    intent(IN) :: dt
  integer,intent(IN) :: numCells  
  integer, intent(IN),dimension(2,MDIM)::blkLimitsGC,blkLimits
#ifdef FIXEDBLOCKSIZE
  real, intent(IN), DIMENSION(GRID_ILO_GC:GRID_IHI_GC, &
                              GRID_JLO_GC:GRID_JHI_GC, &
                              GRID_KLO_GC:GRID_KHI_GC  ) :: &
                                                  tempArea, tempGrav1d_o, &
                                                  tempGrav1d, &
                                                  tempDtDx, tempFict
  real, intent(IN), DIMENSION(NFLUXES,GRID_ILO_GC:GRID_IHI_GC, &
                        GRID_JLO_GC:GRID_JHI_GC, &
                        GRID_KLO_GC:GRID_KHI_GC) :: tempFlx
#else
  real, intent(IN), DIMENSION(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                              blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                              blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)  ) :: &
                                                  tempArea, tempGrav1d_o, &
                                                  tempGrav1d, &
                                                  tempDtDx, tempFict
  real, intent(IN), DIMENSION(NFLUXES,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                        blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                        blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: tempFlx

#endif

  real, pointer ::  solnData(:,:,:,:) 

  integer :: i, j, k, kk, n
  integer :: imin, imax, iskip, jmin, jmax, jskip, kmin, kmax, kskip
  real    :: xnflx2(hy_numXn), xnflx1(hy_numXn), dtdx(numCells)

  real    :: rhoflx1, rhoflx2, uflx1, uflx2, pav1, pav2, utflx1, &
             utflx2, uttflx1, uttflx2, eflx1, eflx2, eintflx1, eintflx2

  real    :: aold_t_geomfc, grav1d_o, grav1d, fict1d, rho_o, velx_o
  real    :: aux1, hunu, hutnu, huttnu, henu, heintnu, ekin,     &
             einternal, ekin_old, etot
  real    :: hxnnu, vely_o, velz_o, inv_new_dens, inv_dens


!!$  integer, save :: irhoflx,iuflx,ipflx,iutflx,iuttflx,           &
!!$                   ieflx,ieintflx,inucflx_begin,                 &
!!$                   idens,ivelx,ively,ivelz,ipres,iener,          &
!!$                   itemp,igamc,igame,inuc_begin, ieint
!!$
!!$  logical, save :: firstcall = .true.



! update in x direction
  
  if (xyzswp == SWEEP_X) then
     
     select case (rangeSwitch)
     case (UPDATE_INTERIOR)
        imin  = blkLimits(LOW,IAXIS)+ 1
        imax  = blkLimits(HIGH,IAXIS) - 1
        iskip = 1
     case (UPDATE_BOUND)
        imin  = blkLimits(LOW,IAXIS)
        imax  = blkLimits(HIGH,IAXIS)
        iskip = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)
     case default
        imin  = blkLimits(LOW,IAXIS)
        imax  = blkLimits(HIGH,IAXIS)
        iskip = 1
     end select
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = imin, imax, iskip
              
              dtdx(i)  = tempDtDx(i,j,k)
              rhoflx1  = tempFlx(RHO_FLUX,i,j,k)
              rhoflx2  = tempFlx(RHO_FLUX,i+1,j,k)
              
!!$              uflx1    = tempFlx(U_FLUX,i,j,k)
!!$              uflx2    = tempFlx(U_FLUX,i+1,j,k)
              
!!$              pav1     = tempFlx(P_FLUX,i,j,k)
!!$              pav2     = tempFlx(P_FLUX,i+1,j,k)
              
!!$              utflx1   = tempFlx(UT_FLUX,i,j,k)
!!$              utflx2   = tempFlx(UT_FLUX,i+1,j,k)
              
!!$              uttflx1  = tempFlx(UTT_FLUX,i,j,k)
!!$              uttflx2  = tempFlx(UTT_FLUX,i+1,j,k)
              
!!$              eflx1    = tempFlx(E_FLUX,i,j,k)
!!$              eflx2    = tempFlx(E_FLUX,i+1,j,k)
!!$
!!$              eintflx1 = tempFlx(EINT_FLUX,i,j,k)
!!$              eintflx2 = tempFlx(EINT_FLUX,i+1,j,k)

              do kk = 1,hy_numXn
                 xnflx1(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i,j,k)
                 xnflx2(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i+1,j,k)
              end do

!!$              aold_t_geomfc = tempArea(i,j,k)
!!$              grav1d_o  = tempGrav1d_o(i,j,k)
!!$              grav1d    = tempGrav1d(i,j,k)
!!$              fict1d    = tempFict(i,j,k)
              
              rho_o     = solnData(DENS_VAR,i,j,k)
!!$              velx_o    = solnData(VELX_VAR,i,j,k)

! update the density               
              solnData(DENS_VAR,i,j,k) = solnData(DENS_VAR,i,j,k) -  &
                   dtdx(i) * (rhoflx2-rhoflx1)
              
! update the mass fractions
              inv_dens = 1.e0/solnData(DENS_VAR,i,j,k)
              do n = 1, NSPECIES
                 aux1  = - dtdx(i) * (xnflx2(n) - xnflx1(n))
                 hxnnu = rho_o*solnData(SPECIES_BEGIN-1+n,i,j,k) + aux1
                 solnData(SPECIES_BEGIN-1+n,i,j,k) = hxnnu*inv_dens
              end do

! update the passive scalars
              do n = NSPECIES+1, hy_numXn                           
                 aux1  = - dtdx(i) * (xnflx2(n) - xnflx1(n))
                 hxnnu = solnData(SPECIES_BEGIN-1+n,i,j,k) + aux1
                 solnData(SPECIES_BEGIN-1+n,i,j,k) = hxnnu
              end do    
              
              solnData(DENS_VAR,i,j,k) =                              &
                   max(hy_smlrho, solnData(DENS_VAR,i,j,k))

              ! NO OTHER VARIABLES TO UPDATE IN THIS SPECIAL IMPLEMENTATION FOR NEITEST
!!$              inv_new_dens = 1.e0/solnData(DENS_VAR,i,j,k)
!!$
!!$! first store the old kinetic energy
!!$              ekin_old = 0.5e00 * ((solnData(VELX_VAR,i,j,k))**2 +    &
!!$                   (solnData(VELY_VAR,i,j,k))**2 +                    &
!!$                   (solnData(VELZ_VAR,i,j,k))**2)
!!$
!!$! update the velocities               
!!$              aux1 = - dtdx(i) * ((uflx2 - uflx1) +                  &
!!$                   aold_t_geomfc * (pav2 - pav1))                    &
!!$               +0.5e0*dt                                             &
!!$               *( (solnData(DENS_VAR,i,j,k) + rho_o)*fict1d           &
!!$                 +(rho_o*grav1d_o + solnData(DENS_VAR,i,j,k)*grav1d)  &
!!$                )
!!$
!!$              hunu  = (rho_o * solnData(VELX_VAR,i,j,k) + aux1)
!!$!              solnData(VELX_VAR,i,j,k) = hunu*inv_new_dens
!!$              
!!$              aux1  = - dtdx(i) * (utflx2 - utflx1)
!!$              hutnu = (rho_o * solnData(VELY_VAR,i,j,k) + aux1)
!!$!              solnData(VELY_VAR,i,j,k) = hutnu*inv_new_dens
!!$              
!!$              aux1   = - dtdx(i) * (uttflx2 - uttflx1)
!!$              huttnu = (rho_o * solnData(VELZ_VAR,i,j,k) + aux1)
!!$!              solnData(VELZ_VAR,i,j,k) = huttnu*inv_new_dens
!!$
!!$! get the internal energy
!!$              einternal = solnData(EINT_VAR,i,j,k)
!!$
!!$! update the total energy
!!$              aux1 = - dtdx(i) * (eflx2 - eflx1)                             &
!!$                   + dt*0.5e00                                               &
!!$                   *( rho_o*velx_o*grav1d_o                                  &
!!$                     +solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)*grav1d &
!!$                    )
!!$
!!$              henu = (rho_o * solnData(ENER_VAR,i,j,k) + aux1)
!!$              etot = henu*inv_new_dens
!!$
!!$! update internal energy
!!$              aux1 = - dtdx(i) * ( (eintflx2 - eintflx1)              &
!!$                   - 0.5e0*(velx_o + solnData(VELX_VAR,i,j,k))*        &
!!$                   aold_t_geomfc * (pav2 - pav1) )
!!$              
!!$              heintnu   = rho_o * einternal + aux1
!!$              einternal = heintnu*inv_new_dens
!!$              einternal = max(einternal, hy_smallp*inv_new_dens)
!!$
!!$! compute the new kinetic energy       
!!$              ekin = 0.5e00 * ((solnData(VELX_VAR,i,j,k))**2 +         &
!!$                     (solnData(VELY_VAR,i,j,k))**2 +                   &
!!$                     (solnData(VELZ_VAR,i,j,k))**2)
!!$
!!$! test whether we should use the internal energy from the evolution
!!$              if (einternal .LT. hy_eintSwitch*ekin) then
!!$                 etot = einternal + ekin
!!$              else
!!$                 einternal = etot - ekin
!!$              endif
!!$              
!!$!              solnData(ENER_VAR,i,j,k) = etot
!!$!              solnData(EINT_VAR,i,j,k) = einternal
              
           end do
        end do
     end do
     
  end if
               
!===============================================================================
#if NDIM > 1

! update in y direction
  
  if ((NDIM >= 2) .and. (xyzswp == SWEEP_Y)) then

     select case (rangeSwitch)
     case (UPDATE_INTERIOR)
        jmin  = blkLimits(LOW,JAXIS) + 1 ! NGUARD*K2D + 2
        jmax  = blkLimits(HIGH,JAXIS) - 1 ! NGUARD*K2D + NYB - 1
        jskip = 1
     case (UPDATE_BOUND)
        jmin  = blkLimits(LOW,JAXIS) ! NGUARD*K2D + 1
        jmax  = blkLimits(HIGH,JAXIS) ! NGUARD*K2D + NYB
        jskip = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) !NYB - 1
     case default
        jmin  = blkLimits(LOW,JAXIS) ! NGUARD*K2D + 1
        jmax  = blkLimits(HIGH,JAXIS) ! NGUARD*K2D + NYB
        jskip = 1
     end select

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = jmin, jmax, jskip
           do i = blkLimits(LOW,IAXIS) ,blkLimits(HIGH,IAXIS)

              dtdx(j) = tempDtDx(i,j,k)

              rhoflx1  = tempFlx(RHO_FLUX,i,j,k)
              rhoflx2  = tempFlx(RHO_FLUX,i,j+1,k)

              uflx1    = tempFlx(U_FLUX,i,j,k)
              uflx2    = tempFlx(U_FLUX,i,j+1,k)

              pav1     = tempFlx(P_FLUX,i,j,k)
              pav2     = tempFlx(P_FLUX,i,j+1,k)

              utflx1   = tempFlx(UT_FLUX,i,j,k)
              utflx2   = tempFlx(UT_FLUX,i,j+1,k)

              uttflx1  = tempFlx(UTT_FLUX,i,j,k)
              uttflx2  = tempFlx(UTT_FLUX,i,j+1,k)

              eflx1    = tempFlx(E_FLUX,i,j,k)
              eflx2    = tempFlx(E_FLUX,i,j+1,k)

              eintflx1 = tempFlx(EINT_FLUX,i,j,k)
              eintflx2 = tempFlx(EINT_FLUX,i,j+1,k)

              do kk = 1,hy_numXn
                 xnflx1(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i,j,k)
                 xnflx2(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i,j+1,k)
              end do


              aold_t_geomfc = tempArea(i,j,k)
              grav1d_o      = tempGrav1d_o(i,j,k)
              grav1d        = tempGrav1d(i,j,k)
              fict1d        = tempFict(i,j,k)

              rho_o  = solnData(DENS_VAR,i,j,k)
              vely_o = solnData(VELY_VAR,i,j,k)

! update the density
              solnData(DENS_VAR,i,j,k) = solnData(DENS_VAR,i,j,k) -  &
                   dtdx(j) * (rhoflx2-rhoflx1)

! update the mass fractions
              inv_dens = 1.e0/solnData(DENS_VAR,i,j,k)
              do n = 1, NSPECIES
                 aux1  = - dtdx(j) * (xnflx2(n) - xnflx1(n))
                 hxnnu = rho_o*solnData(SPECIES_BEGIN-1+n,i,j,k) + aux1
                 solnData(SPECIES_BEGIN-1+n,i,j,k) = hxnnu*inv_dens
              end do

! update the passive scalars
              do n = NSPECIES+1, hy_numXn
                 aux1  = - dtdx(j) * (xnflx2(n) - xnflx1(n))
                 hxnnu = solnData(SPECIES_BEGIN-1+n,i,j,k) + aux1
                 solnData(SPECIES_BEGIN-1+n,i,j,k) = hxnnu
              end do
 
              solnData(DENS_VAR,i,j,k) =                              &
                   max(hy_smlrho, solnData(DENS_VAR,i,j,k))

              inv_new_dens = 1.e0/solnData(DENS_VAR,i,j,k)

! first store the old kinetic energy
              ekin_old = 0.5e00 * ((solnData(VELX_VAR,i,j,k))**2 +    &
                   (solnData(VELY_VAR,i,j,k))**2 +                    &
                   (solnData(VELZ_VAR,i,j,k))**2)
              
! update the velocities               
              aux1 = - dtdx(j) * ((uflx2 - uflx1) +                  &
                   aold_t_geomfc * (pav2 - pav1))                    &
               +dt*0.5e00                                            &
               *( (solnData(DENS_VAR,i,j,k) + rho_o)*fict1d           &
                 +(rho_o*grav1d_o + solnData(DENS_VAR,i,j,k)*grav1d)  &
                )

              hunu  = (rho_o * solnData(VELY_VAR,i,j,k) + aux1)
!              solnData(VELY_VAR,i,j,k) = hunu*inv_new_dens

              aux1  = - dtdx(j) * (utflx2 - utflx1)
              hutnu = (rho_o * solnData(VELX_VAR,i,j,k) + aux1)
!              solnData(VELX_VAR,i,j,k) = hutnu*inv_new_dens

              aux1   = - dtdx(j) * (uttflx2 - uttflx1)
              huttnu = (rho_o * solnData(VELZ_VAR,i,j,k) + aux1)
!              solnData(VELZ_VAR,i,j,k) = huttnu*inv_new_dens

! get the internal energy
              einternal = solnData(EINT_VAR,i,j,k)
               
! update the total energy
              aux1 = - dtdx(j) * (eflx2 - eflx1)                             &
                   + dt*0.5e00                                               &
                   *( rho_o*vely_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)*grav1d &
                    )

              henu = (rho_o * solnData(ENER_VAR,i,j,k) + aux1)
              etot = henu*inv_new_dens

! update internal energy
              aux1 = - dtdx(j) * ( (eintflx2 - eintflx1)              &
                   - 0.5e0*(vely_o + solnData(VELY_VAR,i,j,k))*        &
                   aold_t_geomfc * (pav2 - pav1) )

              heintnu   = rho_o * einternal + aux1
              einternal = heintnu*inv_new_dens
              einternal = max(einternal, hy_smallp*inv_new_dens)

! compute the new kinetic energy
              ekin = 0.5e00 * ((solnData(VELX_VAR,i,j,k))**2 +         &
                     (solnData(VELY_VAR,i,j,k))**2 +                   &
                     (solnData(VELZ_VAR,i,j,k))**2)

! test whether we should use the internal energy from the evolution
              if (einternal .LT. hy_eintSwitch*ekin) then
                 etot = einternal + ekin
              else
                 einternal = etot - ekin
              endif
              
!              solnData(ENER_VAR,i,j,k) = etot
!              solnData(EINT_VAR,i,j,k) = einternal
              
              
           end do
        end do
     end do
     
  end if

!===============================================================================

! update in z direction

  if ((NDIM == 3) .and. (xyzswp == SWEEP_Z)) then
     
     select case (rangeSwitch)
     case (UPDATE_INTERIOR)
        kmin  = blkLimits(LOW,KAXIS) + 1 !NGUARD*K3D + 2
        kmax  = blkLimits(HIGH,KAXIS) - 1 ! NGUARD*K3D + NZB - 1
        kskip = 1
     case (UPDATE_BOUND)
        kmin  = blkLimits(LOW,KAXIS) ! NGUARD*K3D + 1
        kmax  = blkLimits(HIGH,KAXIS) ! NGUARD*K3D + NZB
        kskip = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS) ! NZB - 1
     case default
        kmin  = blkLimits(LOW,KAXIS) ! NGUARD*K3D + 1
        kmax  = blkLimits(HIGH,KAXIS) ! NGUARD*K3D + NZB
        kskip = 1
     end select

     do k = kmin, kmax, kskip
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS) ! NGUARD*K2D+1, NGUARD*K2D+NYB
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS) ! NGUARD+1, NGUARD+NXB

              dtdx(k)  = tempDtDx(i,j,k)

              rhoflx1  = tempFlx(RHO_FLUX,i,j,k)
              rhoflx2  = tempFlx(RHO_FLUX,i,j,k+1)

              uflx1    = tempFlx(U_FLUX,i,j,k)
              uflx2    = tempFlx(U_FLUX,i,j,k+1)

              pav1     = tempFlx(P_FLUX,i,j,k)
              pav2     = tempFlx(P_FLUX,i,j,k+1)

              utflx1   = tempFlx(UT_FLUX,i,j,k)
              utflx2   = tempFlx(UT_FLUX,i,j,k+1)

              uttflx1  = tempFlx(UTT_FLUX,i,j,k)
              uttflx2  = tempFlx(UTT_FLUX,i,j,k+1)

              eflx1    = tempFlx(E_FLUX,i,j,k)
              eflx2    = tempFlx(E_FLUX,i,j,k+1)

              eintflx1 = tempFlx(EINT_FLUX,i,j,k)
              eintflx2 = tempFlx(EINT_FLUX,i,j,k+1)

              do kk = 1,hy_numXn
                 xnflx1(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i,j,k)
                 xnflx2(kk) = tempFlx(SPECIES_FLUX_BEGIN+kk-1,i,j,k+1)
              end do

              aold_t_geomfc = tempArea(i,j,k)
              grav1d_o      = tempGrav1d_o(i,j,k)
              grav1d        = tempGrav1d(i,j,k)
              fict1d        = tempFict(i,j,k)

              rho_o  = solnData(DENS_VAR,i,j,k)
              velz_o = solnData(VELZ_VAR,i,j,k)

! update the density
              solnData(DENS_VAR,i,j,k) = solnData(DENS_VAR,i,j,k) -  &
                   dtdx(k) * (rhoflx2-rhoflx1)

! update the mass fractions
              inv_dens = 1.e0/solnData(DENS_VAR,i,j,k)
              do n = 1, NSPECIES
                 aux1  = - dtdx(k) * (xnflx2(n) - xnflx1(n))
                 hxnnu = rho_o*solnData(SPECIES_BEGIN-1+n,i,j,k) + aux1
                 solnData(SPECIES_BEGIN-1+n,i,j,k) = hxnnu*inv_dens
              end do

! update the passive scalars
              do n = NSPECIES+1, hy_numXn
                 aux1  = - dtdx(k) * (xnflx2(n) - xnflx1(n))
                 hxnnu = solnData(SPECIES_BEGIN-1+n,i,j,k) + aux1
                 solnData(SPECIES_BEGIN-1+n,i,j,k) = hxnnu
              end do 

              solnData(DENS_VAR,i,j,k) =                              &
                   max(hy_smlrho, solnData(DENS_VAR,i,j,k))

              inv_new_dens = 1.e0/solnData(DENS_VAR,i,j,k)

! first store the old kinetic energy
              ekin_old = 0.5e00 * ((solnData(VELX_VAR,i,j,k))**2 +    &
                   (solnData(VELY_VAR,i,j,k))**2 +                    &
                   (solnData(VELZ_VAR,i,j,k))**2)
 
! update the velocities
              aux1 = - dtdx(k) * ((uflx2 - uflx1) +                  &
                   aold_t_geomfc * (pav2 - pav1))                    &
               +dt*0.5e00                                            &
               *( (solnData(DENS_VAR,i,j,k) + rho_o)*fict1d           &
                 +(rho_o*grav1d_o + solnData(DENS_VAR,i,j,k)*grav1d)  &
                )

              hunu  = (rho_o * solnData(VELZ_VAR,i,j,k) + aux1)
!              solnData(VELZ_VAR,i,j,k) = hunu*inv_new_dens
              
              aux1  = - dtdx(k) * (utflx2 - utflx1)
              hutnu = (rho_o * solnData(VELX_VAR,i,j,k) + aux1)
!              solnData(VELX_VAR,i,j,k) = hutnu*inv_new_dens
               
              aux1   = - dtdx(k) * (uttflx2 - uttflx1)
              huttnu = (rho_o * solnData(VELY_VAR,i,j,k) + aux1)
!              solnData(VELY_VAR,i,j,k) = huttnu*inv_new_dens

! get the internal energy
              einternal = solnData(EINT_VAR,i,j,k)

! update the total energy
              aux1 = - dtdx(k) * (eflx2 - eflx1)                             &
                   + dt*0.5e00                                               &
                   *( rho_o*velz_o*grav1d_o                                  &
                     +solnData(DENS_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)*grav1d &
                    )

              henu = (rho_o * solnData(ENER_VAR,i,j,k) + aux1)
              etot = henu*inv_new_dens
               
! update internal energy
              aux1 = - dtdx(k) * ( (eintflx2 - eintflx1)              &
                   - 0.5e0*(velz_o + solnData(VELZ_VAR,i,j,k))*        &
                   aold_t_geomfc * (pav2 - pav1) )
               
              heintnu   = rho_o * einternal + aux1
              einternal = heintnu*inv_new_dens
              einternal = max(einternal, hy_smallp*inv_new_dens)
               
! compute the new kinetic energy
              ekin = 0.5e00 * ((solnData(VELX_VAR,i,j,k))**2 +         &
                     (solnData(VELY_VAR,i,j,k))**2 +                   &
                     (solnData(VELZ_VAR,i,j,k))**2)
               
! test whether we should use the internal energy from the evolution
              if (einternal .LT. hy_eintSwitch*ekin) then
                 etot = einternal + ekin
              else
                 einternal = etot - ekin
              endif
              
!              solnData(ENER_VAR,i,j,k) = etot
!              solnData(EINT_VAR,i,j,k) = einternal
               
               
            end do
         end do
      end do
      
   end if
#endif
   
!===============================================================================
 end subroutine hy_ppm_updateSoln


