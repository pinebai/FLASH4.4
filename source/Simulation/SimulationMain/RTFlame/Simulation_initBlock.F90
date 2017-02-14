!!****if* source/Simulation/SimulationMain/RTFlame/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(integer(in) :: blockid)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   blockid : ID of block in current processor
!!
!! AUTOGENROBODOC
!!
!!
!!***

! Dean Townsley 2009
!
! init block for R-T Flame in a channel

subroutine Simulation_initBlock(blockID)
  
  use Simulation_data
  use Flame_interface, ONLY : Flame_getWidth
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkBoundBox, Grid_getDeltas, Grid_putPointData, &
    Grid_getCellCoords, Grid_getBlkCornerID, Grid_getMinCellSize
  use Eos_interface, ONLY : Eos
  use hse_interface, ONLY : flame_hse

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  
  integer, intent(in) :: blockID

  integer :: i, j, k, n

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: cornerID, stride, cell

  real :: dens, temp, ye, sumy, flam, velx, vely, xv, yv

  real, dimension(EOS_NUM) :: eosData

  real, allocatable, dimension(:) :: iCoords, jCoords, kCoords
  integer :: isizeGC, jsizeGC, ksizeGC
  real :: min_dx

  real :: x_f, fi_x_start, fi_x_end, flamewidth, xp, w, fi_buf_size
  integer :: fi_start, fi_end, fi_size, gi

  real, allocatable, dimension(:) :: fi_flam, fi_dens, fi_temp, fi_ye, fi_sumy

!==============================================================================

  ! get essential info about this block
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  call Grid_getBlkCornerID(BlockID,cornerID,stride)

  isizeGC = blkLimitsGC(HIGH,IAXIS)
  allocate(iCoords(isizeGC))
  jsizeGC = blkLimitsGC(HIGH,JAXIS)
  allocate(jCoords(jsizeGC))
  ksizeGC = blkLimitsGC(HIGH,KAXIS)
  allocate(kCoords(ksizeGC))
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,iCoords,isizeGC)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,jCoords,jsizeGC)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,kCoords,ksizeGC)

  call Grid_getMinCellSize(min_dx)

  call Flame_getWidth(flamewidth)
  
  ! loop over all zones and init
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        ! check if we are near the flame interface
        ! if we are, then we need to calculate a local piece of HSE for perturbed
        ! flame surface
        ! position of flame surface at this j,k
        fi_buf_size = sim_spert_ampl1+sim_spert_ampl2+flamewidth
        if (NDIM == 1) then
           x_f = sim_x0
        else if (NDIM == 2) then
           x_f = sim_x0 + sim_spert_ampl1 * cos((jCoords(j)/sim_spert_wl1+sim_spert_phase1)*2.0*PI) &
                        + sim_spert_ampl2 * cos((jCoords(j)/sim_spert_wl2+sim_spert_phase2)*2.0*PI)
        else if (NDIM == 3) then
           x_f = sim_x0 + sim_spert_ampl1 * cos((jCoords(j)/sim_spert_wl1+sim_spert_phase1)*2.0*PI) &
                                          * cos((kCoords(k)/sim_spert_wl1+sim_spert_phase1)*2.0*PI) &
                        + sim_spert_ampl2 * cos((jCoords(j)/sim_spert_wl2+sim_spert_phase2)*2.0*PI) &
                                          * cos((kCoords(k)/sim_spert_wl2+sim_spert_phase2)*2.0*PI)
        endif
        fi_x_start = x_f - 2.0*fi_buf_size
        fi_x_end   = x_f + 2.0*fi_buf_size
        fi_start = int(fi_x_start / min_dx)
        fi_end = int(fi_x_end / min_dx) + 1

        if ( cornerID(IAXIS) > fi_end .or. (cornerID(IAXIS)+NXB*stride(IAXIS)) < fi_start ) then
           ! no cells in this block are in flame interface region
           fi_size = 0
        else
           fi_size = fi_end - fi_start + 1
           allocate(fi_flam(fi_size))
           allocate(fi_dens(fi_size))
           allocate(fi_temp(fi_size))
           allocate(fi_ye(fi_size))
           allocate(fi_sumy(fi_size))

           call flame_hse(fi_flam,fi_dens,fi_temp,fi_ye,fi_sumy, &
                          sim_dens_u,sim_temp_u,sim_ye_u,sim_sumy_u, &
                          sim_dens_b,sim_temp_b,sim_ye_b,sim_sumy_b, &
                          x_f, (fi_start-0.5)*min_dx, min_dx, sim_grav, fi_size)
        endif

        ! now initialize data
        ! away from interface, use pre-computed HSE
        ! near interface (within buf_size) use interface region
        ! in crossover region (further than buf_size, less that twice) use
        !    linear crossover
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           dens=0.0
           temp=0.0
           ye  =0.0
           sumy=0.0
           flam=0.0
           do n = 1,stride(IAXIS)
              ! global iaxis index for this (sub) point
              gi = cornerID(IAXIS)+stride(IAXIS)*(i-NGUARD-1)+n-1
              ! x position of this (sub-)point
              xp = (gi-0.5)*min_dx
              if ( xp < fi_x_start .or. xp > fi_x_end) then
                 ! away from flame interface (most of domain)
                 ! just use the overal HSE
                 dens = dens+sim_dens_i(gi)
                 temp = temp+sim_temp_i(gi)
                 ye   = ye+  sim_ye_i(gi)
                 sumy = sumy+sim_sumy_i(gi)
                 flam = flam+sim_flam_i(gi)
              else if (xp <= fi_x_start + fi_buf_size ) then
                 ! lower crossover region
                 w = real(xp-fi_x_start)/fi_buf_size
                 dens = dens + w*fi_dens(gi-fi_start+1) + (1.0-w)*sim_dens_i(gi)
                 temp = temp + w*fi_temp(gi-fi_start+1) + (1.0-w)*sim_temp_i(gi)
                 ye   = ye   + w*fi_ye(gi-fi_start+1) + (1.0-w)*sim_ye_i(gi)
                 sumy = sumy + w*fi_sumy(gi-fi_start+1) + (1.0-w)*sim_sumy_i(gi)
                 flam = flam + w*fi_flam(gi-fi_start+1) + (1.0-w)*sim_flam_i(gi)
              else if (xp >= fi_x_end - fi_buf_size) then
                 ! upper crossover region
                 w = real(fi_x_end-xp)/fi_buf_size
                 dens = dens + w*fi_dens(gi-fi_start+1) + (1.0-w)*sim_dens_i(gi)
                 temp = temp + w*fi_temp(gi-fi_start+1) + (1.0-w)*sim_temp_i(gi)
                 ye   = ye   + w*fi_ye(gi-fi_start+1) + (1.0-w)*sim_ye_i(gi)
                 sumy = sumy + w*fi_sumy(gi-fi_start+1) + (1.0-w)*sim_sumy_i(gi)
                 flam = flam + w*fi_flam(gi-fi_start+1) + (1.0-w)*sim_flam_i(gi)
              else
                 ! flame region
                 dens = dens + fi_dens(gi-fi_start+1)
                 temp = temp + fi_temp(gi-fi_start+1)
                 ye   = ye   + fi_ye(gi-fi_start+1)
                 sumy = sumy + fi_sumy(gi-fi_start+1)
                 flam = flam + fi_flam(gi-fi_start+1)
              endif
           enddo
           dens = dens/stride(IAXIS)
           temp = temp/stride(IAXIS)
           ye = ye/stride(IAXIS)
           sumy = sumy/stride(IAXIS)
           flam = flam/stride(IAXIS)

           eosData(EOS_DENS) = dens
           eosData(EOS_TEMP) = temp
           eosData(EOS_ABAR) = 1.0/sumy
           eosData(EOS_ZBAR) = ye*eosData(EOS_ABAR)
           call Eos(MODE_DENS_TEMP,1,eosData)

           ! compute velocity field
           velx=0.0
           vely=0.0
           if (NDIM >= 2) then
              if (iCoords(i) <= sim_x0 .and. iCoords(i+1) > sim_x0) then
                 velx = sim_vel_pert_amp*cos(2*PI*jCoords(j)/sim_vel_pert_wavelength1)
! was here for testing TFI, turn off for now
!              else if ( iCoords(i) > sim_x0+fi_buf_size ) then
!                 xv = ( iCoords(i) - ( sim_x0+fi_buf_size+7.5e5 ) ) / 7.5e5
!                 yv = jCoords(j) / 7.5e5
!                 vely = xv * 5.0e5 / 2.0 / PI * &
!                        exp((1.0 - xv**2 - yv**2)/2.0)
!                 velx = -yv * 5.0e5 / 2.0 / PI * &
!                        exp((1.0 - xv**2 - yv**2)/2.0)
              endif
           endif

           cell(IAXIS) = i
           cell(JAXIS) = j
           cell(KAXIS) = k
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, cell, dens)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, cell, temp)
           call Grid_putPointData(blockId, CENTER, FLAM_MSCALAR, EXTERIOR, cell, flam)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, cell, velx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, cell, vely)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, cell, 0.0)

           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, cell, eosData(EOS_EINT)+0.5*(velx**2+vely**2))
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, cell, eosData(EOS_EINT))
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, cell, eosData(EOS_PRES))
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, cell, eosData(EOS_GAMC))
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, cell, &
                                       eosData(EOS_PRES)/(dens*eosData(EOS_EINT))+1.0)
        enddo

        ! deallocate space used for this j,k
        if (fi_size > 0 ) then
           deallocate(fi_flam)
           deallocate(fi_dens)
           deallocate(fi_temp)
           deallocate(fi_ye)
           deallocate(fi_sumy)
        endif

     enddo
  enddo
  
  deallocate(iCoords)
  deallocate(jCoords)
  deallocate(kCoords)

  return
  
end subroutine Simulation_initBlock
