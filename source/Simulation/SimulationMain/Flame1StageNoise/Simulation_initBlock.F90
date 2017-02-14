!!****if* source/Simulation/SimulationMain/Flame1StageNoise/Simulation_initBlock
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

!!  Initialization of flame in pressure equilibrium
!!  The flame front can be planar or spherical in multiple dimensions
!!  depending on the value of the pseudo_1d parameter.
!!  See description of parameters in Config file for more info.
!!
! Dean Townsley 2008

subroutine Simulation_initBlock(blockID)
  
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkBoundBox, Grid_getDeltas, Grid_putPointData, &
    Grid_getCellCoords
  use Flame_interface, ONLY : Flame_getProfile, Flame_rhJump
  use fl_effData, ONLY: fl_effDeltae

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  
  integer, intent(in) :: blockID

  integer :: i, j, k

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: cell
  integer :: isizeGC, jsizeGC, ksizeGC
  real, allocatable, dimension(:) :: iCoords, jCoords, kCoords

  real, dimension(EOS_NUM) :: eosData
  real :: flam, velx
  real :: fsurf_x_position, fsurf_distance
  real :: ye, yi

!==============================================================================

  ! get essential info about this block - index limits and cell coordinates
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  isizeGC = blkLimitsGC(HIGH,IAXIS)
  allocate(iCoords(isizeGC))
  jsizeGC = blkLimitsGC(HIGH,JAXIS)
  allocate(jCoords(jsizeGC))
  ksizeGC = blkLimitsGC(HIGH,KAXIS)
  allocate(kCoords(ksizeGC))
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,iCoords,isizeGC)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,jCoords,jsizeGC)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,kCoords,ksizeGC)

  !-----------------------------------------------
  ! loop over all zones and init
  !-----------------------------------------------
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           if (.not. sim_ignite) then
              ! no burned material, only unburned
              eosData(:) = sim_eosData_u(:)
              flam = 0.0
              velx=0.0
           else
              !-----------------------------------------------
              ! initialize, including a burned region
              !-----------------------------------------------

              ! find distance from flame surface (positive is in front of flame)
              if (sim_pseudo1d) then
                 ! planar flame surface with normal tilted up by angle theta from x direction
                 fsurf_x_position = sim_fracPerturb*(sim_xmax-sim_xmin) &
                                    - tan(PI*sim_theta/180.0)*(jCoords(j)-sim_yctrPerturb)
                 fsurf_distance = iCoords(i) - fsurf_x_position
              else
                 ! n-dimensional sphere centered at specified point with radius
                 !   frac_perturb* (x domain extent)
                 fsurf_distance = sqrt( (iCoords(i)-sim_xctrPerturb)**2 + &
                                        (jCoords(j)-sim_yctrPerturb)**2 + &
                                        (kCoords(k)-sim_zctrPerturb)**2 ) &
                                  - sim_fracPerturb*(sim_xmax-sim_xmin)
              endif

              ! determine local state in this zone
              if ( fsurf_distance > 1.5*sim_laminarWidth ) then
                 ! unburned material
                 eosData(:) = sim_eosData_u(:)
                 flam = 0.0
              else if ( fsurf_distance < -1.5*sim_laminarWidth ) then
                 ! fully burned
                 eosData(:) = sim_eosData_b(:)
                 flam = 1.0
              else
                 ! partially burned
                 call Flame_getProfile(fsurf_distance, flam)

                 ! calculate propertise for partially burned material
                 ! note, in fact ye_f and ye_a should be equal
                 yi = 1.0/sim_eosData_u(EOS_ABAR)*(1.0-flam) + (1.0/sim_eosData_b(EOS_ABAR))*flam
                 ye = sim_eosData_u(EOS_ZBAR)/sim_eosData_u(EOS_ABAR)*(1.0-flam) + (sim_eosData_b(EOS_ZBAR)/sim_eosData_b(EOS_ABAR))*flam
                 eosData(:) = sim_eosData_u(:)
                 eosData(EOS_ABAR) = 1.0/yi
                 eosData(EOS_ZBAR) = ye*eosData(EOS_ABAR)

                 ! put this in pressure equilibrium with unburned material
                 call Flame_rhJump(sim_eosData_u, eosData, flam*fl_effDeltae, 0.0, MODE_DENS_TEMP)

              endif

              ! init velocity field, nonzero only makes sense with a planar flame front
              if (sim_pseudo1d) then
                 velx = sim_flamespeed*sim_eosData_u(EOS_DENS)* &
                            (1.e0/sim_eosData_b(EOS_DENS) - 1.e0/eosData(EOS_DENS))
              else
                 velx = 0.0
              endif

           endif ! sim_ignite
           

           !-----------------------------------------------
           !  Now store all this info on the grid
           !-----------------------------------------------
           cell(IAXIS) = i
           cell(JAXIS) = j
           cell(KAXIS) = k
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, cell, eosData(EOS_DENS))
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, cell, eosData(EOS_TEMP))

           call Grid_putPointData(blockId, CENTER, FLAM_MSCALAR, EXTERIOR, cell, flam)

           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, cell, velx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, cell, 0.0)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, cell, 0.0)

           !  usually I would just call the EOS, but we happen to have all this data
           !  so we'll just put it in.
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, cell, eosData(EOS_EINT)+0.5*velx**2)
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, cell, eosData(EOS_EINT))
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, cell, eosData(EOS_PRES))
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, cell, eosData(EOS_GAMC))
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, cell, &
                                       eosData(EOS_PRES)/(eosData(EOS_DENS)*eosData(EOS_EINT))+1.0)
        enddo
     enddo
  enddo
  
  deallocate(iCoords)
  deallocate(jCoords)
  deallocate(kCoords)

  return
  
end subroutine Simulation_initBlock






