!!****if* source/Simulation/SimulationMain/PreHeatPlug/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!! 
!!  2014/11/10  Add the material for filled gas
!!  2015/4/30   Add the material for LEH window and Washer (named as Wash)
!!              WashThickness > TargetThickness and WashRadius > TargetRadius are needed
!!
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putPointData
  use Driver_interface, ONLY: Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT

  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer, intent(in) :: blockId
  
  integer :: i, j, k, n
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  real :: tradActual
  real :: rho, tele, trad, tion, zbar, abar
  integer :: species

#ifndef CHAM_SPEC
  integer :: CHAM_SPEC = 1, TARG_SPEC = 2
  integer :: GAS_SPEC = 3                   !2014/11/10, by Po-Yu
  integer :: LEH_SPEC = 4, WASH_SPEC=5    !2015/4/30, by Po-Yu
#endif


  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., &
       xcent, blkLimitsGC(HIGH, IAXIS))
  allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., &
       ycent, blkLimitsGC(HIGH, JAXIS))
  allocate(zcent(blkLimitsGC(HIGH, KAXIS)))
  call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., &
       zcent, blkLimitsGC(HIGH, KAXIS))

  !------------------------------------------------------------------------------

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           species = CHAM_SPEC
           if (sim_initGeom == "slab") then
              if(NDIM == 1) then
                 if ( xcent(i) <= sim_targetHeight + sim_vacuumHeight .and. &
                      xcent(i) >= sim_vacuumHeight ) then
                    species = TARG_SPEC
                 end if                 
              elseif(NDIM == 2 .or. NDIM == 3) then
                 !if ( xcent(i) <= sim_targetRadius .and. &
                 !     ycent(j) <= sim_targetHeight + sim_vacuumHeight .and. &
                 !     ycent(j) >= sim_vacuumHeight ) then
                 !   species = TARG_SPEC
                 !end if
                 !see p151 @ notebook#17 to get the definition of different zones
                 if ( ycent(j) >= sim_vacuumHeight-sim_windowsThickness .and. ycent(j) < sim_vacuumHeight .and. &
                      !xcent(i) <= sim_targetRadius ) then
                     !species = TARG_SPEC
                      xcent(i) <= sim_washRadius ) then   !2015/4/30, by Po-Yu
                     species = LEH_SPEC !zone 1
                 elseif ( ycent(j) >= sim_vacuumHeight .and. ycent(j) < sim_vacuumHeight+sim_targetThickness ) then
                    if ( xcent(i) <= sim_windowsRadius ) then
                        species = GAS_SPEC !zone 2
                        !species = CHAM_SPEC
                    elseif( xcent(i) <=sim_targetRadius ) then
                        species = TARG_SPEC !zone 3
                    elseif( xcent(i) <=sim_washRadius ) then  !2015/4/30, by Po-Yu
                        !species = WASH_SPEC !zone 4    use this line only when different material is used for washer                
                        species = TARG_SPEC !zone 4                  
                    endif
                 elseif ( ycent(j) >= sim_vacuumHeight .and. ycent(j) < sim_vacuumHeight+sim_washThickness ) then !2015/4/30, by Po-Yu
                     if ( xcent(i) <= sim_targetRadius-sim_targetThickness ) then
                         species = GAS_SPEC !zone 5
                         !species = CHAM_SPEC
                     elseif( xcent(i) <=sim_targetRadius ) then
                         species = TARG_SPEC    !zone 6
                     elseif (xcent(i) <= sim_washRadius ) then
                         !species = WASH_SPEC  !zone 7  use this line only when different material is used for washer                
                         species = TARG_SPEC  !zone 7
                     endif                                          
                 !elseif ( ycent(j) >= sim_vacuumHeight+sim_targetThickness .and. & 
                 elseif ( ycent(j) >= sim_vacuumHeight+sim_washThickness .and. &  !2015/4/30, by Po-Yu
                          !ycent(j) < sim_vacuumHeight+sim_targetHeight-sim_targetThickness ) then
                          ycent(j) < sim_vacuumHeight+sim_targetHeight-sim_plugThickness ) then !2015/9/12, by Po-Yu
                     if ( xcent(i) <= sim_targetRadius-sim_targetThickness ) then
                         species = GAS_SPEC !zone 8
                         !species = CHAM_SPEC
                     elseif( xcent(i) <=sim_targetRadius ) then
                         species = TARG_SPEC    !zone 9
                     endif
                 !elseif ( ycent(j) >= sim_vacuumHeight+sim_targetHeight-sim_targetThickness .and. &
                 elseif ( ycent(j) >= sim_vacuumHeight+sim_targetHeight-sim_plugThickness .and. &   !2015/9/12, by Po-Yu
                          ycent(j) < sim_vacuumHeight+sim_targetHeight .and. &
                         xcent(i) <= sim_targetRadius ) then
                     species = TARG_SPEC    !zone 10
                 end if
              end if
           else
               if (sqrt(xcent(i)**2+ycent(j)**2+zcent(k)**2)<= sim_targetRadius) then
                  species = TARG_SPEC
              end if
           end if
           if(species == TARG_SPEC) then
              rho = sim_rhoTarg
              tele = sim_teleTarg
              tion = sim_tionTarg
              trad = sim_tradTarg
           elseif(species == GAS_SPEC) then
              rho = sim_rhoGas
              tele = sim_teleGas
              tion = sim_tionGas
              trad = sim_tradGas               
           elseif(species == LEH_SPEC) then !2015/4/30, by Po-Yu
              rho = sim_rhoLEH
              tele = sim_teleLEH
              tion = sim_tionLEH
              trad = sim_tradLEH                             
           elseif(species == WASH_SPEC) then !2015/4/30, by Po-Yu
              rho = sim_rhoWash
              tele = sim_teleWash
              tion = sim_tionWash
              trad = sim_tradWash              
           else
              rho = sim_rhoCham
              tele = sim_teleCham
              tion = sim_tionCham
              trad = sim_tradCham
           end if

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, tele)

#ifdef FLASH_3T
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, tion)
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, tele)

           ! Set up radiation energy density:
           call RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, tradActual)
#endif
           if (NSPECIES > 0) then
              ! Fill mass fractions in solution array if we have any SPECIES defined.
              ! We put nearly all the mass into either the Xe material if XE_SPEC is defined,
              ! or else into the first species.
              do n = SPECIES_BEGIN,SPECIES_END
                 if (n==species) then
                    call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, 1.0e0-(NSPECIES-1)*sim_smallX)
                 else
                    call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, sim_smallX)
                 end if
              enddo
           end if

        enddo
     enddo
  enddo

  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)

  return

end subroutine Simulation_initBlock
