!!****if* source/physics/Hydro/HydroMain/split/PPM/chomboCompatible/hy_ppm_chomboLikeUpdateSoln
!! NAME
!!
!!  hy_ppm_chomboLikeUpdateSoln
!!
!!
!! SYNOPSIS
!!
!!  hy_ppm_chomboLikeUpdateSoln(integer, intent(IN)  :: rangeSwitch, 
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
!!REORDER(4): solnData, tempFlx

#include "constants.h"
#include "PPM.h"
#include "Flash.h"

!Subroutine updates solution data using conserved fluxes in a Chombo like way.
subroutine hy_ppm_chomboLikeUpdateSoln(rangeSwitch,                        &
                         xyzswp, dt,                          &
                         blkLimits,blkLimitsGC,numCells,               &
                         tempArea, tempGrav1d_o, tempGrav1d,  &
                         tempDtDx, tempFict,                  &
                         tempFlx,  solnData)

  use Hydro_data, ONLY: hy_smlrho, hy_smallp, hy_eintSwitch
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
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
  real :: flx1,flx2  !Left and right cell fluxes
  real :: dtdx

  !One-to-one mapping between Fdef and Udef.
  integer, dimension(NPROP_FLUX) :: Fdef, Udef 
  integer, dimension(3,MDIM) :: loopIdx !First dimension: LOW,HIGH,SKIP
  integer :: i, j, k, v, fIdx, uIdx


  if (rangeSwitch /= UPDATE_BOUND) then
     call Driver_abortFlash("Only designed to update internal boundary cells")
  end if


  !By default set the lower and upper cell indices and skip to 1
  loopIdx(1:3,1:MDIM) = 1

  !Set the lower and upper cell indices
  loopIdx(1:2,1:NDIM) = blkLimits(1:2,1:NDIM)

  !Set the fixed mapping between named flux and solution variables.
  Fdef(1) = E_FLUX;    Udef(1) = ENER_VAR
  Fdef(2) = EINT_FLUX; Udef(2) = EINT_VAR
  Fdef(3) = P_FLUX;    Udef(3) = PRES_VAR
  Fdef(4) = RHO_FLUX;  Udef(4) = DENS_VAR !This must be at densIdx index

  !Set the cell skip and the mapping between flux and solution velocities
  Fdef(5) = U_FLUX
  Fdef(6) = UT_FLUX
  Fdef(7) = UTT_FLUX

  if (xyzswp == SWEEP_X) then
     loopIdx(3,IAXIS) = blkLimits(2,IAXIS) - blkLimits(1,IAXIS)
     Udef(5) = VELX_VAR
     Udef(6) = VELY_VAR
     Udef(7) = VELZ_VAR
  else if (xyzswp == SWEEP_Y) then
     loopIdx(3,JAXIS) = blkLimits(2,JAXIS) - blkLimits(1,JAXIS)
     Udef(5) = VELY_VAR
     Udef(6) = VELX_VAR
     Udef(7) = VELZ_VAR
  else if (xyzswp == SWEEP_Z) then
     loopIdx(3,KAXIS) = blkLimits(2,KAXIS) - blkLimits(1,KAXIS)
     Udef(5) = VELZ_VAR
     Udef(6) = VELX_VAR
     Udef(7) = VELY_VAR
  else
     call Driver_abortFlash("Unknown sweep axis")
  end if


  !Loop over the appropriate boundary cells for each sweep axis.
  do k = loopIdx(1,KAXIS), loopIdx(2,KAXIS), loopIdx(3,KAXIS)
     do j = loopIdx(1,JAXIS), loopIdx(2,JAXIS), loopIdx(3,JAXIS)
        do i = loopIdx(1,IAXIS), loopIdx(2,IAXIS), loopIdx(3,IAXIS)

           dtdx  = tempDtDx(i,j,k)

           !Retrieve the appropriate fluxes and then update the
           !solution according to the left and right cell fluxes
           do v = 1, NFLUXES
              if (v <= NPROP_FLUX) then
                 fIdx = Fdef(v)
                 uIdx = Udef(v)
              else
                 fIdx = v
                 uIdx = v - NPROP_FLUX + NPROP_VARS
              end if


              if (uIdx /= PRES_VAR) then
                 if (xyzswp == SWEEP_X) then
                    flx1 = tempFlx(fIdx,i,j,k)
                    flx2 = tempFlx(fIdx,i+1,j,k)
                 else if (xyzswp == SWEEP_Y) then
                    flx1 = tempFlx(fIdx,i,j,k)
                    flx2 = tempFlx(fIdx,i,j+1,k)
                 else if (xyzswp == SWEEP_Z) then
                    flx1 = tempFlx(fIdx,i,j,k)
                    flx2 = tempFlx(fIdx,i,j,k+1)
                 else
                    call Driver_abortFlash("Unknown sweep axis")
                 end if

                 solnData(uIdx,i,j,k) = solnData(uIdx,i,j,k) - &
                      dtdx * (flx2 - flx1)
              end if
           end do

        end do
     end do
  end do

end subroutine hy_ppm_chomboLikeUpdateSoln
