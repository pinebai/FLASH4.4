!!****if* source/physics/Hydro/HydroMain/split/PPM/chomboCompatible/hy_ppm_flux_conserve
!!
!! NAME
!!
!!  hy_ppm_flux_conserve
!!
!! SYNOPSIS
!!
!!  call hy_ppm_flux_conserve(integer(in) :: operation,
!!                            integer(in) :: xyzswp,
!!                            integer(in) :: blkcount,
!!                            integer, dimension(blkCount)(in) :: blklist)
!!
!! DESCRIPTION
!!
!! Performs either a PRE or POST flux conservation operation on solution data.
!!
!! ARGUMENTS
!!
!!   operation : PRE or POST Flux conservation operation
!!
!!   xyzswp : xyz swap 
!!
!!   blkcount : block count 
!!
!!   blklist : block list
!!
!!
!!
!!***

!!REORDER(4): solnData

#include "constants.h"
#include "Flash.h"
#include "PPM.h"
#include "Eos.h"
#include "hy_ppm_chombo_constants.h"

#ifdef DEBUG_ALL 
#define DEBUG_HYDRO
#endif
#define DEBUG_GRID_GCMASK


subroutine hy_ppm_flux_conserve(operation, xyzswp, blkCount, blkList)
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
       Grid_releaseBlkPtr
  use Driver_interface, ONLY : Driver_abortFlash
  use Hydro_data, ONLY : hy_smlrho
  implicit none
  integer, intent(in) :: operation, xyzswp, blkCount
  integer, dimension(blkCount), intent(in) :: blkList

  real, pointer :: solnData(:,:,:,:)
  real :: einternal, ekin, etot
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer ::lb, blk, i, j, k


  if (operation == POST_FLUX_CONSERVE) then
     do lb = 1, blkCount
        blk = blkList(lb)
        call Grid_getBlkPtr(blk,solnData)
        call Grid_getBlkIndexLimits(blk,blkLimits,blkLimitsGC)

        do k = blkLimits(1,KAXIS), blkLimits(2,KAXIS)
           do j = blkLimits(1,JAXIS), blkLimits(2,JAXIS)
              do i = blkLimits(1,IAXIS), blkLimits(2,IAXIS)

                 !Discard the updated EINT calculation.
                 etot = solnData(ENER_VAR,i,j,k)

                 !Compute the new kinetic energy.
                 ekin = 0.5e0 * (solnData(VELX_VAR,i,j,k)**2 &
                      + solnData(VELY_VAR,i,j,k)**2 &
                      + solnData(VELZ_VAR,i,j,k)**2)

                 !For now we always update internal energy from etot
                 einternal = etot - ekin
                 solnData(EINT_VAR,i,j,k) = einternal
              end do
           end do
        end do
        call Grid_releaseBlkPtr(blk,solnData)
     end do
  end if
end subroutine hy_ppm_flux_conserve
