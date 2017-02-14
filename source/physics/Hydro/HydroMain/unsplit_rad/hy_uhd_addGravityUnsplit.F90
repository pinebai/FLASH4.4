!!****if* source/physics/Hydro/HydroMain/unsplit_rad/hy_uhd_addGravityUnsplit
!!
!!
!! NAME
!!
!!  hy_uhd_addGravityUnsplit
!!
!! SYNOPSIS
!!
!!  call hy_uhd_addGravityUnsplit( integer(IN) :: blockID,
!!                            integer(IN) :: blkLimits(LOW:HIGH,MDIM),
!!                            integer(IN) :: dataSize(MDIM),
!!                            real(IN)    :: dt,
!!                            real(IN)    :: gravX,
!!                            real(IN)    :: gravY,
!!                            real(IN)    :: gravZ )
!!
!!
!! DESCRIPTION
!!
!!  Adds the second part of the gravitational force to the momenta and energy.
!!  The first half is added in unsplitUpdate.  This centers the gravitational
!!  force appropriately on the current timestep rather than extrapolating.
!!
!! ARGUMENTS
!!
!!  blockID   - local block ID
!!  blkLimits - block limits 
!!  dataSize  - array size of gravX,Y,Z in non-fixed block size mode in UG
!!  dt        - timestep
!!  gravX     - gravity source term in x-direction
!!  gravY     - gravity source term in y-direction
!!  gravZ     - gravity source term in z-direction
!!
!! NOTES
!!
!!  The ENER_VAR component of UNK is now always in specific
!!  (i.e., energy per mass) form, on entry as well as on return from this routine.
!!  This is changed from the behavior in FLASH4.2.2 and earlier.
!!
!!***

!!REORDER(4):U

Subroutine hy_uhd_addGravityUnsplit&
     (blockID,blkLimits,dataSize,dt,gravX,gravY,gravZ)

  use Hydro_data,      ONLY : hy_useGravity,        &
                              hy_useGravHalfUpdate

  use Grid_interface,  ONLY : Grid_getBlkPtr,       &
                              Grid_releaseBlkPtr

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! ---- Argument List ----------------------------------
  integer, intent(IN) :: blockID
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
  integer, dimension(MDIM), intent(IN) :: dataSize
  real,    intent(IN) :: dt

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: &
       gravX,gravY,gravZ
#else
  real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(IN) :: &
       gravX,gravY,gravZ
#endif
  !! -----------------------------------------------------
  real :: hdt
  integer :: i,j,k
  real, pointer, dimension(:,:,:,:) :: U
  real, dimension(3) :: velNew

  hdt = 0.5 * dt

  !! Get block pointer for storages of Riemann states
  call Grid_getBlkPtr(blockID,U,CENTER)

  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           velNew(1:3) = U(VELX_VAR:VELZ_VAR,i,j,k)&
                + hdt*(/gravX(i,j,k),gravY(i,j,k),gravZ(i,j,k)/)

           U(ENER_VAR,i,j,k) = U(ENER_VAR,i,j,k) &
                + hdt*dot_product(velNew(1:3),(/gravX(i,j,k),gravY(i,j,k),gravZ(i,j,k)/))

           U(VELX_VAR:VELZ_VAR,i,j,k) = velNew(1:3)
        enddo
     enddo
  enddo

  !! Release block pointer for storages of Riemann states
  call Grid_releaseBlkPtr(blockID,U,CENTER)

end Subroutine hy_uhd_addGravityUnsplit
