!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_accExtrnPotRow
!!
!! NAME
!!
!!  grv_accExtrnPotRow
!!
!!
!! SYNOPSIS
!!
!!  call grv_accExtrnPotRow
!!           integer(in) :: pos(2),
!!           integer(in) :: sweepDir,
!!           integer(in) :: blockID,
!!           integer(in) :: numCells,
!!           real(inout) :: grav(numCells)
!!  )
!!
!! DESCRIPTION
!!
!!   Interpolates in the external gravitational field and returns
!!   one component gravitational acceleration for a row of grid cells.
!!
!! ARGUMENTS
!!
!!   pos      - Row indices transverse to the sweep direction
!!   sweepDir - The sweep direction:  allowed values are 
!!              SWEEP_X, SWEEP_Y, and SWEEP_Z. These values are defined
!!              in constants.h.
!!   blockID  - The local identifier of the block
!!   numCells - Number of cells to update in grav array
!!   grav     - Array to receive result
!!
!! RESULT
!!
!!   Vector grav containing a single component of the gravitational 
!!   acceleration given by the external field.
!!
!!***

subroutine grv_accExtrnPotRow(pos, sweepDir, blockID, numCells, grav)
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords
  use Driver_interface, ONLY : Driver_abortFlash
  use Gravity_data, ONLY : grv_bhExtrnPotN, &
    grv_bhExtrnPotCoord, grv_bhExtrnPotPot, grv_bhExtrnPotAcc, &
    grv_bhExtrnPotDel, grv_useExternalPotential, grv_bhExtrnPotIType, &
    grv_bhEPTypeR, grv_bhEPTypeX, grv_bhEPTypeY, grv_bhEPTypeZ, &
    grv_bhExtrnPotCenterX, grv_bhExtrnPotCenterY, grv_bhExtrnPotCenterZ

  implicit none
#include "constants.h"
#include "Flash.h"
  integer, dimension(2), intent(in) :: pos
  integer, intent(in)               :: sweepDir, blockID,  numCells
  real, intent(inout)               :: grav(numCells)
  integer :: ii, ind
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_KHI_GC) :: zCenter
  real,dimension(GRID_JHI_GC) :: yCenter
  real,dimension(GRID_IHI_GC) :: xCenter
#else
  real,allocatable,dimension(:) ::xCenter,yCenter,zCenter
#endif
  real :: delr, tmpdr2, Fr, p
  integer :: sizeX,sizeY,sizeZ
  logical :: gcell = .true.

  grav(:) = 0.0
  if (.not. grv_useExternalPotential) return

#ifndef FIXEDBLOCKSIZE
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    sizeX=blkLimitsGC(HIGH,IAXIS)
    sizeY=blkLimitsGC(HIGH,JAXIS)
    sizeZ=blkLimitsGC(HIGH,KAXIS)
#else
    sizeX=GRID_IHI_GC
    sizeY=GRID_JHI_GC
    sizeZ=GRID_KHI_GC
#endif

  select case (grv_bhExtrnPotIType)
    case(grv_bhEPTypeR)
#ifndef FIXEDBLOCKSIZE
      allocate(xCenter(sizeX))
      allocate(yCenter(sizeY))
      allocate(zCenter(sizeZ))
#endif
      zCenter = 0.
      yCenter = 0.
      if (NDIM == 3) then 
         call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCenter, sizeZ)
         zCenter = zCenter - grv_bhExtrnPotCenterZ
      endif
      if (NDIM >= 2) then
         call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCenter, sizeY)
         yCenter = yCenter - grv_bhExtrnPotCenterY
      endif
      call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCenter, sizeX)
      xCenter = xCenter - grv_bhExtrnPotCenterX

      if (sweepDir .eq. SWEEP_X) then                       ! x-component
         tmpdr2 = yCenter(pos(1))*yCenter(pos(1)) + zCenter(pos(2))*zCenter(pos(2)) 
         do ii = 1, numCells
            delr = sqrt(xCenter(ii)*xCenter(ii) + tmpdr2)
            ind = floor(delr/grv_bhExtrnPotDel) + 1
            p = (delr - grv_bhExtrnPotCoord(ind)) / (grv_bhExtrnPotCoord(ind+1) - grv_bhExtrnPotCoord(ind))
            Fr = (1.-p)*grv_bhExtrnPotAcc(ind) + p*grv_bhExtrnPotAcc(ind+1)
            grav(ii) = Fr*xCenter(ii) / delr
         enddo
     
      else if (sweepDir .eq. SWEEP_Y) then          ! y-component
         tmpdr2 = xCenter(pos(1))*xCenter(pos(1)) + zCenter(pos(2))*zCenter(pos(2)) 
         do ii = 1, numCells
            delr = sqrt(yCenter(ii)*yCenter(ii) + tmpdr2)
            ind = floor(delr/grv_bhExtrnPotDel) + 1
            p = (delr - grv_bhExtrnPotCoord(ind)) / (grv_bhExtrnPotCoord(ind+1) - grv_bhExtrnPotCoord(ind))
            Fr = (1.-p)*grv_bhExtrnPotAcc(ind) + p*grv_bhExtrnPotAcc(ind+1)
            grav(ii) = Fr*yCenter(ii) / delr
         enddo
     
      else if (sweepDir .eq. SWEEP_Z) then          ! z-component
         tmpdr2 = xCenter(pos(1))*xCenter(pos(1)) + yCenter(pos(2))*yCenter(pos(2)) 
         do ii = 1, numCells
            delr = sqrt(zCenter(ii)*zCenter(ii) + tmpdr2)           
            ind = floor(delr/grv_bhExtrnPotDel) + 1
            p = (delr - grv_bhExtrnPotCoord(ind)) / (grv_bhExtrnPotCoord(ind+1) - grv_bhExtrnPotCoord(ind))
            Fr = (1.-p)*grv_bhExtrnPotAcc(ind) + p*grv_bhExtrnPotAcc(ind+1)
            grav(ii) = Fr*zCenter(ii) / delr
         enddo
      endif
#ifndef FIXEDBLOCKSIZE
      deallocate(xCenter)
      deallocate(yCenter)
      deallocate(zCenter)
#endif
    case(grv_bhEPTypeX)
      call Driver_abortFlash("grv_readExtrnPotential: planex symmetry not supported yet")
    case(grv_bhEPTypeY)
      call Driver_abortFlash("grv_readExtrnPotential: planey symmetry not supported yet")
    case(grv_bhEPTypeZ)
     !! Get local ccordinates first
#ifndef FIXEDBLOCKSIZE
     allocate(zCenter(sizeZ),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate kCoord in ExternalPotential.F90")
#endif
     call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCenter, sizeZ)
     zCenter = zCenter - grv_bhExtrnPotCenterZ
     
     !! do loop over z-coord
     if (sweepDir == SWEEP_Z) then
       do ii = 1, numCells
         ind = floor(abs(zCenter(ii))/grv_bhExtrnPotDel) + 1
         if (ind > grv_bhExtrnPotN) &
           & call Driver_abortFlash("grv_accExtrnPotRow: out of external pot. array")
         p = (abs(zCenter(ii)) - grv_bhExtrnPotCoord(ind)) &
         & / (grv_bhExtrnPotCoord(ind+1) - grv_bhExtrnPotCoord(ind))
         if (zCenter(ii) < 0.0) then
           grav(ii) =  (1.-p)*grv_bhExtrnPotAcc(ind) + p*grv_bhExtrnPotAcc(ind+1)
         else
           grav(ii) = -(1.-p)*grv_bhExtrnPotAcc(ind) - p*grv_bhExtrnPotAcc(ind+1)
         endif
       enddo
     else
       grav(:) = 0.0
     endif
     
#ifndef FIXEDBLOCKSIZE
     deallocate(zCenter)
#endif

    case default
      call Driver_abortFlash("grv_readExtrnPotential: unrecognized potential type")
  end select

end subroutine grv_accExtrnPotRow


