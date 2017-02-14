!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree/gr_markJeans
!!
!! NAME
!!  gr_markJeans
!!
!! SYNOPSIS
!!  call gr_markJeans()
!!
!! PURPOSE
!!  Refine all blocks so that the size of a grid cell will be always smaller than
!!  the Jeans length.
!!
!! ARGUMENTS
!!
!! NOTES
!! 
!!  This routine is specific to the Barnes-Hut-Tree test.
!!
!!***

subroutine gr_markJeans()

  use tree, ONLY : refine, derefine, stay, lrefine, lnblocks, nodetype
  use Grid_data, ONLY : gr_geometry
  use Grid_interface
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Multispecies_interface, ONLY : Multispecies_getSumInv  
  use Simulation_data, ONLY : jeans_ref, jeans_deref
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"  

  real :: delx2, dely2, delz2, del2
  real :: lambda_J2, lJc, newton, gam
  integer :: i,j,k, sizeX, sizeY, sizeZ
  real, dimension(NSPECIES) :: massFrac  
  integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  integer :: blockID, blockCount
  integer :: blockList(MAXBLOCKS)
  logical       :: gcell = .true.
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC) :: xCoord
  real, dimension(GRID_JHI_GC) :: yCoord
  real, dimension(GRID_KHI_GC) :: zCoord
#else
  real,allocatable, dimension(:) :: xCoord,yCoord,zCoord
  integer :: istat
#endif
!-------------------------------------------------------------------------------

  !print *, "entering gr_markJeans"
  call PhysicalConstants_get ("Newton", newton)

! Loop over all leaf-node blocks.
  lJc = 15./4./PI/newton

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  do blockID = 1, blockCount

     ! Get cell coordinates for this block
     call Grid_getBlkIndexLimits(blockList(blockID),blkLimits,blkLimitsGC)
     sizeX=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     sizeY=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     sizeZ=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
     allocate(xCoord(sizeX),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate xCoord in gr_markJeans.F90")
     allocate(yCoord(sizeY),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate yCoord in gr_markJeans.F90")
     allocate(zCoord(sizeZ),stat=istat)
     if (istat .ne. 0) call Driver_abortFlash("could not allocate zCoord in gr_markJeans.F90")
#endif
     ! x coordinates
     call Grid_getCellCoords(IAXIS,blockList(blockID),&
          CENTER,gcell,xCoord,sizeX)
#if NDIM > 1
     ! y coordinates
     call Grid_getCellCoords(JAXIS,blockList(blockID),&
          CENTER,gcell,yCoord,sizeY)
#if NDIM > 2
     ! z coordinates
     call Grid_getCellCoords(KAXIS,blockList(blockID),&
          CENTER,gcell,zCoord,sizeZ)
#endif
#endif
  
    call Grid_getBlkPtr(blockList(blockID),solnData,CENTER)


    !refine(blockList(blockID)) = .false.
    !derefine(blockList(blockID)) = .true.
    
    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
               
           
          delx2 = 0.25*(xCoord(i+1) - xCoord(i-1))**2

          dely2 = 0.25*(yCoord(j+1) - yCoord(j-1))**2
          if ((gr_geometry == SPHERICAL) .or. (gr_geometry == POLAR)) then
            dely2 = dely2 * xCoord(i)**2
          end if

          delz2 = 0.25*(zCoord(k+1) - zCoord(k-1))**2
          if (gr_geometry == SPHERICAL) then
            delz2 = delz2 * (xCoord(i) * sin(yCoord(j)))**2
          else if (gr_geometry == CYLINDRICAL) then
            delz2 = delz2 * xCoord(i)**2
          end if

          massFrac(1) = solnData(SPECIES_BEGIN,i,j,k)
          massFrac(2) = solnData(SPECIES_BEGIN+1,i,j,k)

          call Multispecies_getSumInv(GAMMA, gam, massFrac(:))
          gam       = 1/gam

          lambda_J2 = (gam-1.0)*lJc*solnData(EINT_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)
          del2 = max(delx2, max(dely2, delz2))

          
          !print *, "JEANS: ", gam,lambda_J2, sqrt(lambda_J2), sqrt(del2), blockList(blockID), i,j,k
          if (lambda_J2 < del2*jeans_ref**2 ) then
            ! Jeans criterion needs refinement
            refine(blockList(blockID)) = .true.
            derefine(blockList(blockID)) = .false.
          else if (lambda_J2 < del2*jeans_deref**2 .and. (.not. refine(blockList(blockID)))) then
            ! Jeans criterion suppress derefinement
            stay(blockList(blockID)) = .true.
            derefine(blockList(blockID)) = .false.
          else if (lambda_J2 > del2*jeans_deref**2 .and. (.not. refine(blockList(blockID))) .and. &
                                                         (.not. stay(blockList(blockID)))) then
            ! Jeans criterion allows derefinement
            derefine(blockList(blockID)) = .true.
          endif

          
        end do
      end do
    end do
    !print *, "refine to True: ", sqrt(lambda_J2), sqrt(del2), blockList(blockID), i,j,k

    !print *, "JEANS: ", blockList(blockID), refine(blockList(blockID)), derefine(blockList(blockID))

#ifndef FIXEDBLOCKSIZE
     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)
#endif

  enddo

!-------------------------------------------------------------------------------

  return
end subroutine gr_markJeans
