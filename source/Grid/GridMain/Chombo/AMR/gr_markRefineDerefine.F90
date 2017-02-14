!!****if* source/Grid/GridMain/Chombo/AMR/gr_markRefineDerefine
!!
!! NAME
!!  gr_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  gr_markRefineDerefine(integer(IN) :: iref,
!!                        real(IN) :: refine_cutoff,
!!                        real(IN) :: derefine_cutoff,
!!                        real(IN) :: refine_filter)
!!  
!!  DESCRIPTION
!!  
!!    Blocks are marked for refining or derefining.
!!    This version uses the second derivative calculations on the specified variable to 
!!    determine if the block needs more resoultion (refine) or less resolution (derefine)
!!    de/refine_cutoff are the thresholds for triggering the corresponding action.
!!    Once the blocks have been marked, the control is passed to Paramesh to update refinement.
!!
!!  ARGUMENTS 
!!
!!    iref - index of the refinement variable in data structure "unk"
!!
!!    refine_cutoff - the threshold value for triggering refinement 
!!
!!    derefine_cutoff - the threshold for triggereing derefinement
!!
!!    refine_filter - makes sure that error calculations to determine refinement
!!                    don't diverge numerically 
!! 
!!  NOTES
!!  
!!    See Grid_markRefineDerefine
!!
!!  SEE ALSO
!!  
!!    Grid_markRefineDerefine
!!
!!***

#include "constants.h"
#include "flash_bool.h"
#include "Flash.h"

!Just in case.  There should be no reference to these defines!
#undef NXB
#undef NYB
#undef NZB

!!REORDER(4): ptrBlk, ptrTagBlk

subroutine gr_markRefineDerefine(&
     iref,refine_cutoff,derefine_cutoff,refine_filter)
  use iso_c_binding, ONLY : c_int
  use Driver_interface, ONLY : Driver_abortFlash
  use flash_ftypes, ONLY : box_info_t
  use chombo_f_c_interface, ONLY : ch_get_box_info
  use Grid_data, ONLY: gr_geometry, gr_domainBC, gr_tagRadius
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getBlkPhysicalSize, Grid_getBlkIndexLimits
  implicit none
  include "Flash_mpi.h"

  integer, intent(IN) :: iref
  real, intent(IN) :: refine_cutoff, derefine_cutoff, refine_filter
  integer, parameter :: SQNDIM = NDIM*NDIM

  real, dimension(MDIM) :: bsize
  integer :: nxb, nyb, nzb, blockID
  real delx,dely,delz
  real dely_f, delz_f

  real, allocatable, dimension(:,:,:,:) :: delu, delua
  integer,dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer ,dimension(MAXBLOCKS) :: blkList
  integer :: blkCount
  real delu2(SQNDIM), delu3(SQNDIM), delu4(SQNDIM)
  real num,denom,error

  integer lb,i,j,k
  integer kstart,kend,jstart,jend,istart,iend
  !
  integer :: ii, jj, kk, ip, jp, kp
  integer :: ndim2

  real, pointer :: solnData(:,:,:) 
  real, pointer :: ptrBlk(:,:,:,:), ptrTagBlk(:,:,:,:)
  type(box_info_t) :: boxInfo
  integer(c_int) :: blkID, gds
  integer :: ir, jr, kr
  gds = CENTER

  if (gr_geometry == SPHERICAL .or. &
       gr_geometry == POLAR .or. &
       gr_geometry == CYLINDRICAL) then
     call Driver_abortFlash("Geometry not yet considered")
  end if


  ndim2 = NDIM*NDIM
  call Grid_getListOfBlocks(ALL_BLKS, blkList, blkCount)

  do lb = 1,blkCount
     blockID = blkList(lb)
     call Grid_getBlkPtr(blockID, ptrBlk, CENTER)
     call Grid_getBlkPtr(blockID, ptrTagBlk, SCRATCH_CTR)
     solnData => ptrBlk(iref,:,:,:)

     call Grid_getBlkPhysicalSize(blockID, bsize)
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     nxb = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
     nyb = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
     nzb = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1

     allocate(delu(MDIM, &
          blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
     allocate(delua(MDIM, &
          blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))


     delx = 0.5e0*float(nxb)/bsize(1)
#if N_DIM >= 2
     dely = 0.5e0*float(nyb)/bsize(2)
     dely_f = dely
#endif
#if N_DIM == 3
     delz = 0.5e0*float(nzb)/bsize(3)
     delz_f = delz
#endif

     ! Compute first derivatives

     do k = 1+K3D,nzb+(K3D*((2*NGUARD)-1))
        do j = 1+K2D,nyb+(K2D*((2*NGUARD)-1))
           do i = 2,nxb+(2*NGUARD)-1

              ! d/dx
              delu(1,i,j,k) = solnData(i+1,j,k) - solnData(i-1,j,k)
              delu(1,i,j,k) = delu(1,i,j,k)*delx

              delua(1,i,j,k) = abs(solnData(i+1,j,k)) + &
                   abs(solnData(i-1,j,k))
              delua(1,i,j,k) = delua(1,i,j,k)*delx

#if N_DIM >= 2
              ! d/dy
              delu(2,i,j,k) = solnData(i,j+1,k) - solnData(i,j-1,k)
              delu(2,i,j,k) = delu(2,i,j,k)*dely_f

              delua(2,i,j,k) = abs(solnData(i,j+1,k)) + &
                   abs(solnData(i,j-1,k))
              delua(2,i,j,k) = delua(2,i,j,k)*dely_f
#endif

#if N_DIM == 3
              ! d/dz
              delu(3,i,j,k) = solnData(i,j,k+1) -  solnData(i,j,k-1)
              delu(3,i,j,k) = delu(3,i,j,k)*delz_f

              delua(3,i,j,k) = abs(solnData(i,j,k+1)) + &
                   abs(solnData(i,j,k-1))
              delua(3,i,j,k) = delua(3,i,j,k)*delz_f
#endif

           end do
        end do
     end do

     ! Compute second derivatives

     ! Test if at a block boundary

     ! Two guardcells
     kstart = 2*K3D+1
     kend   = nzb+(K3D*((2*NGUARD)-2))
     jstart = 2*K2D+1
     jend   = nyb+(K2D*((2*NGUARD)-2))
     istart = 3
     iend   = nxb+(2*NGUARD)-2
     ! One guardcell
     !            kstart = 2*K3D+1+K3D
     !            kend   = nzb+(K3D*((2*NGUARD)-2))-K3D
     !            jstart = 2*K2D+1+K2D
     !            jend   = nyb+(K2D*((2*NGUARD)-2))-K2D
     !            istart = NGUARD
     !            iend   = nxb+(2*NGUARD)-3
     ! No guardcells
     !            kstart = K3D*NGUARD+1
     !            kend   = nzb+K3D*NGUARD
     !            jstart = K2D*NGUARD+1
     !            jend   = nyb+K2D*NGUARD
     !            istart = NGUARD+1
     !            iend   = nxb+NGUARD
     !If the value is -20 or lower then this face represents an external boundary
     !CD - In PARAMESH "neigh" gives the wrap around neighbor for periodic BCs.
     blkID = blockID
     call ch_get_box_info(blkID, gds, boxInfo)

     if ( (boxInfo % isNextToLowDomain(1) == FLASH_TRUE) .and. &
          (gr_domainBC(LOW,1) /= PERIODIC) ) istart = NGUARD+1
     if ( (boxInfo % isNextToHighDomain(1) == FLASH_TRUE) .and. &
          (gr_domainBC(HIGH,1) /= PERIODIC) ) iend   = NGUARD+nxb     
#if N_DIM >= 2
     if ( (boxInfo % isNextToLowDomain(2) == FLASH_TRUE) .and. &
          (gr_domainBC(LOW,2) /= PERIODIC) ) jstart = NGUARD*K2D+1
     if ( (boxInfo % isNextToHighDomain(2) == FLASH_TRUE) .and. &
          (gr_domainBC(HIGH,2) /= PERIODIC) ) jend   = NGUARD*K2D+nyb
#endif
#if N_DIM == 3
     if ( (boxInfo % isNextToLowDomain(3) == FLASH_TRUE) .and. &
          (gr_domainBC(LOW,3) /= PERIODIC) ) kstart = NGUARD*K3D+1
     if ( (boxInfo % isNextToHighDomain(3) == FLASH_TRUE) .and. &
          (gr_domainBC(HIGH,3) /= PERIODIC) ) kend   = NGUARD*K3D+nzb
#endif

     do k = kstart,kend
        do j = jstart,jend
           do i = istart,iend

              ! d/dxdx
              delu2(1) = delu(1,i+1,j,k) - delu(1,i-1,j,k)
              delu2(1) = delu2(1)*delx

              delu3(1) = abs(delu(1,i+1,j,k)) + abs(delu(1,i-1,j,k))
              delu3(1) = delu3(1)*delx

              delu4(1) = delua(1,i+1,j,k) + delua(1,i-1,j,k)
              delu4(1) = delu4(1)*delx

#if N_DIM >= 2
              ! d/dydx
              delu2(2) = delu(1,i,j+1,k) - delu(1,i,j-1,k)
              delu2(2) = delu2(2)*dely_f

              delu3(2) = abs(delu(1,i,j+1,k)) + abs(delu(1,i,j-1,k))
              delu3(2) = delu3(2)*dely_f

              delu4(2) = delua(1,i,j+1,k) + delua(1,i,j-1,k)
              delu4(2) = delu4(2)*dely_f

              ! d/dxdy
              delu2(3) = delu(2,i+1,j,k) - delu(2,i-1,j,k)
              delu2(3) = delu2(3)*delx

              delu3(3) = abs(delu(2,i+1,j,k)) + abs(delu(2,i-1,j,k))
              delu3(3) = delu3(3)*delx

              delu4(3) = delua(2,i+1,j,k) + delua(2,i-1,j,k)
              delu4(3) = delu4(3)*delx

              ! d/dydy
              delu2(4) = delu(2,i,j+1,k) - delu(2,i,j-1,k)
              delu2(4) = delu2(4)*dely_f

              delu3(4) = abs(delu(2,i,j+1,k)) +  &
                   &                          abs(delu(2,i,j-1,k))
              delu3(4) = delu3(4)*dely_f

              delu4(4) = delua(2,i,j+1,k) + delua(2,i,j-1,k)
              delu4(4) = delu4(4)*dely_f
#endif

#if N_DIM == 3
              ! d/dzdx
              delu2(5) = delu(1,i,j,k+1) - delu(1,i,j,k-1)
              delu2(5) = delu2(5)*delz_f

              delu3(5) = abs(delu(1,i,j,k+1)) + abs(delu(1,i,j,k-1))
              delu3(5) = delu3(5)*delz_f

              delu4(5) = delua(1,i,j,k+1) + delua(1,i,j,k-1)
              delu4(5) = delu4(5)*delz_f

              ! d/dzdy
              delu2(6) = delu(2,i,j,k+1) - delu(2,i,j,k-1)
              delu2(6) = delu2(6)*delz_f

              delu3(6) = abs(delu(2,i,j,k+1)) + abs(delu(2,i,j,k-1))
              delu3(6) = delu3(6)*delz_f

              delu4(6) = delua(2,i,j,k+1) + delua(2,i,j,k-1)
              delu4(6) = delu4(6)*delz_f

              ! d/dxdz
              delu2(7) = delu(3,i+1,j,k) - delu(3,i-1,j,k)
              delu2(7) = delu2(7)*delx

              delu3(7) = abs(delu(3,i+1,j,k)) + abs(delu(3,i-1,j,k))
              delu3(7) = delu3(7)*delx

              delu4(7) = delua(3,i+1,j,k) + delua(3,i-1,j,k)
              delu4(7) = delu4(7)*delx

              ! d/dydz
              delu2(8) = delu(3,i,j+1,k) - delu(3,i,j-1,k)
              delu2(8) = delu2(8)*dely_f

              delu3(8) = abs(delu(3,i,j+1,k)) + abs(delu(3,i,j-1,k))
              delu3(8) = delu3(8)*dely_f

              delu4(8) = delua(3,i,j+1,k) + delua(3,i,j-1,k)
              delu4(8) = delu4(8)*dely_f

              ! d/dzdz
              delu2(9) = delu(3,i,j,k+1) - delu(3,i,j,k-1)
              delu2(9) = delu2(9)*delz_f

              delu3(9) = abs(delu(3,i,j,k+1)) + abs(delu(3,i,j,k-1))
              delu3(9) = delu3(9)*delz_f

              delu4(9) = delua(3,i,j,k+1) + delua(3,i,j,k-1)
              delu4(9) = delu4(9)*delz_f
#endif

              ! compute the error
              num = 0.
              denom = 0.

              do kk = 1, ndim2
                 num = num + delu2(kk)**2
                 denom = denom + (delu3(kk) + &
                      (refine_filter*delu4(kk)))**2
              end do

              ! mz -- compare the square of the error
              if (denom .eq. 0.0 .AND. num .ne. 0.0) then
                 error = HUGE(1.0)
              else if (denom .ne. 0.0) then
                 error = num/denom
              end if
              error = sqrt(error)

              if (error.gt.refine_cutoff) then
                 !We have found a cell that should be tagged.  This cell
                 !can be an internal cell or an external cell.  To ensure
                 !that we have advance warning about a shock we also tag
                 !a radius of cells around the cell at i,j,k.  We only tag
                 !those cells that exist within the internal region of a block.
                 do kr = k-K3D*gr_tagRadius, k+K3D*gr_tagRadius
                    do jr = j-K2D*gr_tagRadius, j+K2D*gr_tagRadius
                       do ir = i-K1D*gr_tagRadius, i+K1D*gr_tagRadius
                          if ( ir >= blkLimits(1,IAXIS) .and. &
                               ir <= blkLimits(2,IAXIS) .and. &
                               jr >= blkLimits(1,JAXIS) .and. &
                               jr <= blkLimits(2,JAXIS) .and. &
                               kr >= blkLimits(1,KAXIS) .and. &
                               kr <= blkLimits(2,KAXIS) ) then
                             ptrTagBlk(TAGC_SCRATCH_CENTER_VAR,ir,jr,kr) = FLASH_TRUE
                          end if
                       end do
                    end do
                 end do
              end if

           end do
        end do
     end do

     deallocate(delu)
     deallocate(delua)

     call Grid_releaseBlkPtr(blockID, ptrTagBlk, SCRATCH_CTR)
     call Grid_releaseBlkPtr(blockID, ptrBlk, CENTER)

  end do

end subroutine gr_markRefineDerefine
