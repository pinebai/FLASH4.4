!!****if* source/physics/sourceTerms/Turb/TurbMain/Turb_calc
!!
!! NAME
!!
!!  Turb_calc
!!
!! SYNOPSIS
!!
!!  call Turb_calc(integer, INTENT(in)  :: num_blocks,
!!                 integer, INTENT(in), DIMENSION(num_blocks)  :: blocklist)
!!
!! DESCRIPTION
!!
!! Aaron Jackson 2010
!! Calculate the laplacian of the velocity field and store
!! in TURB_VAR, LAPY_VAR, and LAPZ_VAR for each block, then
!! fill guard cells.
!! Calculate the curl of the laplacian for each block, then
!! fill guard cells again.
!! This routine is OP2 in Colin et al. (2000)
!!
!!
!! ARGUMENTS
!!
!!   num_blocks : number of blocks 
!!
!!   blocklist : list of blocks
!!
!!
!!
!!***

#include "Flash.h"
#include "constants.h"
subroutine Turb_calc(num_blocks, blockList)    

  use Turb_data, only : turb_useTurb, turb_stepSize, turb_c2, turb_fillGC
  use Grid_interface, only : Grid_fillGuardCells, Grid_getBlkPtr, &
        Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getDeltas, &
        Grid_getMaxRefinement, &
        Grid_getBlkRefineLevel
  use Turb_interface, only : Turb_laplacian, Turb_curlMag
  use Driver_interface, only : Driver_abortFlash
  use Timers_interface, only : Timers_start, Timers_stop
       
  implicit none
  integer, INTENT(in)                        :: num_blocks
  integer, INTENT(in), DIMENSION(num_blocks) :: blockList

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: istat
  integer :: n, bid, blevel, maxlevel

  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:,:,:) :: laplX, laplY, laplZ, &
                                         curl, lapl, velComp
  integer :: tmp_gcMaskSize
  logical, dimension(NUNK_VARS) :: tmp_gcMask
  logical :: tmp_gcDoEos, tmp_gcMakeConsist

  real, dimension(MDIM) :: deltas

  integer :: sizeI, sizeJ, sizeK
  integer :: local_stepSize, dlevel
  integer :: i,j,k

  if (.NOT. turb_useTurb) return  !  return immediately if disabled

  call Timers_start("turb")

#if defined(FLASH_GRID_PARAMESH)
  call Grid_getMaxRefinement(maxlevel)
#else
  maxlevel = 1
#endif

  tmp_gcMaskSize = NUNK_VARS
  ! fill guard cells
  tmp_gcDoEos = .false.
  tmp_gcMakeConsist = .true.
  tmp_gcMask(:) = .false.
  tmp_gcMask(VELX_VAR) = .true.
  if (NDIM >= 2) tmp_gcMask(VELY_VAR) = .true.
  if (NDIM > 2) tmp_gcMask(VELZ_VAR) = .true.

  call Grid_fillGuardCells(CENTER, ALLDIR, doEos=tmp_gcDoEos, &
                      maskSize=tmp_gcMaskSize, mask=tmp_gcMask, &
                      makeMaskConsistent=tmp_gcMakeConsist)

  ! calculate the laplacian for all blocks
  do n = 1,num_blocks
     bid = blockList(n)

     call Grid_getBlkPtr(bid,solnData)
     call Grid_getBlkRefineLevel(bid,blevel)

     ! make sure we don't divide by zero
     dlevel = maxlevel - blevel
     local_stepSize = turb_stepSize / 2**dlevel

     ! if local_stepSize is zero, then we are at a refinement level
     ! too course to be used by the flame.
     if (local_stepSize > 0) then

        call Grid_getBlkIndexLimits(bid, blkLimits, blkLimitsGC)

        sizeI=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        sizeJ=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        sizeK=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

        allocate(velComp(sizeI,sizeJ,sizeK),STAT=istat)
        if (istat /= 0) call Driver_abortFlash("Cannot allocate velComp in Turb_calc")
        allocate(lapl(sizeI,sizeJ,sizeK), STAT=istat)
        if (istat /= 0) call Driver_abortFlash("Cannot allocate lapl in Turb_calc")
        lapl(:,:,:) = 0.e0

        velComp( blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) ) =  &
                    solnData( VELX_VAR,                          &
                           blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                           blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                           blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )

        call Turb_laplacian(lapl, velComp, local_stepSize, bid)

        ! we use TURB_VAR to hold LAPX
        ! store lapl in TURB_VAR
        solnData(TURB_VAR,  &
                 blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )   = &
                    lapl( blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )
        if (NDIM >= 2) then 
           velComp( blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                    blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                    blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) ) =  &
                    solnData( VELY_VAR,                             &
                           blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                           blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                           blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )

           call Turb_laplacian(lapl, velComp, local_stepSize, bid)
           solnData(LAPY_VAR,  &
                    blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                    blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                    blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )   = &
                       lapl( blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )
        endif
        if (NDIM > 2) then
           velComp( blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                    blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                    blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) ) =  &
                    solnData( VELZ_VAR,                             &
                              blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                              blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                              blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )

           call Turb_laplacian(lapl, velComp, local_stepSize, bid)
           solnData(LAPZ_VAR,  &
                    blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                    blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                    blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )   = &
                       lapl( blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )

        endif

        deallocate(velComp)
        deallocate(lapl)

     endif

     call Grid_releaseBlkPtr(bid, solnData)

  enddo

  ! see if guard cell filling is necessary
  ! this assumes a 5-point stencil and 2 operators
  if (turb_fillGC) then
     ! fill guard cells
     tmp_gcDoEos = .false.
     tmp_gcMask(:) = .false.
     tmp_gcMask(TURB_VAR) = .true.
     if (NDIM >= 2) tmp_gcMask(LAPY_VAR) = .true.
     if (NDIM > 2) tmp_gcMask(LAPZ_VAR) = .true.

     call Grid_fillGuardCells(CENTER, ALLDIR, doEos=tmp_gcDoEos, &
                         maskSize=tmp_gcMaskSize, mask=tmp_gcMask, &
                         makeMaskConsistent=tmp_gcMakeConsist)

  endif

  do n = 1,num_blocks
     bid = blockList(n)

     call Grid_getBlkPtr(bid,solnData)
     call Grid_getBlkIndexLimits(bid, blkLimits, blkLimitsGC)
     call Grid_getBlkRefineLevel(bid,blevel)
     call Grid_getDeltas(bid, deltas)

     ! only need to perform the second operator at the flame level
     if (blevel .eq. maxlevel) then

        sizeI=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        sizeJ=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        sizeK=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

        allocate(laplX(sizeI,sizeJ,sizeK), STAT=istat)
        if (istat /= 0) call Driver_abortFlash("Cannot allocate laplX in Flame_step")
        allocate(laplY(sizeI,sizeJ,sizeK), STAT=istat)
        if (istat /= 0) call Driver_abortFlash("Cannot allocate laplY in Flame_step")
        allocate(laplZ(sizeI,sizeJ,sizeK), STAT=istat)
        if (istat /= 0) call Driver_abortFlash("Cannot allocate laplZ in Flame_step")
        allocate(curl(sizeI,sizeJ,sizeK),STAT=istat)
        if (istat /= 0) call Driver_abortFlash("Cannot allocate curl in Flame_step")
        curl(:,:,:) = 0.e0

        laplX( blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )   = &
                     solnData(TURB_VAR,  &
                          blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )
        if (NDIM >= 2) then
           laplY( blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )   = &
                        solnData(LAPY_VAR,  &
                             blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )
        else
           laplY(:,:,:) = 0.e0
        endif
        if (NDIM > 2) then
           laplZ( blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )   = &
                         solnData(LAPZ_VAR,  &
                             blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                             blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                             blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) )
        else
           laplZ(:,:,:) = 0.e0
        endif

        call Turb_curlMag(curl, laplX, laplY, laplZ, turb_stepSize, bid)

        ! calculate OP2 from Colin et al. 2000

        curl( blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
              blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
              blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS) )   =  &
                 turb_c2 * (turb_stepSize * deltas(IAXIS))**3 *  &
                 curl( blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS) )

        ! store u' on Grid using TURB
        solnData(TURB_VAR,  &
                 blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                 blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                 blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS) )   = &
                    curl( blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),  &
                          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),  &
                          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS) )

        deallocate(laplX)
        deallocate(laplY)
        deallocate(laplZ)
        deallocate(curl)

     else

        ! set TURB_VAR to zero
        solnData(TURB_VAR,  &
                 blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), &
                 blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), &
                 blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS) ) = 0.0e0

     endif

     call Grid_releaseBlkPtr(bid, solnData)

  enddo

  call Timers_stop("turb")

  return
end subroutine Turb_calc
