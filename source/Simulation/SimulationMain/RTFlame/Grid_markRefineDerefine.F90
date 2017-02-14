!!****if* source/Simulation/SimulationMain/RTFlame/Grid_markRefineDerefine
!!
!! NAME
!!
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Grid_markRefineDerefine()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! AUTOGENROBODOC
!!
!!
!!***

! Dean Townsley 2009
!
! See source/Grid/Grid_markRefineDerefine.F90
!     for API documentation
! and source/Grid/GridMain/paramesh/Grid_markRefineDerefine.F90
!     for example implementation
!
! This version is specific for the RT Flame in a channel simulation
!   two options are selected by value of refine_uniform_region parameter
!   (1) if .false. act like standard refinement
!   (2) if .true. refine a region of the selected size to full
!       refinement (see diagram below)


!#define KEEP_REDUNDANT_EOS_CALLS

subroutine Grid_markRefineDerefine()

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var
  use tree, ONLY : newchild, refine, derefine, stay, lrefine, &
                   lrefine_max, lrefine_min
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getListOfBlocks, &
    Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getCellCoords, &
    Grid_fillGuardCells, Grid_getBlkBoundBox
  use Simulation_data, ONLY : sim_refine_uniform_region, &
       sim_refine_lead, sim_refine_region_size, sim_refine_buf, &
       sim_refine_region_stepdown_size, &
       sim_ParticleRefineRegion, sim_ParticleRefineRegionLevel, &
       sim_ParticleRefineRegionBottom, sim_ParticleRefineRegionTop
  use Eos_interface, ONLY: Eos_wrapped

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"

  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref

  logical :: gcMask(NUNK_VARS) ! indicates for which variables the refinement criteria
  ! that are implemented in this file need values in guardcells
  
  integer, dimension(2,MDIM) :: blkLimits
  integer, dimension(2,MDIM) :: blkLimitsGC
  integer ,dimension(MAXBLOCKS) :: blkList
  integer :: blkCount, iblk, isizeGC, j, k, error
  real, dimension(2,MDIM) :: boundBox

  real, dimension(:,:,:,:), pointer :: solnData

  real, dimension(MAXBLOCKS) :: blk_min_x, blk_max_x
  real :: flame_local_max_x, flame_max_x

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC) :: iCenter
#else
  integer :: istat
  real, allocatable, dimension(:) :: iCenter
#endif

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  if (.not. sim_refine_uniform_region ) then
     !----------------
     ! Standard refinement conditions

     gcMask = .FALSE.
     do i = 1,gr_numRefineVars
        iref = gr_refine_var(i)
        if (iref > 0) gcMask(iref) = .TRUE.
     end do

     if(any(gcMask)) then
        call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.true.,&
             maskSize=NUNK_VARS, mask=gcMask, makeMaskConsistent=.true.)
     end if

#ifdef KEEP_REDUNDANT_EOS_CALLS
     call Grid_getListOfBlocks(LEAF, blkList,blkCount)
     do iblk = 1, blkCount
        call Grid_getBlkIndexLimits(blkList(iblk), blkLimits, blkLimitsGC)
        call Eos_wrapped(MODE_DENS_EI,blkLimits, blkList(iblk))
     end do
#endif

     do l = 1,gr_numRefineVars
        iref = gr_refine_var(l)
        ref_cut = gr_refine_cutoff(l)
        deref_cut = gr_derefine_cutoff(l)
        ref_filter = gr_refine_filter(l)
        call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
     end do

#ifdef FLASH_GRID_PARAMESH2
     ! Make sure lrefine_min and lrefine_max are obeyed - KW
     if (gr_numRefineVars .LE. 0) then
        call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
     end if
#endif

  else

     !--------------------------------------
     ! put down a region of uniform refinement near the flame

     call Grid_getListOfBlocks(LEAF, blkList,blkCount)

     !-----------------------------------
     !  Find max height of burned material
     !-----------------------------------
     flame_local_max_x = -HUGE(1.0)
     do iblk = 1, blkCount
        call Grid_getBlkIndexLimits(blkList(iblk), blkLimits, blkLimitsGC)

        isizeGC=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1

#ifndef FIXEDBLOCKSIZE
        allocate(iCenter(isizeGC),STAT=istat)
        if (istat /= 0) print *,' ERROR allocating iCenter in Grid_markRefineDerefine'
#endif
        call Grid_getCellCoords(IAXIS, blkList(iblk), CENTER, .true., iCenter, isizeGC)

        blk_min_x(blkList(iblk)) = iCenter(blkLimits(LOW,IAXIS))
        blk_max_x(blkList(iblk)) = iCenter(blkLimits(HIGH,IAXIS))

        call Grid_getBlkPtr(blkList(iblk),solnData,CENTER)

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 if ( solnData(FLAM_MSCALAR,i,j,k) > 0.01  .and. iCenter(i) > flame_local_max_x ) flame_local_max_x = iCenter(i)
              enddo
           enddo
        enddo

        call Grid_releaseBlkPtr(blkList(iblk),solnData)

#ifndef FIXEDBLOCKSIZE
        deallocate(iCenter,STAT=istat)
        if (istat /= 0) print *,' ERROR deallocating iCenter in Grid_markRefineDerefine'
#endif
     enddo

     ! now take max over whole domain
     call MPI_AllReduce(flame_local_max_x, flame_max_x, 1, MPI_Double_Precision, MPI_Max, MPI_Comm_World, error)

     !-----------------------------------
     !  Mark refined region 
     !
     !                        <buf>                              <--- refine_lead ---> <buf>
     !------------------------------------------------------------------------
     !                       |     | XXXXXXXXXXXXXXXXXXXXXXXXXX | XXXXXXXXXXXXXXXXXXX |     |
     !   |    |xxxxxxxxxxxxxxxxxxxx|                       flame_max_x
     !   <---stepdown_size--> <------------------ refine_region_size ---------------->
     !
     !  We enforce the entire region inside the buf's to be fully refined, but only
     !  derefine blocks fully outside of the buffer
     !  the assumption is that flame_max_x is approximately moving steadily upward

     do iblk = 1, blkCount
        ! refine if in refine region and not at max refine
        if ( blk_min_x(blkList(iblk)) <= flame_max_x + sim_refine_lead .and. &
             blk_max_x(blkList(iblk)) >  flame_max_x + sim_refine_lead - sim_refine_region_size + sim_refine_buf .and. &
             lrefine(blkList(iblk)) < lrefine_max ) then
            refine(blkList(iblk)) = .true.
        endif
        ! refine if in stepdown region and not at max refine - 1
        if ( blk_min_x(blkList(iblk)) <= flame_max_x + sim_refine_lead - sim_refine_region_size .and. &
             blk_max_x(blkList(iblk)) >  flame_max_x + sim_refine_lead - sim_refine_region_size &
                                             - sim_refine_region_stepdown_size + sim_refine_buf .and. &
             lrefine(blkList(iblk)) < lrefine_max-1 ) then
            refine(blkList(iblk)) = .true.
        endif
        ! derefine if outside (behind) refine region and above max refine - 1
        if ( ( blk_max_x(blkList(iblk)) <=  flame_max_x + sim_refine_lead - sim_refine_region_size ) .and. &
             lrefine(blkList(iblk)) > lrefine_max-1 ) then
           derefine(blkList(iblk)) = .true.
        endif
        ! derefine completely if outside both refine and stepdown region
        if ( ( blk_min_x(blkList(iblk)) >  flame_max_x + sim_refine_lead + sim_refine_buf .or. &
               blk_max_x(blkList(iblk)) <=  flame_max_x + sim_refine_lead &
                                             - sim_refine_region_size - sim_refine_region_stepdown_size ) .and. &
             lrefine(blkList(iblk)) > lrefine_min ) then
            derefine(blkList(iblk)) = .true.
        endif
     enddo

  endif

  if (sim_ParticleRefineRegion) then
     !--------------------------------------
     ! refinement switches to allow accounting for particle load
     ! just refines to level sim_ParticleRefineRegionLevel over the region
     ! between sim_ParticleRefineRegionBottom and sim_ParticleRefineRegionTop
     call Grid_getListOfBlocks(LEAF, blkList,blkCount)
     do iblk = 1, blkCount
        call Grid_getBlkBoundBox(blkList(iblk), boundBox)
        if ( boundBox(2,IAXIS) >= sim_ParticleRefineRegionBottom .and. &
             boundBox(1,IAXIS) <= sim_ParticleRefineRegionTop) then
             if ( lrefine(blkList(iblk)) < sim_ParticleRefineRegionLevel .and. &
                  lrefine(blkList(iblk)) < lrefine_max) then
                 refine(blkList(iblk)) = .true.
                 derefine(blkList(iblk)) = .false.
             elseif (  lrefine(blkList(iblk)) == sim_ParticleRefineRegionLevel .and. &
                       lrefine(blkList(iblk)) <= lrefine_max ) then
                derefine(blkList(iblk)) = .false.
             endif
        endif
     enddo
  endif

  return
end subroutine Grid_markRefineDerefine














