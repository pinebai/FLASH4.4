!!****if* source/Grid/GridAllocatableScratches/Grid_ascModule
!!
!! NAME
!!  Grid_ascModule
!!
!! SYNOPSIS
!!
!!  use Grid_ascModule
!!
!! DESCRIPTION 
!!  
!!  Module with data and some accessor routines for allocatable scratches.
!!
!!   
!!***

!!REORDER(5):gr_ascArrayScratch, gr_ascArrayScratch_ctr, gr_ascArrayScratch_facevar[xyz]

Module Grid_ascModule

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, save :: hy_iguard = NGUARD
  integer, save :: hy_jguard = NGUARD 
  integer, save :: hy_kguard = NGUARD

#ifdef FIXEDBLOCKSIZE
  integer, save :: hy_ilo = GRID_ILO
  integer, save :: hy_ihi = GRID_IHI
  integer, save :: hy_jlo = GRID_JLO
  integer, save :: hy_jhi = GRID_JHI
  integer, save :: hy_klo = GRID_KLO
  integer, save :: hy_khi = GRID_KHI
#else
  integer, save :: hy_ilo = NGUARD+1
  integer, save :: hy_ihi
  integer, save :: hy_jlo = NGUARD*K2D+1
  integer, save :: hy_jhi
  integer, save :: hy_klo = NGUARD*K3D+1
  integer, save :: hy_khi
#endif

!  Type define

  integer,save,dimension(MAXBLOCKS) :: LeafNoToBId, BIdToLeafNo


  real,allocatable, target,dimension(:,:,:,:,:) :: gr_ascArrayScratch
  real,allocatable, target,dimension(:,:,:,:,:) :: gr_ascArrayCenter
  real,allocatable, target,dimension(:,:,:,:,:) :: gr_ascArrayScratch_ctr
  real,allocatable, target,dimension(:,:,:,:,:) :: gr_ascArrayScratch_facevarx
  real,allocatable, target,dimension(:,:,:,:,:) :: gr_ascArrayScratch_facevary
  real,allocatable, target,dimension(:,:,:,:,:) :: gr_ascArrayScratch_facevarz
  real,allocatable, target,dimension(:,:,:,:,:) :: gr_ascArrayFacevarx
  real,allocatable, target,dimension(:,:,:,:,:) :: gr_ascArrayFacevary
  real,allocatable, target,dimension(:,:,:,:,:) :: gr_ascArrayFacevarz

  real,allocatable, target,dimension(:,:,:,:,:,:) :: gr_ascArray5Center
  real,allocatable, target,dimension(:,:,:,:,:,:) :: gr_ascArray5Scratch_ctr





contains
  subroutine Grid_ascStart()
    use Grid_getBlkIndexLimits_mod, ONLY: Grid_getBlkIndexLimits
    use Driver_interface, ONLY: Driver_abortFlash
#ifdef DEBUG_ALL
    use Logfile_interface, ONLY: Logfile_stampMessage
#endif

    integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

    hy_iguard = NGUARD
    hy_jguard = NGUARD*K2D
    hy_kguard = NGUARD*K3D
    BidToLeafNo(:) = -1

#ifndef FIXEDBLOCKSIZE
    call Grid_getBlkIndexLimits(1,blkLimits,blkLimitsGC)
    hy_ihi = blkLimits(HIGH,IAXIS)
    hy_jhi = blkLimits(HIGH,JAXIS)
    hy_khi = blkLimits(HIGH,KAXIS)
#endif


  end subroutine Grid_ascStart

  subroutine Grid_ascAllocMem(gds,var1,nvars, nGuardCtr,nGuardFaceN,nGuardFaceT, &
                                leafBlocks, arrayRank, highSize)

    integer, intent(IN) :: gds,var1,nvars
    integer, OPTIONAL, intent(IN) :: nGuardCtr,nGuardFaceN,nGuardFaceT
    integer, OPTIONAL, intent(IN), dimension(:) :: leafBlocks
    integer, OPTIONAL, intent(IN) :: arrayRank
    integer, OPTIONAL, intent(IN) :: highSize

    integer :: lbnd(MDIM),ubnd(MDIM)
    integer :: lvar,    uvar
    integer :: nGuardCtrLoc,nGuardFaceNLoc,nGuardFaceTLoc
    integer :: nLeafBlocks,nBlocks
    integer :: rankLoc, highSizeLoc
    integer :: lb
    
    if (present(nGuardCtr)) then
       nGuardCtrLoc = nGuardCtr
    else
       nGuardCtrLoc = NGUARD
    end if
    if (present(nGuardFaceN)) then
       nGuardFaceNLoc = nGuardFaceN
    else
       nGuardFaceNLoc = NGUARD
    end if
    if (present(nGuardFaceT)) then
       nGuardFaceTLoc = nGuardFaceT
    else
       nGuardFaceTLoc = NGUARD
    end if

    if (present(leafBlocks)) then
       nLeafBlocks = size(leafBlocks)
       LeafNoToBid(1:nLeafBlocks) = leafBlocks(:)
       do lb=1,nLeafBlocks
          BidToLeafNo(leafBlocks(lb)) = lb
       end do
       nBlocks = nLeafBlocks
    else
       nLeafBlocks = -1
       do lb=1,MAXBLOCKS
          LeafNoToBid(lb) = lb
          BidToLeafNo(lb) = lb
       end do
       nBlocks = MAXBLOCKS
    end if

    if (present(arrayRank)) then
       rankLoc = arrayRank
    else
       rankLoc = 4
       if (present(highSize)) rankLoc = 5
    end if
    if (present(highSize)) then
       highSizeLoc = highSize
    else
       highSizeLoc = NDIM
    end if

#ifdef DEBUG_ALL
    if (rankLoc < 4 .OR. rankLoc > 5) then
       print*,'Grid_ascAllocMem: arrayRank must be 4 or 5.'
       call Driver_abortFlash('Grid_ascAllocMem: arrayRank musr be 4 or 5.')
    end if
    if (rankLoc == 5 .AND. gds .NE. CENTER .AND. gds .NE. SCRATCH_CTR) then
       print*,'Grid_ascAllocMem: arrayRank of 5 is only supported for gds=CENTER or SCRATCH_CTR.'
       call Driver_abortFlash('Grid_ascAllocMem: arrayRank of 5 is only supported for gds=CENTER.')
    end if
    if (rankLoc < 5 .AND. present(highSize)) then
       if (highSizeLoc .NE. 1) then
          print*,'Grid_ascAllocMem: arrayRank of 5 is only supported for gds=CENTER.'
          call Driver_abortFlash('Grid_ascAllocMem: arrayRank of 5 is only supported for gds=CENTER.')
       else
          call Logfile_stampMessage('[Grid_ascAllocMem] The highSize argument is superfluous.')
       end if
    end if
#endif


    lbnd(:) = (/hy_ilo,hy_jlo,hy_klo/)
    ubnd(:) = (/hy_ihi,hy_jhi,hy_khi/)

    lvar = var1
    uvar = var1 + nvars - 1


    select case (gds)
     case(SCRATCH)
        ubnd(1:NDIM) = ubnd(1:NDIM) + nGuardCtrLoc + 1
!!$        print *, 'TRIED TO GET SOMETHING OTHER THAN SCRATCH_CTR OR FACE[XYZ].'
!!$        call Driver_abortFlash("[Grid_ascModule] Unsupported gds=SCRATCH")
     case(CENTER)
        lbnd(1:NDIM) = lbnd(1:NDIM) - nGuardCtrLoc
        ubnd(1:NDIM) = ubnd(1:NDIM) + nGuardCtrLoc
     case(SCRATCH_CTR)
        lbnd(1:NDIM) = lbnd(1:NDIM) - nGuardCtrLoc
        ubnd(1:NDIM) = ubnd(1:NDIM) + nGuardCtrLoc
     case(SCRATCH_FACEX)
        lbnd(2:NDIM) = lbnd(2:NDIM) - nGuardFaceTLoc
        lbnd(1)      = lbnd(1)      - nGuardFaceNLoc
        ubnd(2:NDIM) = ubnd(2:NDIM) + nGuardFaceTLoc
        ubnd(1)      = ubnd(1)      + nGuardFaceNLoc !!+ 1
     case(SCRATCH_FACEY)
        lbnd(1:NDIM:2) = lbnd(1:NDIM:2) - nGuardFaceTLoc
        lbnd(2)        = lbnd(2)        - nGuardFaceNLoc
        ubnd(1:NDIM:2) = ubnd(1:NDIM:2) + nGuardFaceTLoc
        ubnd(2)        = ubnd(2)        + nGuardFaceNLoc !!+ K2D
     case(SCRATCH_FACEZ)
        lbnd(1:2)    = lbnd(1:2)    - nGuardFaceTLoc
        lbnd(3)      = lbnd(3)      - nGuardFaceNLoc
        ubnd(1:2)    = ubnd(1:2)    + nGuardFaceTLoc
        ubnd(3)      = ubnd(3)      + nGuardFaceNLoc !!+ K3D
     case(FACEX)
        lbnd(2:NDIM) = lbnd(2:NDIM) - nGuardFaceTLoc
        lbnd(1)      = lbnd(1)      - nGuardFaceNLoc
        ubnd(2:NDIM) = ubnd(2:NDIM) + nGuardFaceTLoc
        ubnd(1)      = ubnd(1)      + nGuardFaceNLoc + 1
     case(FACEY)
        lbnd(1:NDIM:2) = lbnd(1:NDIM:2) - nGuardFaceTLoc
        lbnd(2)        = lbnd(2)        - nGuardFaceNLoc
        ubnd(1:NDIM:2) = ubnd(1:NDIM:2) + nGuardFaceTLoc
        ubnd(2)        = ubnd(2)        + nGuardFaceNLoc + K2D
     case(FACEZ)
        lbnd(1:2)    = lbnd(1:2)    - nGuardFaceTLoc
        lbnd(3)      = lbnd(3)      - nGuardFaceNLoc
        ubnd(1:2)    = ubnd(1:2)    + nGuardFaceTLoc
        ubnd(3)      = ubnd(3)      + nGuardFaceNLoc + K3D
     case DEFAULT
        print *, 'TRIED TO GET SOMETHING OTHER THAN SCRATCH_CTR OR FACE[XYZ].'
        call Driver_abortFlash("[Grid_ascModule] Unsupported gds")
    end select


    if (rankLoc==4) then 
     select case (gds)
     case(SCRATCH)
        allocate(gr_ascArrayScratch(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(CENTER)
        allocate(gr_ascArrayCenter(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(SCRATCH_CTR)
        allocate(gr_ascArrayScratch_ctr(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(SCRATCH_FACEX)
        allocate(gr_ascArrayScratch_facevarx(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(SCRATCH_FACEY)
        allocate(gr_ascArrayScratch_facevary(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(SCRATCH_FACEZ)
        allocate(gr_ascArrayScratch_facevarz(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(FACEX)
        allocate(gr_ascArrayFacevarx(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(FACEY)
        allocate(gr_ascArrayFacevary(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(FACEZ)
        allocate(gr_ascArrayFacevarz(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case DEFAULT
        call Driver_abortFlash("[Grid_ascAllocMem] Unsupported gds")
     end select
    else if (rankLoc==5) then
    select case (gds)
!!$     case(SCRATCH)
!!$        allocate(gr_ascArray5Scratch(lvar:uvar, &
!!$             lbnd(1):ubnd(1), &
!!$             lbnd(2):ubnd(2), &
!!$             lbnd(3):ubnd(3), &
!!$             highSizeLoc, nBlocks))
     case(CENTER)
        allocate(gr_ascArray5Center(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             highSizeLoc, nBlocks))
     case(SCRATCH_CTR)
        allocate(gr_ascArray5Scratch_ctr(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             highSizeLoc, nBlocks))
!!$     case(SCRATCH_FACEX)
!!$        allocate(gr_ascArray5Scratch_facevarx(lvar:uvar, &
!!$             lbnd(1):ubnd(1), &
!!$             lbnd(2):ubnd(2), &
!!$             lbnd(3):ubnd(3), &
!!$             highSizeLoc, nBlocks))
!!$     case(SCRATCH_FACEY)
!!$        allocate(gr_ascArray5Scratch_facevary(lvar:uvar, &
!!$             lbnd(1):ubnd(1), &
!!$             lbnd(2):ubnd(2), &
!!$             lbnd(3):ubnd(3), &
!!$             highSizeLoc, nBlocks))
!!$     case(SCRATCH_FACEZ)
!!$        allocate(gr_ascArray5Scratch_facevarz(lvar:uvar, &
!!$             lbnd(1):ubnd(1), &
!!$             lbnd(2):ubnd(2), &
!!$             lbnd(3):ubnd(3), &
!!$             highSizeLoc, nBlocks))
     case DEFAULT
        call Driver_abortFlash("[Grid_ascAllocMem] Unsupported gds for arrayRank=5")
     end select
    end if
  end subroutine Grid_ascAllocMem

  subroutine Grid_ascDeallocMem(gds, arrayRank)
    use Driver_interface, ONLY: Driver_abortFlash

    integer, OPTIONAL, intent(IN) :: gds
    integer, OPTIONAL, intent(IN) :: arrayRank

    integer :: rankLoc

    if (present(arrayRank)) then
       rankLoc = arrayRank
    else
       rankLoc = 0
    end if

    if (present(gds)) then
       if (rankLoc==0 .OR. rankLoc==4) then
          select case (gds)
          case(SCRATCH)
             deallocate(gr_ascArrayScratch)
          case(CENTER)
             if (allocated(gr_ascArrayCenter)) deallocate(gr_ascArrayCenter)
          case(SCRATCH_CTR)
             deallocate(gr_ascArrayScratch_ctr)
          case(SCRATCH_FACEX)
             deallocate(gr_ascArrayScratch_facevarx)
          case(SCRATCH_FACEY)
             deallocate(gr_ascArrayScratch_facevary)
          case(SCRATCH_FACEZ)
             deallocate(gr_ascArrayScratch_facevarz)
          case(FACEX)
             deallocate(gr_ascArrayFacevarx)
          case(FACEY)
             deallocate(gr_ascArrayFacevary)
          case(FACEZ)
             deallocate(gr_ascArrayFacevarz)
          case DEFAULT
             call Driver_abortFlash("[Grid_ascDeallocMem] Unsupported gds")
          end select
       end if
       if (rankLoc==5) then
          select case (gds)
          case(CENTER)
             deallocate(gr_ascArray5Center)
          case(SCRATCH_CTR)
             deallocate(gr_ascArray5Scratch_ctr)
          case DEFAULT
             call Driver_abortFlash("[Grid_ascDeallocMem] Unsupported gds for arrayRank=5")
          end select
       else if (rankLoc==0) then
          select case (gds)
          case(CENTER)
             if (allocated(gr_ascArray5Center)) deallocate(gr_ascArray5Center)
          case(SCRATCH_CTR)
             if (allocated(gr_ascArray5Scratch_ctr)) deallocate(gr_ascArray5Scratch_ctr)
          end select
       end if
    else
       if (rankLoc==0 .OR. rankLoc==4) then
          if (allocated(gr_ascArrayScratch))          deallocate(gr_ascArrayScratch)
          if (allocated(gr_ascArrayCenter))      deallocate(gr_ascArrayCenter)
          if (allocated(gr_ascArrayScratch_ctr))      deallocate(gr_ascArrayScratch_ctr)
          if (allocated(gr_ascArrayScratch_facevarx)) deallocate(gr_ascArrayScratch_facevarx)
          if (allocated(gr_ascArrayScratch_facevary)) deallocate(gr_ascArrayScratch_facevary)
          if (allocated(gr_ascArrayScratch_facevarz)) deallocate(gr_ascArrayScratch_facevarz)
          if (allocated(gr_ascArrayFacevarx))    deallocate(gr_ascArrayFacevarx)
          if (allocated(gr_ascArrayFacevary))    deallocate(gr_ascArrayFacevary)
          if (allocated(gr_ascArrayFacevarz))    deallocate(gr_ascArrayFacevarz)
       end if
       if (rankLoc==0 .OR. rankLoc==5) then
          if (allocated(gr_ascArray5Center))      deallocate(gr_ascArray5Center)
          if (allocated(gr_ascArray5Scratch_ctr))      deallocate(gr_ascArray5Scratch_ctr)
       end if
    end if

  end subroutine Grid_ascDeallocMem

end Module Grid_ascModule
