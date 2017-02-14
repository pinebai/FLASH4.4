!!****if* source/physics/Hydro/HydroMain/hy_memScratchData
!!
!! NAME
!!  hy_memScratchData
!!
!! SYNOPSIS
!!
!!  use hy_memScratchData
!!
!! DESCRIPTION 
!!  
!!  Data and some basic routines for allocatable scratch arrays private to the Hydro unit
!!  
!!  
!! 
!!***

!!REORDER(5):hy_memArrayScratch, hy_memArrayScratch_ctr, hy_memArrayScratch_facevar[xyz]

Module hy_memScratchData

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, save :: hy_iguard = NGUARD
  integer, save :: hy_jguard = NGUARD 
  integer, save :: hy_kguard = NGUARD

#ifdef FIXEDBLOCKSIZE
  integer, save :: hy_iloGc = GRID_ILO_GC
  integer, save :: hy_ihiGc = GRID_IHI_GC
  integer, save :: hy_jloGc = GRID_JLO_GC
  integer, save :: hy_jhiGc = GRID_JHI_GC
  integer, save :: hy_kloGc = GRID_KLO_GC
  integer, save :: hy_khiGc = GRID_KHI_GC

  integer, save :: hy_ilo = GRID_ILO
  integer, save :: hy_ihi = GRID_IHI
  integer, save :: hy_jlo = GRID_JLO
  integer, save :: hy_jhi = GRID_JHI
  integer, save :: hy_klo = GRID_KLO
  integer, save :: hy_khi = GRID_KHI
#else
  integer, save :: hy_iloGc = 1
  integer, save :: hy_ihiGc
  integer, save :: hy_jloGc = 1
  integer, save :: hy_jhiGc
  integer, save :: hy_kloGc = 1
  integer, save :: hy_khiGc

  integer, save :: hy_ilo = NGUARD+1
  integer, save :: hy_ihi
  integer, save :: hy_jlo = NGUARD*K2D+1
  integer, save :: hy_jhi
  integer, save :: hy_klo = NGUARD*K3D+1
  integer, save :: hy_khi
#endif

!  Type define

  integer,save,dimension(MAXBLOCKS) :: LeafNoToBId, BIdToLeafNo


  real,allocatable, target,dimension(:,:,:,:,:) :: hy_memArrayScratch
  real,allocatable, target,dimension(:,:,:,:,:) :: hy_memArrayCenter
  real,allocatable, target,dimension(:,:,:,:,:) :: hy_memArrayScratch_ctr
  real,allocatable, target,dimension(:,:,:,:,:) :: hy_memArrayScratch_facevarx
  real,allocatable, target,dimension(:,:,:,:,:) :: hy_memArrayScratch_facevary
  real,allocatable, target,dimension(:,:,:,:,:) :: hy_memArrayScratch_facevarz

  real,allocatable, target,dimension(:,:,:,:,:,:) :: hy_memArray5Center
  real,allocatable, target,dimension(:,:,:,:,:,:) :: hy_memArray5Scratch_ctr


  real ,save :: hy_imin,hy_imax,hy_jmin,hy_jmax,hy_kmin,hy_kmax




contains
  subroutine hy_memScratchInitialize()
    use Grid_interface, ONLY: Grid_getBlkIndexLimits
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
    hy_ihiGc = blkLimitsGC(HIGH,IAXIS)
    hy_jhiGc = blkLimitsGC(HIGH,JAXIS)
    hy_khiGc = blkLimitsGC(HIGH,KAXIS)
#endif


  end subroutine hy_memScratchInitialize

  subroutine hy_memAllocScratch(gds,var1,nvars, nGuardCtr,nGuardFaceN,nGuardFaceT, &
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
       print*,'hy_memAllocScratch: arrayRank must be 4 or 5.'
       call Driver_abortFlash('hy_memAllocScratch: arrayRank musr be 4 or 5.')
    end if
    if (rankLoc == 5 .AND. gds .NE. CENTER .AND. gds .NE. SCRATCH_CTR) then
       print*,'hy_memAllocScratch: arrayRank of 5 is only supported for gds=CENTER or SCRATCH_CTR.'
       call Driver_abortFlash('hy_memAllocScratch: arrayRank of 5 is only supported for gds=CENTER.')
    end if
    if (rankLoc < 5 .AND. present(highSize)) then
       if (highSizeLoc .NE. 1) then
          print*,'hy_memAllocScratch: arrayRank of 5 is only supported for gds=CENTER.'
          call Driver_abortFlash('hy_memAllocScratch: arrayRank of 5 is only supported for gds=CENTER.')
       else
          call Logfile_stampMessage('[hy_memAllocScratch] The highSize arguemnt is superfluous.')
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
!!$        call Driver_abortFlash("[hy_memScratchData] Unsupported gds=SCRATCH")
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
     case DEFAULT
        print *, 'TRIED TO GET SOMETHING OTHER THAN SCRATCH_CTR OR FACE[XYZ].'
        call Driver_abortFlash("[hy_memScratchData] Unsupported gds")
    end select


    if (rankLoc==4) then 
     select case (gds)
     case(SCRATCH)
        allocate(hy_memArrayScratch(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(CENTER)
        allocate(hy_memArrayCenter(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(SCRATCH_CTR)
        allocate(hy_memArrayScratch_ctr(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(SCRATCH_FACEX)
        allocate(hy_memArrayScratch_facevarx(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(SCRATCH_FACEY)
        allocate(hy_memArrayScratch_facevary(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case(SCRATCH_FACEZ)
        allocate(hy_memArrayScratch_facevarz(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             nBlocks))
     case DEFAULT
        call Driver_abortFlash("[hy_memAllocScratch] Unsupported gds")
     end select
    else if (rankLoc==5) then
    select case (gds)
!!$     case(SCRATCH)
!!$        allocate(hy_memArray5Scratch(lvar:uvar, &
!!$             lbnd(1):ubnd(1), &
!!$             lbnd(2):ubnd(2), &
!!$             lbnd(3):ubnd(3), &
!!$             highSizeLoc, nBlocks))
     case(CENTER)
        allocate(hy_memArray5Center(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             highSizeLoc, nBlocks))
     case(SCRATCH_CTR)
        allocate(hy_memArray5Scratch_ctr(lvar:uvar, &
             lbnd(1):ubnd(1), &
             lbnd(2):ubnd(2), &
             lbnd(3):ubnd(3), &
             highSizeLoc, nBlocks))
!!$     case(SCRATCH_FACEX)
!!$        allocate(hy_memArray5Scratch_facevarx(lvar:uvar, &
!!$             lbnd(1):ubnd(1), &
!!$             lbnd(2):ubnd(2), &
!!$             lbnd(3):ubnd(3), &
!!$             highSizeLoc, nBlocks))
!!$     case(SCRATCH_FACEY)
!!$        allocate(hy_memArray5Scratch_facevary(lvar:uvar, &
!!$             lbnd(1):ubnd(1), &
!!$             lbnd(2):ubnd(2), &
!!$             lbnd(3):ubnd(3), &
!!$             highSizeLoc, nBlocks))
!!$     case(SCRATCH_FACEZ)
!!$        allocate(hy_memArray5Scratch_facevarz(lvar:uvar, &
!!$             lbnd(1):ubnd(1), &
!!$             lbnd(2):ubnd(2), &
!!$             lbnd(3):ubnd(3), &
!!$             highSizeLoc, nBlocks))
     case DEFAULT
        call Driver_abortFlash("[hy_memAllocScratch] Unsupported gds for arrayRank=5")
     end select
    end if
  end subroutine hy_memAllocScratch

  subroutine hy_memDeallocScratch(gds, arrayRank)
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
             deallocate(hy_memArrayScratch)
          case(CENTER)
             if (allocated(hy_memArrayCenter)) deallocate(hy_memArrayCenter)
          case(SCRATCH_CTR)
             deallocate(hy_memArrayScratch_ctr)
          case(SCRATCH_FACEX)
             deallocate(hy_memArrayScratch_facevarx)
          case(SCRATCH_FACEY)
             deallocate(hy_memArrayScratch_facevary)
          case(SCRATCH_FACEZ)
             deallocate(hy_memArrayScratch_facevarz)
          case DEFAULT
             call Driver_abortFlash("[hy_memDeallocScratch] Unsupported gds")
          end select
       end if
       if (rankLoc==5) then
          select case (gds)
          case(CENTER)
             deallocate(hy_memArray5Center)
          case(SCRATCH_CTR)
             deallocate(hy_memArray5Scratch_ctr)
          case DEFAULT
             call Driver_abortFlash("[hy_memDeallocScratch] Unsupported gds for arrayRank=5")
          end select
       else if (rankLoc==0) then
          select case (gds)
          case(CENTER)
             if (allocated(hy_memArray5Center)) deallocate(hy_memArray5Center)
          case(SCRATCH_CTR)
             if (allocated(hy_memArray5Scratch_ctr)) deallocate(hy_memArray5Scratch_ctr)
          end select
       end if
    else
       if (rankLoc==0 .OR. rankLoc==4) then
          if (allocated(hy_memArrayScratch))          deallocate(hy_memArrayScratch)
          if (allocated(hy_memArrayCenter))      deallocate(hy_memArrayCenter)
          if (allocated(hy_memArrayScratch_ctr))      deallocate(hy_memArrayScratch_ctr)
          if (allocated(hy_memArrayScratch_facevarx)) deallocate(hy_memArrayScratch_facevarx)
          if (allocated(hy_memArrayScratch_facevary)) deallocate(hy_memArrayScratch_facevary)
          if (allocated(hy_memArrayScratch_facevarz)) deallocate(hy_memArrayScratch_facevarz)
       end if
       if (rankLoc==0 .OR. rankLoc==5) then
          if (allocated(hy_memArray5Center))      deallocate(hy_memArray5Center)
          if (allocated(hy_memArray5Scratch_ctr))      deallocate(hy_memArray5Scratch_ctr)
       end if
    end if

  end subroutine hy_memDeallocScratch

end Module hy_memScratchData
