!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhInit
!!
!! NAME
!!
!!  gr_bhInit
!!
!! 
!! SYNOPSIS
!!
!!  call gr_bhInit()
!!
!!
!! DESCRIPTION
!!
!!  Initialize the tree Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the gr_bhData module.
!!
!!
!!
!!***

subroutine gr_bhInit()

  use Grid_data, ONLY : gr_geometry, gr_isolatedBoundaries, &
       gr_meshComm, gr_meshMe, gr_meshNumProcs
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getMinCellSizes, &
    Grid_getListOfBlocks, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Gravity_interface, ONLY : Gravity_bhGetNodeStruct
  use TreeRay_interface, ONLY : TreeRay_bhGetNodeStruct
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use gr_bhData, ONLY : gr_bhTreeLevels, &
    gr_bhTreeArray, gr_bhTreeBS, gr_bhTreeLimangle2, &
    gr_bhTreeSafeBox, gr_bhTreeSafeBoxHalf2, &
    gr_bhLocCoords, gr_bhTreeLoff, gr_bhTreeMyPE, &
    gr_bhTreeBCen, gr_bhTreeLrefine, gr_bhTreeCellSize, &
    gr_bhTreeNodeSize, gr_bhTreeNodeSize2, gr_bhLocRecvTreeLevels, &
    gr_bhLocSentTreeLevels, gr_bhTreeNodetype, &
    gr_bhLx, gr_bhLy, gr_bhLz, gr_bhLxHalf, gr_bhLyHalf, gr_bhLzHalf, &
    gr_bhComm, gr_bhTreeNumProcs, gr_bhTreeLnblocks, &
    gr_bhTreeParentTree, gr_bhLocParentTree, gr_bhTreeChild, &
    gr_bhTreeBlocklist, gr_bhTreeLimangle2i, &
    gr_bhGravFac, gr_bhTreeLrefineMax, gr_bhTreeLimangle, &
    gr_bhTreeNFirstLev, gr_bhTreeFirstLevBlocks, &
    gr_bhUseGravity, gr_bhUseTreeRay, GR_TREE_IM, GR_TREE_IX, GR_TREE_IY, & 
    GR_TREE_IZ, GR_TREE_NSIZE, GR_TREE_BNSIZE, GR_TREE_BASE_NSIZE, &
    GR_TREE_BASE_BNSIZE , gr_bhPhysMACComm, gr_bhPhysMACTW, & 
    gr_bhUseUnifiedTW, gr_bhTWMaxQueueSize, gr_bhPriorityQueue, &
    gr_bhTypePQElement, gr_bhTreeMinCellSize, gr_bhTreeMinCellSize2, &
    gr_bhPQNull, gr_bhTWType, GR_TREE_TWSTD, GR_TREE_TWUNI, GR_TREE_TWPQ, &
    gr_bhDensVar, gr_bhGpotVar, gr_bhUseRelAccErr, gr_bhLoadBalancing, &
    gr_bhAcceptAccurateOld, gr_bhMaxBlkWeight
  use tree, ONLY : nodetype, lrefine, child, mchild, maxblocks_tr, mfaces
  use Logfile_interface, ONLY : Logfile_stamp
  use gr_bhLocalInterface, ONLY : gr_bhCalcBlockTreePos
  use gr_bhInterface, ONLY : gr_bhGetTreeNodeSize
 
  implicit none
#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer :: istat
  integer :: i, j, k, ni, nj, nk, hi, hj, hk, maxl
  integer :: nblockx, nblocky, nblockz

  real :: xmax, xmin, ymax, ymin, zmax, zmin
  real :: mcs(MDIM)
  character(len=MAX_STRING_LENGTH) :: strBuff
  integer :: blockID, blockCount
  integer :: blockList(MAXBLOCKS)
  real, POINTER, DIMENSION(:,:,:,:) :: solnData


  gr_bhTreeMyPE = gr_meshMe
  gr_bhTreeNumProcs = gr_meshNumProcs
  gr_bhComm = gr_meshComm
  call RuntimeParameters_get("gr_bhPhysMACTW", gr_bhPhysMACTW)
  call RuntimeParameters_get("gr_bhPhysMACComm", gr_bhPhysMACComm)
  call RuntimeParameters_get("gr_bhTreeLimangle", gr_bhTreeLimangle)
  call RuntimeParameters_get("gr_bhTreeSafeBox", gr_bhTreeSafeBox)
  call RuntimeParameters_get("gr_bhUseUnifiedTW", gr_bhUseUnifiedTW)
  call RuntimeParameters_get("gr_bhTWMaxQueueSize", gr_bhTWMaxQueueSize)
  call RuntimeParameters_get("lrefine_max", gr_bhTreeLrefineMax)
  call RuntimeParameters_get("useGravity", gr_bhUseGravity)
  call RuntimeParameters_get("useTreeRay", gr_bhUseTreeRay)

! highly experimental features, switched off in this version
!  call RuntimeParameters_get("gr_bhLoadBalancing", gr_bhLoadBalancing)
!  call RuntimeParameters_get("gr_bhMaxBlkWeight", gr_bhMaxBlkWeight)
!  call RuntimeParameters_get("gr_bhAcceptAccurateOld", gr_bhAcceptAccurateOld)

  if (gr_bhUseGravity) then
    call RuntimeParameters_get("grv_bhUseRelAccErr", gr_bhUseRelAccErr)
  endif

  gr_bhTreeLimangle2 = gr_bhTreeLimangle**2
  gr_bhTreeLimangle2i = 1.0/gr_bhTreeLimangle2
  gr_bhTreeSafeBoxHalf2 = (0.5*gr_bhTreeSafeBox)**2
  gr_bhDensVar = -1
  gr_bhGpotVar = -1

  ! physical MAC cannot be used for communitaion without being used for Tree Walk
  if (gr_bhPhysMACComm .and. .not. gr_bhPhysMACTW) then
    call Driver_abortFlash ('[gr_bhInit]: gr_bhPhysMACTW must be set true if gr_bhPhysMACComm is true')
  endif

  ! calculate size and structure of the tree node and the tree bottom node
  GR_TREE_NSIZE  = GR_TREE_BASE_NSIZE
  GR_TREE_BNSIZE = GR_TREE_BASE_BNSIZE 
  call Gravity_bhGetNodeStruct(GR_TREE_IM, GR_TREE_IX, GR_TREE_IY, GR_TREE_IZ &
  & , GR_TREE_NSIZE, GR_TREE_BNSIZE)
  call TreeRay_bhGetNodeStruct(GR_TREE_IM, GR_TREE_IX, GR_TREE_IY, GR_TREE_IZ &
  & , GR_TREE_NSIZE, GR_TREE_BNSIZE)
  if (gr_meshMe == MASTER_PE) &
  & print *, "NSIZE = ", GR_TREE_NSIZE, GR_TREE_BNSIZE

  ! Select Tree Walk Type
  gr_bhTWType = GR_TREE_TWSTD
  if (gr_bhUseUnifiedTW) gr_bhTWType = GR_TREE_TWUNI

  ! Check if we support the requested grid geometry.
  if ((NDIM /= 3) .or. (gr_geometry /= CARTESIAN)) then
     call Driver_abortFlash ('[gr_bhInit] ERROR: tree Poisson solver works only in 3D Cartesian geometry')
  endif

  ! check if NBX = NBY = NBZ
  if ((NXB /= NYB) .or. (NXB /= NZB) .or. (NYB /= NZB)) then
     call Driver_abortFlash ('[gr_bhInit] ERROR: must be NXB = NYB = NZB')
  endif
  
  ! set the tree parameters: gr_bhTreeBS, gr_bhTreeLevels
  gr_bhTreeBS = NXB
  gr_bhTreeLevels = 0 ! gr_bhTreeLevels = ln_2 (NBX)
  i = gr_bhTreeBS
  do
    if (mod(i, 2) == 0) then
      i = i / 2
      gr_bhTreeLevels = gr_bhTreeLevels + 1
    else 
      exit 
    endif
    if (i == 1) exit
  enddo
  do i = 0, gr_bhTreeLevels+1
    gr_bhTreeLoff(i) = 1 + GR_TREE_NSIZE*(8**i - 1)/7 ! offset of the level
  enddo

  ! cell size
  allocate(gr_bhTreeCellSize(1:gr_bhTreeLrefineMax,MDIM), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeCellSize in gr_bhInit.F90")
  call Grid_getMinCellSizes(mcs)
  gr_bhTreeCellSize(gr_bhTreeLrefineMax,:) = mcs
  gr_bhTreeMinCellSize  = min(mcs(IAXIS), mcs(JAXIS), mcs(KAXIS))
  gr_bhTreeMinCellSize2 = gr_bhTreeMinCellSize**2
  do i = gr_bhTreeLrefineMax-1,1,-1
    gr_bhTreeCellSize(i,:) = 2*gr_bhTreeCellSize(i+1,:)
  enddo

  ! Sizes of Tree nodes: maximum of a block edge over all three directions
  allocate(gr_bhTreeNodeSize(1:gr_bhTreeLrefineMax+gr_bhTreeLevels), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeNodeSize in gr_bhInit.F90")
  allocate(gr_bhTreeNodeSize2(1:gr_bhTreeLrefineMax+gr_bhTreeLevels), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeNodeSize2 in gr_bhInit.F90")
  do i = 1, gr_bhTreeLrefineMax
    !gr_bhTreeDiag2(i) = (NXB*gr_bhTreeCellSize(i,IAXIS))**2 &
    !            + (NYB*gr_bhTreeCellSize(i,JAXIS))**2 &
    !            + (NZB*gr_bhTreeCellSize(i,KAXIS))**2
    !gr_bhTreeDiag(i) = sqrt(gr_bhTreeDiag2(i))
    gr_bhTreeNodeSize(i) = max(NXB*gr_bhTreeCellSize(i,IAXIS) &
                     & ,   NYB*gr_bhTreeCellSize(i,JAXIS) &
                     & ,   NZB*gr_bhTreeCellSize(i,KAXIS))
    gr_bhTreeNodeSize2(i) = gr_bhTreeNodeSize(i)**2
  enddo
  do i = gr_bhTreeLrefineMax+1, gr_bhTreeLrefineMax+gr_bhTreeLevels
    gr_bhTreeNodeSize2(i) = gr_bhTreeNodeSize2(i-1)/4.0
    gr_bhTreeNodeSize(i) = sqrt(gr_bhTreeNodeSize2(i))
  enddo

#ifdef DEBUG_SOLVERS
  if (gr_bhTreeMyPE == MASTER_PE) then
    print *, "Tree solver initialized"
    do i = 1, gr_bhTreeLrefineMax+gr_bhTreeLevels
      print *, "level,NodeSize: ", i, gr_bhTreeNodeSize(i), gr_bhGetTreeNodeSize(i)
    enddo
  endif
#endif

  ! BlockTreePos and BlockTreeNodeCen
  call gr_bhCalcBlockTreePos()
      
  allocate(gr_bhTreeArray(0:gr_bhTreeNumProcs-1,1:MAXBLOCKS), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeArray in gr_bhInit.F90")

  allocate(gr_bhLocSentTreeLevels(MAXBLOCKS, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhLocSentTreeLevels in gr_bhInit.F90")
  allocate(gr_bhLocRecvTreeLevels(MAXBLOCKS, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhLocRecvTreeLevels in gr_bhInit.F90")
  
  allocate(gr_bhLocCoords(gr_bhTreeBS+1,MDIM, MAXBLOCKS), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhLocCoords in gr_bhInit.F90")
  allocate(gr_bhTreeBCen(MDIM, MAXBLOCKS,0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeBCen in gr_bhInit.F90")

  allocate(gr_bhTreeParentTree(GR_TREE_NSIZE, MAXBLOCKS, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeParentTrees in gr_bhInit.F90")
  allocate(gr_bhLocParentTree(GR_TREE_NSIZE, MAXBLOCKS), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhLocParentTrees in gr_bhInit.F90")
  
  allocate(gr_bhTreeNodetype(maxblocks_tr, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeNodetype in gr_bhInit.F90")
  allocate(gr_bhTreeLrefine(maxblocks_tr, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeLrefine in gr_bhInit.F90")
  allocate(gr_bhTreeChild(2, mchild, maxblocks_tr, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeChild in gr_bhInit.F90")
  allocate(gr_bhTreeLnblocks(0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeLnblocks in gr_bhInit.F90")

  ! NOT USED IN THIS MOMENT, STORES LIST OF BLOCKS ON ALL CPUs
  allocate(gr_bhTreeBlocklist(1:1, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeBlockList in gr_bhInit.F90")

  ! computational domain
  call RuntimeParameters_get("xmin", xmin)
  call RuntimeParameters_get("xmax", xmax)
  call RuntimeParameters_get("ymin", ymin)
  call RuntimeParameters_get("ymax", ymax)
  call RuntimeParameters_get("zmin", zmin)
  call RuntimeParameters_get("zmax", zmax)
  gr_bhLx = xmax - xmin
  gr_bhLy = ymax - ymin
  gr_bhLz = zmax - zmin
  gr_bhLxHalf = 0.5*gr_bhLx
  gr_bhLyHalf = 0.5*gr_bhLy
  gr_bhLzHalf = 0.5*gr_bhLz

  call RuntimeParameters_get("nblockx", nblockx)
  call RuntimeParameters_get("nblocky", nblocky)
  call RuntimeParameters_get("nblockz", nblockz)
  gr_bhTreeNFirstLev = nblockx*nblocky*nblockz
  allocate(gr_bhTreeFirstLevBlocks(2, gr_bhTreeNFirstLev), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeFirstLevBlocks in gr_bhInit.F90")

  ! Priority Queue
  allocate(gr_bhPriorityQueue(1:gr_bhTWMaxQueueSize), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhPriorityQueue in gr_bhInit.F90")
  do i = 1, gr_bhTWMaxQueueSize
    gr_bhPriorityQueue(i) = gr_bhPQNull
  enddo

  ! Load balancing
  ! erase array for recording number of interactions, if it exists
#ifdef WRKL_VAR
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  do blockID = 1, blockCount
    call Grid_getBlkPtr(blockID,solnData,CENTER) 
    solnData(WRKL_VAR,:,:,:) = 0.0
    call Grid_releaseBlkPtr(blockID,solnData,CENTER) 
  enddo
#endif

  ! to ensure all allocations are made before the code continues
  call MPI_Barrier(gr_bhComm, i) 

  if (gr_bhTreeMyPE == MASTER_PE) print *, "Tree solver initialized"
  
  return
end subroutine gr_bhInit


