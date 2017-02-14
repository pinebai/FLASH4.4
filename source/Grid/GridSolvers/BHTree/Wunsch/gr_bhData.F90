!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhData
!!
!! NAME
!!
!!  gr_bhData
!!
!!
!! SYNOPSIS
!!
!!  use gr_bhData
!!
!!
!! DESCRIPTION
!!
!!  Variable declarations for the tree Poisson solver.
!!
!!   
!!***

module gr_bhData

  use Grid_interface, ONLY: GRID_PDE_BND_ISOLATED, GRID_PDE_BND_PERIODIC

  implicit none
#include "constants.h"

  logical, save              :: gr_bhPhysMACTW, gr_bhPhysMACComm
  logical, save              :: gr_bhPhysMACTW_step, gr_bhPhysMACComm_step
  logical, save              :: gr_bhUseUnifiedTW, gr_bhUseRelAccErr
  logical, save              :: gr_bhLoadBalancing, gr_bhAcceptAccurateOld
  logical, save              :: gr_bhFirstCall = .true.
  integer, save              :: GR_TREE_NSIZE, GR_TREE_BNSIZE
  integer, save              :: gr_bhDensVar, gr_bhGpotVar
  integer, save              :: gr_bhTreeBS, gr_bhTreeLevels, gr_bhTreeZones
  integer, save              :: gr_bhTreeLoff(0:255)
  integer, save              :: gr_bhTreeMaxlnblocks
  integer, save              :: gr_bhTreeMyPE, gr_bhTreeNumProcs
  integer, save              :: gr_bhComm !replaces MPI_COMM_WORLD
  integer, save              :: gr_bhTreeLrefineMax
  integer, save              :: gr_bhTWType
  integer, save              :: gr_bhTreeNewTW, gr_bhTreeOldAccept
  real, save                 :: gr_bhTreeDcount(1:3)
  real, save                 :: gr_bhGravFac
  real, save                 :: gr_bhTreeLimangle, gr_bhTreeLimangle2, gr_bhTreeLimangle2i
  real, save                 :: gr_bhTreeSafeBox, gr_bhTreeSafeBoxHalf2
  real, save                 :: gr_bhLx, gr_bhLy, gr_bhLz, gr_bhLxHalf, gr_bhLyHalf, gr_bhLzHalf
  logical, save              :: gr_bhUseGravity = .false.
  logical, save              :: gr_bhUseTreeRay = .false.
  real, save                 :: gr_bhOAAvg, gr_bhOAMin, gr_bhOAMax, gr_bhOACnt
  real, save :: gr_bhTotMass

  ! number and indeces of physical modules
  integer, parameter :: GR_TREE_NMODULES = 2
  integer, parameter :: GR_TREE_IGRAVITY = 1
  integer, parameter :: GR_TREE_ITREERAY = 1
  
  ! size of the tree node needed by the tree solver, shared by other physical units
  ! recently, mass and coordinates of the mass centre
  integer, parameter :: GR_TREE_BASE_NSIZE = 4
  ! size of the tree bottom node needed by the tree solver, shared by other physical units
  ! recently, only mass
  integer, parameter :: GR_TREE_BASE_BNSIZE = 1

  ! positions of basic variables in tree node array
  integer, parameter :: GR_TREE_IM = 1
  integer, parameter :: GR_TREE_IX = 2
  integer, parameter :: GR_TREE_IY = 3
  integer, parameter :: GR_TREE_IZ = 4

  ! Tree Walk Types
  integer, parameter :: GR_TREE_TWSTD = 1
  integer, parameter :: GR_TREE_TWUNI = 2
  integer, parameter :: GR_TREE_TWPQ  = 3

  ! Supported boundary constants.
  integer, parameter :: GR_TREE_BND_ISOLATED = GRID_PDE_BND_ISOLATED !!  = 0
  integer, parameter :: GR_TREE_BND_PERIODIC = GRID_PDE_BND_PERIODIC !!  = 1

  ! Remember here the boundary condition for which Grid_solvePoisson is called.
  integer, save :: gr_bhBndType(6)

  ! Indeces in gr_bhBlockTreePos
  integer, parameter :: GR_TREE_BTP_N   = 13
  integer, parameter :: GR_TREE_BTP_POS = 1
  integer, parameter :: GR_TREE_BTP_LEV = 2
  integer, parameter :: GR_TREE_BTP_I   = 3
  integer, parameter :: GR_TREE_BTP_J   = 4
  integer, parameter :: GR_TREE_BTP_K   = 5
  integer, parameter :: GR_TREE_BTP_C1  = 6
  integer, parameter :: GR_TREE_BTP_C8  = 13

  ! array of pointers to trees 2D = (CPUs, blocks)
  type p_tree
    real, dimension(:), pointer :: p
  end type p_tree
  type(p_tree), save, dimension(:,:), allocatable :: gr_bhTreeArray

  ! array of pointers to messages 1D = (CPUs)
  type p_message
    real, dimension(:), pointer :: p
  end type p_message

  ! levels of trees sent to different CPUs (block, toCPU)
  integer, save, allocatable :: gr_bhLocSentTreeLevels(:,:)
  ! levels of trees sent to different CPUs (block, fromCPU)
  integer, save, allocatable :: gr_bhLocRecvTreeLevels(:,:)

  ! Cell sizes (lrefine,MDIM)
  real, save, allocatable :: gr_bhTreeCellSize(:,:) ! is a function of lrefine and MDIM only
  real, save :: gr_bhTreeMinCellSize, gr_bhTreeMinCellSize2
  ! Block sizes (longest edge) (lrefine)
  real, save, allocatable :: gr_bhTreeNodeSize(:), gr_bhTreeNodeSize2(:)
  ! Positions (indeces), levels and children of nodes in the block-tree array
  integer, save, allocatable :: gr_bhBlockTreePos(:,:)
  ! Geometrical centres of nodes in the block-tree array
  ! - with respect to the block centre, distances meassured in grid cells
  real, save, allocatable :: gr_bhBlockTreeNodeCen(:,:)

  ! Positions of all block centres (dim, block, cpu)
  real, save, allocatable :: gr_bhTreeBCen(:,:,:)

  ! Coordinates of all local blocks (gr_bhTreeBS+1, dim, block)
  real, save, allocatable :: gr_bhLocCoords(:,:,:)
  
  ! mass and mass centre positions of parent blocks (mass + MC position, block, cpu)
  ! filled also for LEAF blocks (in BuildTreeBlock) to communicate date for cases parent and child are at different CPU
  real, save, allocatable :: gr_bhTreeParentTree(:,:,:)
  real, save, allocatable :: gr_bhLocParentTree(:,:)

  ! amr tree structure - needed for parent blocks
  ! crated in treeComBlkProperties
  integer, save, allocatable :: gr_bhTreeNodetype(:,:)
  integer, save, allocatable :: gr_bhTreeLrefine(:,:)
  integer, save, allocatable :: gr_bhTreeChild(:,:,:,:)
  integer, save, allocatable :: gr_bhTreeNeigh(:,:,:,:)
  integer, save, allocatable :: gr_bhTreeLnblocks(:)

  ! stores indeces of leaf and parent blocks in a continous array (block_number, cpu)
  ! created in treeComBlkProperties
  integer, save, allocatable :: gr_bhTreeBlocklist(:, :)
  
  ! blocks on the top of the tree/forest
  integer, save, allocatable :: gr_bhTreeFirstLevBlocks(:, :)
  integer, save :: gr_bhTreeNFirstLev

  ! Priority queue for TreeWalk
  TYPE :: gr_bhTypePQElement
    integer :: tr, cpu, btp, int_type
    real :: dr(MDIM+2)
    real :: perr ! maximum partial error
  END TYPE gr_bhTypePQElement
  integer, save :: gr_bhTWMaxQueueSize, gr_bhPQSize
  TYPE(gr_bhTypePQElement), save, dimension(:), allocatable :: gr_bhPriorityQueue
  TYPE(gr_bhTypePQElement), parameter :: gr_bhPQNull = gr_bhTypePQElement(-1, -1, -1, -1, -1.0, -1.0)
  real, save :: gr_bhPQSumSquare
  integer, parameter :: GR_TREE_INT_SELF = 1
  integer, parameter :: GR_TREE_INT_CC   = 2
  integer, parameter :: GR_TREE_INT_CN   = 3
  integer, parameter :: GR_TREE_INT_CP   = 4

  ! iteration of the solver (needed by TreeRay)
  integer, parameter :: gr_bhMaxSolverIter = 10

  ! Checking accuracy of blocks
  integer, parameter :: gr_bhNTestCells = 8
  integer, dimension(MDIM,gr_bhNTestCells), parameter :: gr_bhTestCells = reshape(source = &
    & (/(/1,1,1/),  (/1,1,NZB/),   (/1,NYB,1/),   (/1,NYB,NZB/), &
    & (/NXB,1,1/), (/NXB,1,NZB/), (/NXB,NYB,1/), (/NXB,NYB,NZB/)/) &
    & , shape = (/MDIM,gr_bhNTestCells/))

  ! Load balancing
  real, save :: gr_bhWrklMin, gr_bhWrklMinGlob = 0.0
  real, save :: gr_bhWrklMax, gr_bhWrklMaxGlob = 1.0
  ! decay "time"-scale of the block workload for 
  ! Ornstein-Uhlenbeck scheme; in number of iteractions
  real, parameter :: gr_bhWrklDecay = 10. , gr_bhWrklDecayI = 1.0/gr_bhWrklDecay
  real, save :: gr_bhMaxBlkWeight
  

end module gr_bhData


