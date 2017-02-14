!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_init
!!
!! NAME
!!
!!  TreeRay_init
!!  
!! SYNOPSIS
!!
!!  TreeRay_init()
!!
!! DESCRIPTION
!!
!!  Initialize unit scope variables in the TreeRay unit, which are typically the 
!!  runtime parameters.  This routine must be called once by Driver_initFlash.F90. 
!!  Calling multiple times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!! 
!!  tr_bhUseTreeRay  BOOLEAN true  Controls turning on/off the compiled TreeRay unit
!!
!! NOTES
!!   
!!  Each implementation of TreeRay has its own runtime parameters.  Be sure to check
!!  the documentation or Config files to see them.
!!
!!***

subroutine TreeRay_init ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype, &
    Driver_getComm, Driver_getNumProcs
  use Grid_interface, ONLY : Grid_getMinCellSize
  use TreeRay_data, ONLY : tr_useGravity, tr_bhUseTreeRay, &
    tr_nside, tr_bhMaxDist, tr_ilNTheta, tr_ilNPhi, tr_ilNNS, tr_ilNR, tr_ilFinePix, tr_ilNI, &
    tr_boundary, tr_bhRayR, tr_bhRayR2, &
    tr_bhRayRi, tr_bhRayR2i, tr_bhMassRays, tr_mHi, &
    tr_bhVolRays, tr_bhSrcfRays, tr_bhEradRays, tr_nPix, tr_bhNR, tr_smallx, &
    tr_mH, tr_lightSpeed, tr_4PIi, tr_boltz, tr_bhRayRadRes, tr_bhMinCellSize, &
    tr_bhNTBLevels, tr_bhLRefineMax, tr_bhTreeLevels, tr_bhRelErr, &
    tr_smlrho, tr_bhErrControl, tr_meshMe, tr_comm, tr_numProcs, tr_nEb, tr_nCd, &
    tr_bhTreeNodeSize, tr_bhCdMaps, tr_nPo4pi 

  use tr_bhLocalInterface, ONLY : tr_bhGenIntersectList
  use gr_bhInterface, ONLY : gr_bhGetTreeNodeSize
  use tr_odInterface, ONLY : tr_odInit
  !use tr_osInterface, ONLY : tr_osInit
  !use tr_rpInterface, ONLY : tr_rpInit

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
   

!  logical, save :: tr_bhUseTreeRay
  integer :: i, ir, ir2, istat
  real :: dx, max_ray_length, theta, phi
  real :: xmin, xmax, ymin, ymax, zmin, zmax

  character(len=MAX_STRING_LENGTH) :: tr_boundary_type

!==============================================================================

  call Driver_getMype(MESH_COMM, tr_meshMe)
  call Driver_getComm(GLOBAL_COMM, tr_comm)
  call Driver_getNumProcs(GLOBAL_COMM, tr_numProcs)

  !! It is a failure to invoke the stub when tr_bhUseTreeRay is set TRUE.
 
  call RuntimeParameters_get ("useTreeRay", tr_bhUseTreeRay)
  call RuntimeParameters_get ("useGravity", tr_useGravity)
  call RuntimeParameters_get ("lrefine_max", tr_bhLRefineMax)
  call RuntimeParameters_get ("tr_nSide", tr_nSide)
  call RuntimeParameters_get ("tr_bhMaxDist", tr_bhMaxDist)

  call RuntimeParameters_get ("tr_bhRayRadRes", tr_bhRayRadRes)
  call RuntimeParameters_get ("tr_bhErrControl", tr_bhErrControl)
  call RuntimeParameters_get ("tr_bhRelErr", tr_bhRelErr)

  call RuntimeParameters_get ("tr_ilNTheta", tr_ilNTheta)
  call RuntimeParameters_get ("tr_ilNPhi", tr_ilNPhi)
  call RuntimeParameters_get ("tr_ilNNS", tr_ilNNS)
  call RuntimeParameters_get ("tr_ilNR", tr_ilNR)
  call RuntimeParameters_get ("tr_ilFinePix", tr_ilFinePix)
  !call RuntimeParameters_get ("tr_ilNI", tr_ilNI)
  call RuntimeParameters_get ("tr_boundary_type", tr_boundary_type)

  call RuntimeParameters_get ("xmin", xmin)
  call RuntimeParameters_get ("xmax", xmax)
  call RuntimeParameters_get ("ymin", ymin)
  call RuntimeParameters_get ("ymax", ymax)
  call RuntimeParameters_get ("zmin", zmin)
  call RuntimeParameters_get ("zmax", zmax)

  call RuntimeParameters_get ("smallx", tr_smallx)
  call RuntimeParameters_get ("smlrho", tr_smlrho)
  call PhysicalConstants_get( 'proton mass', tr_mH)
  call PhysicalConstants_get( 'Boltzmann', tr_boltz)
  call PhysicalConstants_get( 'speed of light', tr_lightSpeed)

  if (.not. tr_bhUseTreeRay) return

  tr_mHi = 1.0 / tr_mH
  tr_4PIi = 1.0/(4*PI)
  ! total number of levels: amr tree + block trees
  i = NXB
  tr_bhTreeLevels = 0
  do
    i = i / 2
    tr_bhTreeLevels = tr_bhTreeLevels + 1
    if (i == 1) exit
  enddo

  tr_bhNTBLevels = tr_bhLRefineMax + tr_bhTreeLevels
  allocate(tr_bhTreeNodeSize(tr_bhNTBLevels), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate tr_bhTreeNodeSize")
  do i = 1, tr_bhNTBLevels
    tr_bhTreeNodeSize(i) = gr_bhGetTreeNodeSize(i)
    if (tr_meshMe .eq. MASTER_PE) &
    & print *, "TreeRay: level,NodeSize: ", i, tr_bhTreeNodeSize(i)
  enddo

  !print*,'TreeRay_init: tr_bhUseTreeRay, nside, maxdist=',tr_bhUseTreeRay, tr_nSide, tr_maxdist

  ! check dimensionality, set boundary conditions
  if (NDIM .lt. 3) call Driver_abortFlash("Currently TreeRay ONLY works in 3D!")
  tr_nPix = tr_nSide * tr_nSide * 12
  tr_ilNI = tr_nPix
  tr_nPo4pi = tr_nPix / (4.0*PI)
  select case (tr_boundary_type)
    case("isolated")
      tr_boundary(1:6) = ISOLATED
    case("periodic")
      tr_boundary(1:6) = PERIODIC
  end select

  ! generate intersection list
  call tr_bhGenIntersectList()

  ! call inits of sub-modules (needed to calculate tr_bhNEB and tr_nCd)
  call tr_odInit()
  !call tr_osInit()
  !call tr_rpInit()

#ifdef TR_OPTICALDEPTH
  ! allocate column density maps
  allocate(tr_bhCdMaps(tr_nCd, 0:tr_nPix-1, 1:NXB, 1:NYB, 1:NZB), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate tr_bhCdMaps")
#endif


end subroutine TreeRay_init
