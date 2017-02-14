!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!
!! 
!! SYNOPSIS
!!
!!  call Gravity_init()
!!
!!
!! DESCRIPTION
!!
!!  Initialize the Barnes-Hut tree Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the tree_common module
!!
!!***

subroutine Gravity_init()

  use Gravity_data, ONLY : grv_bhNewton, grav_boundary, grav_geometry, &
       grav_poisfact, grv_meshMe, grv_meshNumProcs, grv_meshComm, &
       updateGravity, useGravity, grv_bhPiGhalf, grv_bhSContrGeom, &
       grv_bhMeanBlockAccErrInv, grv_bhMAC, grv_bhMACNum, &
       grv_bhMPDegree, grv_bhMPDegree_p2, grv_bhMPDegree_half, &
       grv_bhMPDegree_halfp1, grv_bhMAC_APE, grv_bhMAC_MPE, grv_bhMAC_SS, &
       grv_bhUseRelAccErr, grv_bhAccErr, grv_bhIDMIN, &
       grv_bhN5_B2, grv_bhN5_DMIN, grv_bhN5_NONE, grv_bhIB2, & 
       grv_bhIB3, grv_bhUseEwald, grv_bhA1Dir, grv_bhA2Dir, grv_bhA3Dir, &
       grv_bhPhysMACTW, grv_bhPhysMACComm, grv_bhUseRelAccErr, &
       grav_unjunkPden, grv_bhNTBLevels, grv_bhLRefineMax, grv_bhTreeLevels, &
       grv_useExternalPotential, grv_usePoissonPotential, grv_bhExtrnPotType, &
       grv_bhExtrnPotFile, grv_bhExtrnPotCenterX, grv_bhExtrnPotCenterY, &
       grv_bhExtrnPotCenterZ
  use grv_bhInterface, ONLY : grv_elintF, grv_readExtrnPotential, &
       grv_bhGenEwaldField
  use Grid_interface, ONLY : Grid_getMinCellSizes
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype,&
       Driver_getComm, Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  character(len=MAX_STRING_LENGTH) :: strGeometry
  character(len=MAX_STRING_LENGTH) :: grav_boundary_type, grav_boundary_type_x,&
     grav_boundary_type_y, grav_boundary_type_z
  integer :: maxl, i, istat
  real :: a1, a2, a3, hk, theta
  real :: mcs(MDIM)

  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getComm(MESH_COMM,grv_meshComm)
  call Driver_getNumProcs(MESH_COMM,grv_meshNumProcs)

  ! read some general runtime parameters
  call RuntimeParameters_get("geometry", strGeometry)
  call RuntimeParameters_mapStrToInt(strGeometry, grav_geometry)
  call RuntimeParameters_get("useGravity", useGravity)
  call RuntimeParameters_get("updateGravity", updateGravity)
  call RuntimeParameters_get("grav_unjunkPden", grav_unjunkPden)
  call RuntimeParameters_get ("lrefine_max", grv_bhLRefineMax)


  ! the next two parameters belong to GridSolvers/BHTree module
  call RuntimeParameters_get("gr_bhPhysMACTW", grv_bhPhysMACTW)
  call RuntimeParameters_get("gr_bhPhysMACComm", grv_bhPhysMACComm)

  ! read gravity module runtime parameters
  call RuntimeParameters_get("grv_bhMAC", grv_bhMAC)
  call RuntimeParameters_get("grv_bhMPDegree", grv_bhMPDegree)
  call RuntimeParameters_get("grv_bhUseRelAccErr", grv_bhUseRelAccErr)
  call RuntimeParameters_get("grv_bhAccErr", grv_bhAccErr)

  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
  call RuntimeParameters_get("grav_boundary_type_x", grav_boundary_type_x)
  call RuntimeParameters_get("grav_boundary_type_y", grav_boundary_type_y)
  call RuntimeParameters_get("grav_boundary_type_z", grav_boundary_type_z)

  call RuntimeParameters_get("grv_useExternalPotential", grv_useExternalPotential)
  call RuntimeParameters_get("grv_usePoissonPotential",  grv_usePoissonPotential)
  call RuntimeParameters_get("grv_bhExtrnPotType",       grv_bhExtrnPotType)
  call RuntimeParameters_get("grv_bhExtrnPotFile",       grv_bhExtrnPotFile)
  call RuntimeParameters_get("grv_bhExtrnPotCenterX",    grv_bhExtrnPotCenterX)
  call RuntimeParameters_get("grv_bhExtrnPotCenterY",    grv_bhExtrnPotCenterY)
  call RuntimeParameters_get("grv_bhExtrnPotCenterZ",    grv_bhExtrnPotCenterZ)
  
  call RuntimeParameters_get("grv_bhNewton", grv_bhNewton)
  if (grv_bhNewton <= 0.0) then
    call PhysicalConstants_get("Newton", grv_bhNewton)
  endif

  if (.not. useGravity) return

  grav_poisfact  = 4.0 * PI * grv_bhNewton
  grv_bhPiGhalf  = 0.5 * PI * grv_bhNewton
  ! total number of levels: amr tree + block trees
  i = NXB
  grv_bhTreeLevels = 0
  do
    i = i / 2
    grv_bhTreeLevels = grv_bhTreeLevels + 1
    if (i == 1) exit
  enddo
  grv_bhNTBLevels = grv_bhLRefineMax + grv_bhTreeLevels



  ! physical MAC is used for communication cannot be used 
  ! together with Relative Acceleration Error control
  if (grv_bhPhysMACComm .and. grv_bhUseRelAccErr) then
     call Driver_abortFlash("Gravity_init:"//&
          " physical MAC used for communication cannot be used together with Relative Acceleration Error control")
  endif

  select case (grv_bhMAC)
    case("ApproxPartialErr")
      grv_bhMACNum = grv_bhMAC_APE
    case("MaxPartialErr")
      grv_bhMACNum = grv_bhMAC_MPE
      if (grv_bhUseRelAccErr) then
        call Driver_abortFlash("Gravity_init: MaxPartialErr MAC cannot be used with Relative Acceleration Error control")
      endif
    case("SumSquare")
      call Driver_abortFlash ('[Gravity_init] ERROR: the SumSquare MAC is not implemented in this version')
      !grv_bhMACNum = grv_bhMAC_SS
      !gr_bhTWType = GR_TREE_TWPQ
      !if (grv_bhUseRelAccErr) then
      !  call Driver_abortFlash("Gravity_init: SumSquare MAC cannot be used with Relative Acceleration Error control")
      !endif
    case default
      call Driver_abortFlash("Gravity_init: unrecognized or unsupported MAC criterion")
  end select
  
  ! MPDegree used for ApproxPartialErr criterion
  grv_bhMPDegree_p2     = grv_bhMPDegree + 2
  grv_bhMPDegree_half   = grv_bhMPDegree / 2
  grv_bhMPDegree_halfp1 = grv_bhMPDegree_half + 1

  ! grid cell sizes sorted: a1 > a2 > a3
  call Grid_getMinCellSizes(mcs)
  if (mcs(IAXIS) > mcs(JAXIS)) then
    if (mcs(JAXIS) > mcs(KAXIS)) then
      grv_bhA1Dir = IAXIS
      grv_bhA2Dir = JAXIS
      grv_bhA3Dir = KAXIS
    else 
      if (mcs(IAXIS) > mcs(KAXIS)) then
        grv_bhA1Dir = IAXIS
        grv_bhA2Dir = KAXIS
        grv_bhA3Dir = JAXIS
      else
        grv_bhA1Dir = KAXIS
        grv_bhA2Dir = IAXIS
        grv_bhA3Dir = JAXIS
      endif
    endif
  else 
    if (mcs(IAXIS) > mcs(KAXIS)) then
      grv_bhA1Dir = JAXIS
      grv_bhA2Dir = IAXIS
      grv_bhA3Dir = KAXIS
    else 
      if (mcs(KAXIS) > mcs(JAXIS)) then
        grv_bhA1Dir = KAXIS
        grv_bhA2Dir = JAXIS
        grv_bhA3Dir = IAXIS
      else
        grv_bhA1Dir = JAXIS
        grv_bhA2Dir = IAXIS
        grv_bhA3Dir = KAXIS
      endif
    endif
  endif

  ! calculate eliptic integral used for self-contribution of a cell
  ! cell is approximated by a uniform tri-axial ellipsoid
  ! contribution to the potential is the potential in the centre of that ellipsoid
  a1 = mcs(grv_bhA1Dir)
  a2 = mcs(grv_bhA2Dir)
  a3 = mcs(grv_bhA3Dir)

  if ((a3/a1 - 1.0) < 1d-99) then
    grv_bhSContrGeom = 1.0 ! sphere
  else
    hk = sqrt((a1*a1-a2*a2)/(a1*a1-a3*a3+1d-99))
    theta = acos(a3/a1)
    grv_bhSContrGeom = grv_elintF(theta, hk)/(sin(theta)+1d-99)
  endif


  !! FUTURE : need to add scale factor for cosmology
  !grav_scaleFactor = 1

  select case (grav_boundary_type)
    case("isolated")
      grav_boundary = ISOLATED
      grv_bhUseEwald = .false.
    case("periodic")
      grav_boundary = PERIODIC
      grv_bhUseEwald = .true.
    case("mixed")
      grv_bhUseEwald = .true.
      select case (grav_boundary_type_x)
        case("isolated")
          grav_boundary(1) = ISOLATED
          grav_boundary(2) = ISOLATED
        case("periodic")
          grav_boundary(1) = PERIODIC
          grav_boundary(2) = PERIODIC
        case default
          call Driver_abortFlash("Gravity_init: unrecognized or unsupported gravity boundary type in x direction")
      end select
      select case (grav_boundary_type_y)
        case("isolated")
          grav_boundary(3) = ISOLATED
          grav_boundary(4) = ISOLATED
        case("periodic")
          grav_boundary(3) = PERIODIC
          grav_boundary(4) = PERIODIC
        case default
          call Driver_abortFlash("Gravity_init: unrecognized or unsupported gravity boundary type in y direction")
      end select
      select case (grav_boundary_type_z)
        case("isolated")
          grav_boundary(5) = ISOLATED
          grav_boundary(6) = ISOLATED
        case("periodic")
          grav_boundary(5) = PERIODIC
          grav_boundary(6) = PERIODIC
        case default
          call Driver_abortFlash("Gravity_init: unrecognized or unsupported gravity boundary type in z direction")
      end select
    case default
      call Driver_abortFlash("Gravity_init: unrecognized or unsupported gravity boundary type")
  end select

  ! PERIODIC BOUNDARIES
  call grv_bhGenEwaldField()
  

  ! EXTERNAL POTENTIAL
  call grv_readExtrnPotential()

  return
end subroutine Gravity_init
