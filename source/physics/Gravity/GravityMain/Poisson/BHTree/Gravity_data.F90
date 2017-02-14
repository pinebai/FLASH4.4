!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_data
!!
!! NAME
!!
!!  Gravity_data
!!  
!! SYNOPSIS
!!
!!  use Gravity_data
!!
!! DESCRIPTION
!!
!!  This modules stores the data for the Gravity unit
!!
!!***

module Gravity_data

#include "constants.h"
#include "Flash.h"

  ! Common Gravity module variables
  character(len=MAX_STRING_LENGTH), save :: grav_boundary_type !string boundary condition
  integer, save :: grav_boundary(6)  !integer boundary condition

  integer, save :: grav_geometry  !mesh geometry
  integer, save :: grv_meshMe, grv_meshNumProcs, grv_meshComm
  integer, save :: grv_commSize=1

  logical, save :: useGravity, updateGravity
  logical, save :: grav_temporal_extrp !extrapolate or otherwise rescale
  logical, save :: grav_unjunkPden

  real,    save :: grav_poisfact

  ! Variables specific to the BHTree
  character(len=MAX_STRING_LENGTH), save :: grv_bhMAC

  ! configuration parameters
  logical, save :: grv_bhPhysMACTW, grv_bhPhysMACComm
  logical, save :: grv_bhUseRelAccErr
  logical, save :: grv_bhUseEwald



  ! integer vars (mainly degree of the multipole expansion for APE MAC)
  integer, save :: grv_bhMPDegree, grv_bhMACNum
  integer, save :: grv_bhMPDegree_p2, grv_bhMPDegree_half, grv_bhMPDegree_halfp1

  ! real vars (mainly error in acceleration)
  real,    save :: grv_bhMeanBlockAccErrInv
  real,    save :: grv_bhAccErr, grv_bhAccErrInv
  real,    save :: grv_poisson_max = 0.0, grv_sink_max = 0.0

  ! integer constants (mainly indeces to obtain given quantities from arrays)
  integer, save :: grv_bhNSIZE = 0, grv_bhBNSIZE = 0
  integer, save :: grv_bhNTBLevels, grv_bhTreeLevels, grv_bhLRefineMax
  integer, save :: grv_bhIM, grv_bhIBM, grv_bhIX, grv_bhIY, grv_bhIZ
  integer, save :: grv_bhNODE5, grv_bhIDMIN, grv_bhIB2, grv_bhIB3
  integer, save :: grv_bhA1DIR, grv_bhA2DIR, grv_bhA3DIR
  integer, save :: grv_defaultGpotVar = GPOT_VAR

  ! real constants
  real,    save :: grv_bhPiGhalf
  real,    save :: grv_bhNewton
  real,    save :: grv_bhSContrGeom ! self contribution geometry

  ! variables related to periodic boundaries and the Ewald field
!#ifdef GRAV_TREE_EWALD_V42
! old version of computing the Ewald field
  logical, save :: grv_bhGenEwaldAccV42, grv_bhGenEwaldPotV42
  integer, save :: grv_bhEwaldFieldNxV42, grv_bhEwaldFieldNyV42, grv_bhEwaldFieldNzV42
  integer, save :: grv_bhEwaldNRefV42
  integer, save :: grv_bhDirectionQuadV42
  real,    save :: grv_bhEwaldLMaxV42, grv_bhDyIV42, grv_bhDzIV42, grv_bhMinEFSizeV42
  real, save, allocatable :: grv_bhTreeEwaldAccV42(:,:,:,:,:)
  real, save, allocatable :: grv_bhTreeEwaldPotV42(:,:,:,:)
!#endif

! variables needed for both old and new implementations
  integer, save :: grv_bhEwaldSeriesN
  real,    save :: grv_bhDxI, grv_bhLx,  grv_bhLy,  grv_bhLz

! new computation of the Ewald field
  integer, save :: grv_bhEwald_periodicity, grv_margin
  real, save, allocatable :: grv_bhTreeEwald(:,:,:,:)
  real, save :: grv_bhPotConst, grv_bhLayerAccConst, grv_bhLayerPotConst
  real, save :: grv_bhDLog, grv_bhDLogI
  real, save, allocatable :: grv_bhLogfield(:), grv_bhLogfieldDer(:)
  real, save :: grv_L2inv

  ! External potential
  logical, save :: grv_useExternalPotential, grv_usePoissonPotential
  character(len=MAX_STRING_LENGTH), save :: grv_bhExtrnPotType, grv_bhExtrnPotFile
  integer, parameter :: grv_bhExtrnPotNMax = 50000
  integer, save :: grv_bhExtrnPotN, grv_bhExtrnPotIType
  real, save    :: grv_bhExtrnPotCenterX, grv_bhExtrnPotCenterY, grv_bhExtrnPotCenterZ
  real, save    :: grv_bhExtrnPotDel
  real, save, allocatable :: grv_bhExtrnPotCoord(:),  grv_bhExtrnPotPot(:)
  real, save, allocatable :: grv_bhExtrnPotAcc(:)
  integer, parameter :: grv_bhEPTypeX = 1
  integer, parameter :: grv_bhEPTypeY = 2
  integer, parameter :: grv_bhEPTypeZ = 3
  integer, parameter :: grv_bhEPTypeR = 4

  ! constants representing individual MACs
  integer, parameter :: grv_bhMAC_APE = 1
  integer, parameter :: grv_bhMAC_MPE = 2
  integer, parameter :: grv_bhMAC_SS  = 3

  ! meaning of the fifth field in tree nodes
  integer, parameter :: grv_bhN5_NONE = 0
  integer, parameter :: grv_bhN5_DMIN = 1 ! minimum distance
  integer, parameter :: grv_bhN5_B2   = 2 ! second order moment

  ! B3 moment of a uniform unit cube with rho = 1
  ! obtained numerically, because I didn't find the analytical integral :(
  !real, parameter :: GRAV_B3_CONST = 0.137674463674

end module Gravity_data
