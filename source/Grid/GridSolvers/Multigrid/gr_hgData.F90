!!****if* source/Grid/GridSolvers/Multigrid/gr_hgData
!!
!! NAME
!!  gr_hgData
!!
!! SYNOPSIS
!!
!!  use gr_hgData
!!
!! DESCRIPTION
!!
!!  Data module for the Huang Greengard Multigrid solver
!!
!!***

Module gr_hgData

  use tree, ONLY: maxblocks_tr, nchild

#include "constants.h"  ! for MDIM
#include "Flash.h" ! for PFFT_WITH_MULTIGRID
  ! Data structures for prolongation

  integer, allocatable, dimension(:), save       :: send_prolong_req
  real, allocatable, dimension(:,:,:,:,:), save  :: send_prolong_data
  real, allocatable, dimension(:,:,:),save       :: recv_prolong_data

  integer, save                :: nmax1, nmax2, hg_myPE, nbbuf_prolong
  integer, dimension(8),save   :: n1off, n2off, n3off

  real,allocatable,dimension(:,:,:), save :: Px
  real,allocatable,dimension(:,:,:), save :: Py
  real,allocatable,dimension(:,:,:), save :: Pz
  real, dimension(-2:2,2), save  :: Pns, Pew, Pud  
  
  !Data structures for restriction

  integer, save                :: hg_restrict_n1, hg_restrict_n2, hg_restrict_n3
  integer, save                ::  nbbuf_restrict

  integer, save                :: nxl1, nxl2, nyl1, nyl2, nzl1, nzl2
  integer, save                :: nxr1, nxr2, nyr1, nyr2, nzr1, nzr2

  real, allocatable, save      :: send_restrict_data(:,:,:,:)
  real, allocatable, save      :: recv_restrict_data(:,:,:)
  integer, allocatable, save      :: send_restrict_req(:)  

!***************************************************
  
  integer,save :: gr_hgMeshRefineMax, gr_hgMeshRefineMin
  
  ! Average value of source term (subtracted off when doing periodic boundaries).
  
  real,save    :: gr_hgAvgSource

  ! A place to save nodetypes as PARAMESH has them before we change
  ! them to perform various single-level operations.
  
  integer, allocatable,dimension(:), save :: gr_hgSaveNodetype
  ! A place to save child data before the prolongation messes with it.
  
  logical, allocatable,dimension(:), save :: gr_hgSaveNewchild


  ! Ranges of interior indices for blocks.
  
  integer,save :: hg_ili, hg_iui, hg_jli, hg_jui, hg_kli, hg_kui
  
  ! Ranges of exterior indices for blocks.
  
  integer,save :: hg_ile, hg_iue, hg_jle, hg_jue, hg_kle, hg_kue
  
  ! Boundary conditions.
  
  integer, save, dimension(2*MDIM) :: gr_hgBndTypes

! Supported geometry constants.

  
  integer,save :: gr_hgGeometry

  ! Just doing a mg_quadrant?
  
  logical,save :: gr_hgQuadrant

  integer, save :: gr_hgSolnIndex  

  real, save   :: hg_cx(-2:2,NXB), hg_cy(-2:2,NYB), hg_cz(-2:2,NZB)

  integer, save       :: gr_hgMaxCorrections
  real, save          :: gr_hgMaxResidualNorm
  logical, save       :: gr_hgPrintNorm

  logical, save :: gr_hgCurrentGcReq_extrap
  logical, save :: gr_hgNowActive = .FALSE.

#ifdef PFFT_WITH_MULTIGRID
  logical, save :: gr_hgUsingPfft = .true.
#else
  logical, save :: gr_hgUsingPfft = .false.
#endif

end Module gr_hgData
