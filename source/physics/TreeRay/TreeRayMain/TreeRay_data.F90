!!****if* source/physics/TreeRay/TreeRayMain/TreeRay_data
!!
!! NAME
!!
!!  TreeRay_data
!!  
!! SYNOPSIS
!!
!!  use TreeRay_data
!!
!! DESCRIPTION
!!
!!  This modules stores the data for the TreeRay unit
!!
!!***

module TreeRay_data

#include "constants.h"

  logical, save :: tr_bhUseTreeRay, updateTreeRay, tr_useGravity
  integer, save :: tr_nSide !healpix resolution
  integer, save :: tr_nPix
  integer, save :: tr_meshMe, tr_comm, tr_numProcs
  integer, save :: tr_bhLRefineMax, tr_bhTreeLevels, tr_bhNTBLevels
  real, save, allocatable :: tr_bhTreeNodeSize(:)
  real, save :: tr_bhMaxDist, tr_bhMinCellSize, tr_mH, tr_mHi, tr_4PIi
  real, save :: tr_nPo4pi 
  real, save :: tr_boltz, tr_lightSpeed
  real, save :: tr_smallx, tr_smlrho

  ! indeces in the tree node
  integer, save :: tr_bhIM, tr_bhIBM, tr_bhIX, tr_bhIY, tr_bhIZ

  ! energy bands
  integer, parameter :: TR_MAXNEB = 1
  integer, save :: tr_nEb = 0
  integer, dimension(TR_MAXNEB) :: tr_mapEbSoln
  integer, save :: tr_iEbEUV = -1

  
  ! rays for all directions for each cell of a block
  real, allocatable :: tr_bhMassRays(:,:,:,:,:), tr_bhVolRays(:,:,:,:,:)
  real, allocatable :: tr_bhSrcfRays(:,:,:,:,:,:), tr_bhEradRays(:,:,:,:,:,:)

  ! column densities in all directions for all cells in a processed block
  ! indeces: quantity, ipix, i, j, k
  integer, save :: tr_nCd = 0
  real, allocatable :: tr_bhCdMaps(:,:,:,:,:)

  ! values along the ray
  integer, save :: tr_bhNR
  integer, parameter :: tr_nFineR = 10
  real, save :: tr_bhRayRadRes
  real,allocatable,dimension(:) :: tr_bhRayR, tr_bhRayR2, &
  & tr_bhRayRi, tr_bhRayR2i !, tr_dr3Int_jmh
  integer,allocatable,dimension(:,:,:) :: tr_bhRadNodeMapInd
  real,allocatable,dimension(:,:,:) :: tr_bhRadNodeMapVal

  integer, save :: tr_boundary(6)  !integer boundary condition

  ! intersection list
  integer, save :: tr_ilNTheta, tr_ilNPhi, tr_ilNNS, tr_ilNR, tr_ilFinePix
  integer, save :: tr_ilNI
  real, save    :: tr_ilNSSampFac = 1.6, tr_ilNSSampFacI  
  real, allocatable :: tr_intersectList(:,:,:,:)


  ! error control
  character(len=MAX_STRING_LENGTH),save :: tr_bhErrControl
  real, save :: tr_bhRelErr, tr_bhMaxRelEradErr, tr_bhLocRelErr
  real, save :: tr_bhOldEradTot = 0.0, tr_bhLocEradTot, tr_bhEradTot = 0.0
  real, save :: tr_bhOldMionTot = 0.0, tr_bhLocMionTot, tr_bhMionTot = 0.0



  logical, save :: tr_debprint

end module TreeRay_data

