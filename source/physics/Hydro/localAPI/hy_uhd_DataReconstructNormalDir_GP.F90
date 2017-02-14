!!****if* source/physics/Hydro/localAPI/hy_uhd_DataReconstructNormalDir_GP
!!
!! NAME
!!
!!  hy_uhd_DataReconstructNormalDir_GP
!!
!! SYNOPSIS
!!
!! call hy_uhd_DataReconstructNormalDir_GP(integer(IN) :: dir,
!!                                    real(IN)    :: dt,
!!                                    real(IN)    :: delta,
!!                                    pointer(IN) :: DataMultiD(:,:),
!!                                    pointer(IN) :: DataGravMultiD(:),
!!                                    real(IN)    :: x1,x2,x3,
!!                                    real(IN)    :: FlatCoeff,
!!                                    logical(IN) :: TransUpdateOnly,
!!                                    real(OUT)   :: lambda0(:),
!!                                    real(OUT)   :: leig0(:,:),
!!                                    real(OUT)   :: reig0(:,:),
!!                                    real(OUT)   :: Wp(:),
!!                                    real(OUT)   :: Wm(:),
!!                                    real(OUT)   :: sig(:),
!!                                    real(IN),optional :: dnBnormal,
!!                                    real(IN),optional :: dnGLMnormal,
!!                                    real(IN),optional :: aBnormal(:),
!!                                    real(IN),optional :: aGLMnormal(:),
!!                                    real(OUT),optional :: Sr(:),
!!                                    real(OUT),optional :: Sl(:),
!!                                    real(OUT),optional :: SpcSig(:),
!!
!! ARGUMENTS
!!
!!  dir         - normal direction along which the reconstuction is performed
!!  dt          - timestep
!!  delta       - deltas in each {x,y,z} direction
!!  DataMultiD  - pointer array holding neighboring stencil hydro/MHD data for reconstruction
!!  DataGravMultiD- pointer array holding neighboring stencil gravity data for reconstruction
!!  x1,x2,x3    - grid points
!!  FlatCoeff   - flattening parameters, primarily for PPM
!!  TransUpdateOnly - a switch for a selective transverse flux update in the normal direction
!!  lambda      - eigenvalue
!!  leig        - left eigenvector
!!  reig        - right eigenvector
!!  Wp,Wm       - left(minus) and right(plus) Riemann states
!!  sig         - transverse flux term
!!  dnBnormal   - MHD multidimensional term for constrained-transport (i.e., USM) MHD scheme
!!  dnGLMnormal - additional term for GLM-MHD
!!  aBnormal    - MHD multidimensional term for constrained-transport (i.e., USM) MHD scheme
!!  aGLMnormal  - additional term for GLM-MHD
!!  Sr          - right Riemann state of species and mass scalars 
!!  Sl          - left  Riemann state of species and mass scalars 
!!  SpcSig      - transverse flux term for spieces and mass scalars
!!
!!
!! DESCRIPTION
!!
!!  Dummy for a routine that used GP for reconstruction.
!!
!! NOTES
!!
!!  This is a dummy implementation. It does not do anything.
!!
!!***

Subroutine hy_uhd_DataReconstructNormalDir_GP&
     (dir,dt,delta,DataMultiD,DataGravMultiD,&
      x1,x2,x3,radius, &
      FlatCoeff,TransUpdateOnly,&
      lambda0,leig0,reig0,&
      Wp,Wm,sig,&
      dnBnormal,aBnormal,&     !These are optional
      dnGLMnormal,aGLMnormal,& !These are optional
      Sr,Sl,SpcSig)  

  implicit none

#include "Flash.h"
#include "UHD.h"

 
  !!-----Arguments---------------------------------------------------------
  integer,intent(IN) :: dir
  real,   intent(IN) :: dt,delta
#if NDIM < 3
  real, pointer, dimension(:,:,:)   :: DataMultiD,DataGravMultiD
#elif NDIM == 3
  real, pointer, dimension(:,:,:,:) :: DataMultiD,DataGravMultiD
#endif
  real, pointer, dimension(:) :: x1
  real, pointer, dimension(:) :: x2
  real, pointer, dimension(:) :: x3
  integer, intent(IN) :: radius
  real,    intent(IN) :: FlatCoeff
  logical, intent(IN) :: TransUpdateOnly
  real, dimension(HY_WAVENUM),intent(OUT) :: lambda0
  real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
  real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
  real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
  real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
  !optional arguments for MHD-----
  real,intent(IN), optional :: dnBnormal,dnGLMnormal
  real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
  real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
  real,intent(OUT),dimension(:),          optional :: SpcSig
  !!------------------------------------------------------------------------

  lambda0(:) = 0.0
  leig0(:,:) = 0.0
  reig0(:,:) = 0.0
  Wp(:) = 0.0
  Wm(:) = 0.0
  sig(:) = 0.0

  if (present(Sr)) Sr(:) = 0.0
  if (present(Sl)) Sl(:) = 0.0
  if (present(SpcSig)) SpcSig(:) = 0.0
End Subroutine Hy_uhd_DataReconstructNormalDir_GP
