!!****ih* source/Grid/GridMain/paramesh/Paramesh2/physicaldata
!!
!! NAME
!!
!!    physicaldata
!!
!!
!! SYNOPSIS
!!
!!   include file with physical data 
!!   
!!
!! DESCRIPTION
!!
!!   This file has physical data needed for the code.
!!   Description of data is included in the file close to 
!!   the declarations
!!  
!!***


module physicaldata
  implicit none
#include "constants.h"
#include "Flash.h"

!
!
!
!Pre-Processor Control
!
!#ifndef CRAY
!#define shmem_real_get shmem_get8
!#define SHMEM_REAL_GET SHMEM_GET8
!#define shmem_real_put shmem_put8
!#define SHMEM_REAL_PUT SHMEM_PUT8
!#endif
!
! set pre-processor variable to control use of different timesteps on
! different grid blocks
!#define VAR_DT                                                 !<<< USER EDIT
!
! Does the algorithm use predictor-corrector type timestepping?
!#define PRED_CORR                                              !<<< USER EDIT
!
! Will any grid blocks represent obstacles to flow?
!#define EMPTY_CELLS                                             !<<< USER EDIT
!
! set the model dimension here. If running 2.5D also edit l2p5d below.
! Skip this.  We will set N_DIM directly from compilation command line.
!#define N_DIM 2
!
  integer, parameter :: k1d = 1
  integer, parameter :: ndim = NDIM
!
!-----------------------------------------------------------------
! physicaldata.fh

!!inserted for compatability with paramesh3.0..  i want this always
!!defined for the mesh
      logical, save :: NPGFlag = .false.

!
! set physical dimension of model and number of edges on each grid block
  integer nbedges
!
! l2p5d needs to be declared not matter what the dimension of the
! problem is set to !
  integer,parameter :: l2p5d=0
                                                                ! 1 if 2.5D
                                                                ! 0 if 2 or 3D
!
  parameter(nbedges=ndim*2**(ndim-1))
!
!
!
! an increment variable for the z dimension to enable the same code to
! work for 2D or 3D models.
  integer, parameter::  k3d = K3D
  integer, parameter::  k2d = K2D

!
!
! set size of grid blocks
! The spatial relationship between grid points in parents
! and their children is subtly different depending on whether
! the block dimensions (ie nxb,nyb,nzb) are even or odd. This has 
! significant consequences when defining interpolation within
! restriction or prolongation operations.
!
  integer nxb,nyb,nzb,maxdim

#ifdef NXB
  parameter(nxb = NXB)
#else
  parameter(nxb=8)                                        !<<< USER EDIT  
#endif
  
#ifdef NYB
  parameter(nyb = NYB)
#else
#if NDIM > 1
  parameter(nyb=8)                                    !<<< USER EDIT
#else
  parameter(nyb=1)
#endif
#endif

#ifdef NZB
  parameter(nzb = NZB)
#else
#if NDIM > 2   
  parameter(nzb=8)                                     !<<< USER EDIT
#else
  parameter(nzb=1)
#endif
#endif
  
  
!

!
!
!
  parameter(maxdim=max(nxb,nyb,nzb))
!
!
! these guard cell offsets are required to accomodate differences
! in cases when block dimensions are odd or even
  integer gc_off_x,gc_off_y,gc_off_z
  parameter(gc_off_x=mod(nxb,2))
  parameter(gc_off_y=mod(nyb,2))
  parameter(gc_off_z=mod(nzb,2))
!
! set the maximum number of blocks per processor
#ifdef MAXBLOCKS
  integer,parameter :: maxblocks = MAXBLOCKS
#else
#if NDIM == 3
  parameter (maxblocks = 200)              !<<< USER EDIT
#else /* N_DIM < 3 */
  parameter (maxblocks = 1000)             !<<< USER EDIT 
#endif /*N_DIM*/
#endif /*MAXBLOCKS*/

  integer, parameter :: nvar=NUNK_VARS, nguard=NGUARD
!
!
!
!
!..this include file is very important at this location, as it sets a 
!..parameter (ionmax) that determines the array sizes and do-loop limits 
!..of the mesh, hydro, eos and burn modules. it touches just about everything.
!

!!#include "flash_defines.fh"

!
! set number of unknowns associated with each grid cell
  integer,parameter ::  nvar2=5, nvarsm = 2
!
!
! common block storing the solution for cell-centered quantities.
! unksm stores copies of global variables which DO NOT need guard cells
! AND do not need to be saved from one timestep to the next !!!
!
  integer, parameter ::  il_bnd = GRID_ILO_GC ,iu_bnd = GRID_IHI_GC
  integer, parameter ::  jl_bnd = GRID_JLO_GC ,ju_bnd = GRID_JHI_GC
  integer, parameter ::  kl_bnd = GRID_KLO_GC ,ku_bnd = GRID_KHI_GC

  integer, parameter :: nxlo=NGUARD+1,nylo=NGUARD*K2D+1,nzlo=NGUARD*K3D+1
  integer, parameter :: nxhi=NGUARD+NXB
  integer, parameter :: nyhi=NGUARD*K2D+NYB
  integer, parameter :: nzhi=NGUARD*K3D+NZB
  
  real, DIMENSION(UNK_VARS_BEGIN:UNK_VARS_END,il_bnd:iu_bnd,jl_bnd:ju_bnd,                   &
       &     kl_bnd:ku_bnd,maxblocks), TARGET :: unk
  real, DIMENSION(nvar2,il_bnd:iu_bnd,maxblocks),                     &
       &     TARGET :: unk2
  real, DIMENSION(nvar2,jl_bnd:ju_bnd,maxblocks),                     &
       &     TARGET :: unk3
  real, DIMENSION(nvar2,kl_bnd:ku_bnd,maxblocks),                     &
       &     TARGET :: unk4
  real, DIMENSION(nvarsm,nxlo:nxhi,                                   &
       &     nylo:nyhi,nzlo:nzhi,maxblocks), TARGET :: unksm
  
  
 
!
!
!
! The convention for relating variables associated with cell faces to the
! variables defined at cell centers is as follows:
!
! If iface_off=0 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i-1/2) face of the (i,j,k)-th mesh cell.
! If iface_off=-1 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i+1/2) face of the (i,j,k)-th mesh cell.
!
  integer,parameter ::iface_off=0                                  !<<< USER EDIT
!
!
! The number of data words needed on a cell face is set by nfacevar.
!
!
! 2 added to store strong fields at faces for all components of B
  integer,parameter ::nfacevar=0   !<<< USER EDIT
!
  integer,parameter ::nbndvar=max(1,nfacevar)
!
  integer,parameter ::maxblocksf= 1+(maxblocks-1)*min(1,nfacevar) 
!
! common block storing the solution for cell-face-centered quantities.
  real,target ::                                                          &
       &     facevarx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,              &
       &     kl_bnd:ku_bnd,maxblocksf)                                    &
       &     ,facevary(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+K2D,           &
       &     kl_bnd:ku_bnd,maxblocksf)                                    &
       &     ,facevarz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,               &
       &     kl_bnd:ku_bnd+K3D,maxblocksf)
  
!
! set data length of grid blocks

  integer,parameter ::len_block=iu_bnd*ju_bnd*ku_bnd*NUNK_VARS
!
  integer,parameter ::len_blockfx=(NXB+2*NGUARD+1)*(NYB+2*NGUARD*K2D)*        &
       &                          (NZB+2*NGUARD*K3D)
  integer,parameter ::len_blockfy=(NXB+2*NGUARD)*(NYB+2*NGUARD+1)*            &
       &                          (NZB+2*NGUARD*K3D)
  integer,parameter ::len_blockfz=(NXB+2*NGUARD)*(NYB+2*NGUARD)*              &
       &                          ((NZB+2*NGUARD)*K3D+1)
!
!
! common block for timestep control
  integer,parameter ::maxlevels=20
  common/timecntl/                                                  &
       & time_loc(maxblocks),dtlevel(maxlevels),ldtcomplete(maxblocks)
  real time_loc,dtlevel
  logical ldtcomplete
!
!
#if defined(VAR_DT) || defined(PRED_CORR)
!      parameter(maxblocks_t=(maxblocks-1)*ivar_dt+1)
!      parameter(nvar_t=(NUNK_VARS-1)*ivar_dt+1)
  common/tsolution/                                                 &
       &      t_unk(UNK_VARS_BEGIN:UNK_VARS_END,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,       &
       &       maxblocks)
  real t_unk
#endif
!
#ifdef PRED_CORR
  real ::                                                         &
       &           tfacevarx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,       &
       &                          kl_bnd:ku_bnd,maxblocksf)               &
       &          ,tfacevary(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+K2D,     &
       &                          kl_bnd:ku_bnd,maxblocksf)               &
       &          ,tfacevarz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,         &
       &                          kl_bnd:ku_bnd+K3D,maxblocksf)
  
#endif
!
! To average fluxes set red_f = .25.
! To sum fluxes set red_f = 1.0.
!
! changed -- 2-24-00 
! we are now converting the flux densities to fluxes before the call to
! amr_flux_conserve.  This is to get the proper geometry factors included
! in non-cartesian cases.  The fluxes are converted back after the call.
! 
! red_f is set according to dimension, to sum instead of average
!      parameter(red_f = 0.25)   
!
! NOTE: In the Paramesh2 code included with FLASH (at least FLASH3 and up), the
! parameter red_f is not actually used anywhere. - KW 2012-07-30
#if N_DIM == 1
  real, parameter :: red_f = 0.25
#elif N_DIM == 2
  real, parameter :: red_f = 0.5
#elif N_DIM == 3
  real,parameter :: red_f = 1.0
#endif
!
!-----------------------------------------------------------------
! include header file defining data structure on cell faces
!
!-----------------------------------------------------------------
!
! This file defines a data structure to be used for quantities
! which may need to be defined at grid block interfaces, eg fluxes,
! pressures.
!
!
!
! storage used for fluxes at block boundaries. This is used when conservation
! constraints need to be imposed.
!
! updated 2-15-00 -- allocate storage for the internal energy flux


  integer,parameter ::nfluxvar=NFLUXES

  integer,parameter ::nfluxes=max(1,nfluxvar)

  integer,parameter ::maxblocksfl= 1+(maxblocks-1)*min(1,nfluxvar) 
!
!
!
!..in 1d the flux_y, flux_z, tflux_y, and tflux_z arrays are not used,
!..but do need to be declared. thus, in 1d the parameter maxblocksfl
!..has been replaced with a 1. this significantly reduces the
!..memory footprint for 1d simulations.
!
!..in 2d the flux_z and tflux_z arrays are not used,
!..but do need to be declared. thus, in 2d the parameter maxblocksfl
!..has been replaced with a 1. this significantly reduces the
!..memory footprint for 2d simulations.
!
#if N_DIM == 1
  real, target ::                                                  &
       flux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl),&
       flux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,1),          &
       flux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,1)
!
  real, target ::                                                   &
       tflux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl),&
       tflux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,1),          &
       tflux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,1)

#endif
!
!
#if N_DIM == 2
  real, target ::                                                   &
       flux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl), &
       flux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,maxblocksfl), &
       flux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,1)
  
!
  real, target ::                                                    &
       tflux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl), &
       tflux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,maxblocksfl), &
       tflux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,1)
#endif
!
!
#if N_DIM == 3
  real, target ::                                                  &
       flux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl),&
       flux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,maxblocksfl),&
       flux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,maxblocksfl)
!
  real, target ::                                                    &
       tflux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl), &
       tflux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,maxblocksfl), &
       tflux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,maxblocksfl)
#endif
!
!
!
! storage used for cell edges at block boundaries. 
! This is used when quantities located at cell edge centers need to
! be used consistently at the boundaries between blocks at different
! refinement levels.
!

  integer,parameter ::nedgevar=1         !!! USER EDIT
!
  integer,parameter ::nedges=max(1,nedgevar)

  integer,parameter ::maxblockse= 1+(maxblocks-1)*min(1,nedgevar) 
!
!
!..flash does not presently use these edge storage variables
!..but it might when magnetic fields are included. until then,
!..the maxblockse size declaration has been replaced with 1 in order
!..to reduce the memory footprint.
!
!      common/edges/
!     .  bedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,
!     .    kl_bnd:ku_bnd+1,maxblockse),
!     .  bedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,
!     .    kl_bnd:ku_bnd+1,maxblockse),
!     .  bedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,
!     .    kl_bnd:ku_bnd+1,maxblockse),
!     .  bedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,
!     .    kl_bnd:ku_bnd+1,maxblockse),
!     .  bedge_facez_x(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,
!     .    1:2,maxblockse),
!     .  bedge_facez_y(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,
!     .    1:2,maxblockse),
!
!
  real ::                                                         &
       &  bedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,                       &
       &    kl_bnd:ku_bnd+1,1),                                           &
       &  bedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,                       &
       &    kl_bnd:ku_bnd+1,1),                                           &
       &  bedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,                       &
       &    kl_bnd:ku_bnd+1,1),                                           &
       &  bedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,                       &
       &    kl_bnd:ku_bnd+1,1),                                           &
       &  bedge_facez_x(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,           &
       &    1:2,1),                                                       &
       &  bedge_facez_y(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,           &
       &    1:2,1),                                                       &
       &  recvarx1e(nedges,1:2,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd+1),          &
       &  recvary1e(nedges,il_bnd:iu_bnd+1,1:2,kl_bnd:ku_bnd+1),          &
       &  recvarz1e(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,1:2),          &
       &  recvarx2e(nedges,1:2,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd+1),          &
       &  recvary2e(nedges,il_bnd:iu_bnd+1,1:2,kl_bnd:ku_bnd+1),          &
       &  recvarz2e(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,1:2)
!
!
!..flash does not presently use these edge storage variables
!..but it might when mhd is added in.
!..the maxblockse size declaration has been reduced to 1 in order
!..to help reduce the memory footprint.
!..28nov99 fxt
!
!
!      common/tedges/
!     .  tbedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,
!     .    kl_bnd:ku_bnd+1,maxblockse),
!     .  tbedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,
!     .    kl_bnd:ku_bnd+1,maxblockse),
!     .  tbedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,
!     .    kl_bnd:ku_bnd+1,maxblockse),
!     .  tbedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,
!     .    kl_bnd:ku_bnd+1,maxblockse),
!     .  tbedge_facez_x(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,
!     .    1:2,maxblockse),
!     .  tbedge_facez_y(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,
!     .    1:2,maxblockse)
!
  real ::                                                         &
       &  tbedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,                      &
       &    kl_bnd:ku_bnd+1,1),                                           &
       &  tbedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,                      &
       &    kl_bnd:ku_bnd+1,1),                                           &
       &  tbedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,                      &
       &    kl_bnd:ku_bnd+1,1),                                           &
       &  tbedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,                      &
       &    kl_bnd:ku_bnd+1,1),                                           &
       &  tbedge_facez_x(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,          &
       &    1:2,1),                                                       &
       &  tbedge_facez_y(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,          &
       &    1:2,1)
!
!
!
! workspace arrays used for inter-block communications
  integer,parameter ::nbndmax=max(nbndvar,nfluxes)
  real ::                                                         &
       &     recvarx1(nbndmax,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd)            &
       &    ,recvary1(nbndmax,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd)            &
       &    ,recvarz1(nbndmax,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2)            &
       &    ,bndtempx1(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd)           &
       &    ,bndtempy1(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd)           &
       &    ,bndtempz1(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2)
  
!
!
!
!
! parameters used in communication calls
  integer,parameter ::len_block_bndx=2*(NYB+2*NGUARD*K2D)*                    &
       &                                 (NZB+2*NGUARD*K3D)
  integer,parameter ::len_block_bndy=2*(NXB+2*NGUARD*K2D)*                    &
       &                                 (NZB+2*NGUARD*K3D)
  integer,parameter ::len_block_bndz=2*(NXB+2*NGUARD)*(NYB+2*NGUARD)

  integer,parameter ::len_block_ex=2*(NYB+K2D+2*NGUARD*K2D)*                  &
       &                                 (NZB+K3D+2*NGUARD*K3D)
  integer,parameter ::len_block_ey=2*(NXB+1+2*NGUARD)*                        &
       &                                 (NZB+K3D+2*NGUARD*K3D)
  integer,parameter ::len_block_ez=2*(NXB+1+2*NGUARD)*                        &
       &                                 (NYB+K2D+2*NGUARD)
  

end module physicaldata






