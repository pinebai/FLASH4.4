!!****if* source/physics/IncompNS/IncompNSMain/constdens/ins_ab2rk3
!!
!!
!! NAME
!!
!!  ins_ab2rk3
!!
!!
!! SYNOPSIS
!!
!!  ins_ab2rk3(integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount)
!!             real(IN)    :: timeEndAdv
!!             real(IN)    :: dt)
!!
!!
!! DESCRIPTION
!!
!!  Performs a second order Adams Bashforth or third order Runge-
!!  Kutta step on a fractional step time discretization of the
!!  Incompressible Navier Stokes flow problem.
!!
!!  The blockList and blockCount arguments tell this routine on
!!  which blocks and on how many to operate.  blockList is an
!!  integer array of size blockCount that contains the local
!!  block numbers of blocks on which to advance.
!!
!!  dt gives the timestep through which this update should advance.
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - time level at the end of step
!!  dt         - timestep
!!
!!***

subroutine ins_ab2rk3( blockCount, blockList, timeEndAdv, dt)

#include "Flash.h"

  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
                             Grid_getListOfBlocks, &
                             Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_setInterpValsGcell,&
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes,    &
                             Grid_conserveField,     &
                             Grid_updateRefinement,  &
                             Grid_solvePoisson, Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use Grid_interface, ONLY : Grid_conserveField

  use ins_interface, only  :  ins_vt,&
                           ins_rhs3d,&
                           ins_rhs2d,&
                       ins_predictor,&
                      ins_divergence,&
                       ins_corrector,&
                         ins_fluxfix,&
                       ins_fluxfix_p,&
                   ins_computeQinout,&
                   ins_rescaleVelout,&
                   ins_convectVelout,&
                      ins_UstarStats,&
                  ins_pressgradients

  use IncompNS_data, ONLY : ins_isgs, ins_invRe, ins_intschm, ins_prescoeff, ins_meshMe,&
                            ins_restart, ins_nstep, ins_Qin, ins_Qout, ins_predcorrflg, &
                            ins_convvel, ins_alf, ins_gam, ins_rho, ins_gama, ins_alfa, &
                            ins_rhoa, ins_outflowgridChanged, ins_tlevel, &
                            ins_vardt, rkstep, ins_intschm_type
  use IncompNS_data, ONLY : ins_meshComm, ins_domainBC

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Driver_interface, only : Driver_getNStep

  implicit none

#include "constants.h"
#include "IncompNS.h"
  include "Flash_mpi.h"


  !! ---- Argument List ----------------------------------
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
  real,    INTENT(IN) :: timeEndAdv,dt
  !! -----------------------------------------------------

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox


  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,blockID,ii,jj,kk,ierr,i,j,k

  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC) :: newu
  real, dimension(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC) :: newv
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC+1) :: neww

  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_u
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_v
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_w

  integer :: sx,sy,sz,ex,ey,ez

  real dtdxdz,dtdydz,dtdxdy

  integer TA(2),count_rate
  real*8  ET

  integer TAIB(2),count_rateIB
  real*8  ETIB

  real maxfp,minfp,maxflb,minflb

  real bsize(MDIM),coord(MDIM)
  integer datasize(MDIM)

  integer nxc, nyc, nzc
  real del(MDIM)

  integer, dimension(6) :: bc_types
  integer :: idimn,ibound,eachBoundary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real, dimension(2,6)  :: bc_values = 0.
  real poisfact,alfadt

  integer ist,itmx

  logical, save :: firstcall = .true.

  ! debug VAR:
  integer aa,bb,cc
  real :: meanPres,meanVelx,meanVely,meanVelz
  real :: minu,maxu,minv,maxv,minw,maxw,minp,maxp,mndivv,mxdivv
  real :: vecminaux(5),vecmaxaux(5),vecmin(5),vecmax(5)

  logical :: gridChanged

  character(len=6) :: IndNStep
  integer :: NStep
! --------------------------------------------------------------------------


  CALL SYSTEM_CLOCK(TA(1),count_rate)

  newu = 0.
  newv = 0.
  neww = 0.
  flxint_u = 0.
  flxint_v = 0.
  flxint_w = 0.


  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1

!!$  write(*,*) 'NXYZB=',NXB,NYB,NZB,NGUARD
!!$  write(*,*) 'interp_mask_unk=',interp_mask_unk
!!$  write(*,*) 'interp_mask_facex=',interp_mask_facex
!!$  write(*,*) 'interp_mask_facey=',interp_mask_facey
!!$  write(*,*) 'interp_mask_facez=',interp_mask_facez
!!$  write(*,*) 'interp_mask_unk_res=',interp_mask_unk_res
!!$  write(*,*) 'interp_mask_facex_res=',interp_mask_facex_res
!!$  write(*,*) 'interp_mask_facey_res=',interp_mask_facey_res
!!$  write(*,*) 'interp_mask_facez_res=',interp_mask_facez_res


  do idimn = 1,NDIM
  do ibound = LOW, HIGH
     eachBoundary = 2*(idimn-1)+ibound
     select case (ins_domainBC(ibound,idimn))
     case (PERIODIC)
#ifdef FLASH_GRID_UG
        bc_types(eachBoundary) = PERIODIC
#else
        bc_types(eachBoundary) = GRID_PDE_BND_PERIODIC !MG_BND_PERIODIC
#endif
     case (SLIP_INS,NOSLIP_INS,INFLOW_INS,NEUMANN_INS,MOVLID_INS,OUTFLOW_INS)
#ifdef FLASH_GRID_UG
        bc_types(eachBoundary) = OUTFLOW
#else
        bc_types(eachBoundary) = GRID_PDE_BND_NEUMANN !MG_BND_NEUMANN
#endif
     case default
     if (ins_meshMe .eq. MASTER_PE) then
        write(*,*) 'ins_ab2rk3 Error: Boundary Conditions match for Poisson Solver not defined.'
        write(*,*) 'ins_ab2rk3 Error: LOW-HIGH,AXIS=',ibound,idimn
        write(*,*) 'ins_ab2rk3 Error: ins_domainBC(ibound,idimn) =',ins_domainBC(ibound,idimn)
     endif
     call Driver_abortFlash('ins_ab2rk3 Error: BCs do not have matching Poisson solver BCs')
     end select
  enddo
  enddo

  ! shift timesteps
  do i = -rkstep,-1
     ins_vardt(i) = ins_vardt(i+1)
  end do
  ins_vardt(0) = dt


  ! Select Euler step (for starting) of Adams-Bashforth coefficients
  ! 2nd order Adams Bashforth coefficients (for constant timestep only):
  if (ins_intschm .eq. AB2_SCHM) then
     ins_gam(1) = 1.5
     ins_gam(2) = 0.0
     ins_gam(3) = 0.0
     ins_rho(1) = -0.5
     ins_rho(2) =  0.0
     ins_rho(3) =  0.0

     itmx = 1

  ! 2nd Order Adams-Bashforth coefficents for variable timesteps
  elseif (ins_intschm .eq. AB2_SCHM_V ) then
     ins_gam(1) = 1+ins_vardt(0)/(2.*ins_vardt(-1))
     ins_gam(2) = 0.0
     ins_gam(3) = 0.0
     ins_rho(1) = -ins_vardt(0)/(2.*ins_vardt(-1))
     ins_rho(2) =  0.0
     ins_rho(3) =  0.0

     itmx = 1

  ! 3rd order Runge Kutta coefficients
  elseif (ins_intschm .eq. RK3_SCHM) then
     ins_gam(1) = 8./15.
     ins_gam(2) = 5./12.
     ins_gam(3) = 3./4.
     ins_rho(1) = 0.0
     ins_rho(2) = -17./60.
     ins_rho(3) = -5./12.

     itmx = 3
  else
     if (ins_meshMe .eq. MASTER_PE) then
        write(*,*) 'Unknown Incompressible Flow integrator scheme:'
        write(*,*) 'ins_schm=',ins_intschm
        write(*,*) 'where ins_schm=2 Adams-Bashforth, ins_schm=3 Runge-Kutta'
     endif
     itmx = 0                   ! assign something to avoid compiler warnings
  endif

  ! Euler coefficients (starting from scratch for Adams-Bashforth):
  if ((ins_nstep .eq. 1) .and. (ins_restart .eqv. .false.) .and. &
      (ins_intschm_type .eq. INS_INTSCHM_MULTISTEP)) then
     ins_gam(1) = 1.0; ins_gam(2) = 0.0; ins_gam(3) = 0.0
     ins_rho(1) = 0.0; ins_rho(2) = 0.0; ins_rho(3) = 0.0
     itmx = 1
  endif

  ins_alf = ins_gam + ins_rho

  ! Set Interpolation values for guardcell-filling:
  call Grid_setInterpValsGcell(.true.)

!!$  ! Fill Guardcells for all variables if grid has changed:
!!$  if (ins_outflowgridChanged) then
!!$     gcMask = .TRUE.
!!$     call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
!!$       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
!!$     ins_outflowgridChanged = .false.
!!$  endif

  ! Fill Guardcells with interpolation values defind for ins:
  if(ins_restart .and. firstcall) then
     ! FILL GUARDCELLS FOR RESTART VELOCITIES AND PRESSURE:
     ! ---- ---------- --- ------- ---------- --- --------
     gcMask = .FALSE.
     gcMask(PRES_VAR) = .TRUE.                                ! pressure
     gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! u
     gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! v
#if NDIM == 3
     gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! w
#endif
     ins_predcorrflg = .false.
     call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
          maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

     firstcall=.false.
  endif


  ins_tlevel = timeEndAdv - dt

  ! Timestep Loop:
  do ist = 1,itmx

  ins_gama = ins_gam(ist)
  ins_rhoa = ins_rho(ist)
  ins_alfa = ins_alf(ist)

  ins_tlevel = ins_tlevel + ins_alfa*dt

  ! These two subroutine calls ar used in case of outflow BCs, only when NEUMANN_INS and
  ! OUTFLOW_INS are present.
  ! Compute inflow volume ratio: (Not computed on NOT_BOUNDARY, NEUMANN_INS, OUTFLOW_INS)
  call ins_computeQinout( blockCount, blockList, .true., ins_Qin)

  ! For OUTFLOW_INS condition compute convective velocity
  call ins_convectVelout( blockCount, blockList, ins_convvel)
  !if(ins_meshMe .eq. MASTER_PE) write(*,*) 'After convect',ins_convvel(HIGH,:)

  ! TURBULENT VISCOSITY COMPUTATION:
  ! --------- --------- -----------
#if NDIM == 3
  if (ins_isgs .NE. 0) then
     do lb = 1,blockCount
        blockID = blockList(lb)

        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)

        ! Get blocks coord and bsize
        ! Bounding box:
        call Grid_getBlkBoundBox(blockId,boundBox)
        bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)

        call Grid_getBlkCenterCoords(blockId,coord)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        ! calculate turbulent viscosity
        call ins_vt(ins_isgs,NGUARD,nxc,nyc,nzc,                   &
                    ins_invRe,del(DIR_X),del(DIR_Y),del(DIR_Z),    &
                    coord,bsize,                                   &
                    facexData,&
                    faceyData,&
                    facezData,&
                    solnData)

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

     enddo
     ! apply BC and fill guardcells for turbulent viscosity
     gcMask = .FALSE.
     gcMask(TVIS_VAR) = .TRUE.                            ! only turbulent viscosity
     call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
  endif
#endif

!!$  CALL SYSTEM_CLOCK(TA(1),count_rate)


  ! Compute forcing pressure gradients if required:
  call ins_pressgradients(ins_tlevel,ins_alfa*dt)

#if NDIM < 3
  ! just so we can pass a valid argument in the 2D case to subroutines which will then ignore it anyway:
  allocate(facezData(max(VELC_FACE_VAR,RHDS_FACE_VAR),1,1,1))
#endif

  call Timers_start("RightHandSide_Predictor")
  ! COMPUTE RIGHT HAND SIDE AND PREDICTOR STEP:
  ! ------- ----- ---- ---- --- --------- ----
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     ! compute RHS of momentum equation
     call ins_rhs3d (  facexData(VELC_FACE_VAR,:,:,:),            &
                       faceyData(VELC_FACE_VAR,:,:,:),            &
                       facezData(VELC_FACE_VAR,:,:,:),            &
                       solnData(TVIS_VAR,:,:,:),                  &
                       ins_invRe,                                 &
                       blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                       blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                       blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                       del(DIR_X),del(DIR_Y),del(DIR_Z),newu,newv,neww )
     call Grid_getBlkPtr(blockID,facezData,FACEZ) ! not sure why again... -KW

#elif NDIM ==2
     ! compute RHS of momentum equation
     call ins_rhs2d(  facexData(VELC_FACE_VAR,:,:,:),            &
                      faceyData(VELC_FACE_VAR,:,:,:),            &
                      ins_invRe,                                 &
                      blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                      blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                      del(DIR_X),del(DIR_Y),newu,newv)

#endif
     ! compute intermediate velocities
     call ins_predictor(facexData(VELC_FACE_VAR,:,:,:),&
                        faceyData(VELC_FACE_VAR,:,:,:),&
                        facezData(VELC_FACE_VAR,:,:,:),&
                        newu,newv,neww,                &
                        facexData(RHDS_FACE_VAR,:,:,:),&
                        faceyData(RHDS_FACE_VAR,:,:,:),&
                        facezData(RHDS_FACE_VAR,:,:,:),&
                        solnData(PRES_VAR,:,:,:),      &
                        dt,del(DIR_X),del(DIR_Y),del(DIR_Z),      &
            blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
            blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
            blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
            ins_gama,ins_rhoa,ins_alfa )

     ! save RHS for next step
     facexData(RHDS_FACE_VAR,:,:,:) = newu(:,:,:)
     faceyData(RHDS_FACE_VAR,:,:,:) = newv(:,:,:)


     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


#if NDIM ==3
     facezData(RHDS_FACE_VAR,:,:,:) = neww(:,:,:)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif


  enddo

  call Timers_stop("RightHandSide_Predictor")

!!$   !CALL SYSTEM_CLOCK(TA(2),count_rate)
!!$   !ET=REAL(TA(2)-TA(1))/count_rate
!!$   !write(*,*) 'Predictor time =',ET


  call Timers_start("Gcell_IntermVelocs")
  ! APPLY BC AND FILL GUARDCELLS FOR INTERMEDIATE VELOCITIES:
  ! ----- -- --- ---- ---------- --- ------------ ----------
  gcMask = .FALSE.
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! ustar
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! vstar
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! wstar
#endif
  ins_predcorrflg = .true.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

  ! FIX FLUXES FOR USTAR: (Only for AMR grids)
  ! --- ------ --- -----
#ifdef FLASH_GRID_PARAMESH
  ! Fix fluxes at block boundaries
  call ins_fluxfix(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                   blockCount,blockList)
#endif
  call Timers_stop("Gcell_IntermVelocs")

!!$  ! Compute outflow mass volume ratio: (computed on NEUMANN_INS, OUTFLOW_INS)
!!$  call ins_computeQinout( blockCount, blockList, .false., ins_Qout)
!!$  !write(*,*) 'Qout before ref=',ins_Qout
!!$
!!$  ! Rescale Velocities at outflows for overall conservation:
!!$  call ins_rescaleVelout(  blockCount, blockList, ins_Qin, ins_Qout)

  if (ist .eq. 1) then
  call Timers_start("Grid_updateRefinement")
  call Grid_updateRefinement( ins_nstep, timeEndAdv ,gridChanged)
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  call Timers_stop("Grid_updateRefinement")
  endif

  ! Compute inflow volume ratio: (Not computed on NOT_BOUNDARY, NEUMANN_INS, OUTFLOW_INS)
  call ins_computeQinout( blockCount, blockList, .true., ins_Qin)  ! Shizhao

  ! Compute DivUstar - delta_mass and print to screen, add to ins_Qin:
  call ins_UstarStats( blockCount, blockList, .true., .true.)

  ! Compute outflow mass volume ratio: (computed on NEUMANN_INS, OUTFLOW_INS)
  call ins_computeQinout( blockCount, blockList, .false., ins_Qout)
  !write(*,*) 'Qout after ref=',ins_Qout

  ! Rescale Velocities at outflows for overall conservation:
  call ins_rescaleVelout(  blockCount, blockList, ins_Qin, ins_Qout)

  ! DIVERGENCE OF USTAR:
  ! ---------- -- -----
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM ==3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif


     ! compute divergence of intermediate velocities
     call ins_divergence(facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                       del(DIR_X),del(DIR_Y),del(DIR_Z),&
                       solnData(DUST_VAR,:,:,:) )

     ! Poisson RHS source vector
     solnData(DUST_VAR,                                   &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) =   &
     solnData(DUST_VAR,                                   &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))/(dt*ins_alfa)

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo

!  call Grid_computeVarMean(DUST_VAR,meanPres)
!  if (ins_meshMe .eq. MASTER_PE) write(*,*) 'Mean Div Ustar A=',meanPres*(dt*ins_alfa)

  call Timers_start("Grid_solvePoisson")
  ! SOLUTION OF POISSON EQUATION FOR PRESSURE:
  ! -------- -- ------- -------- --- --------
  poisfact = 1.0
  call Grid_solvePoisson (DELP_VAR, DUST_VAR, bc_types, bc_values, poisfact)
  call Timers_stop("Grid_solvePoisson")

!  call Grid_computeVarMean(PRES_VAR,meanPres)
!  if (ins_meshMe .eq. MASTER_PE) write(*,*) 'Mean Pressure=',meanPres
!  call Grid_computeVarMean(DELP_VAR,meanPres)
!  if (ins_meshMe .eq. MASTER_PE) write(*,*) 'Mean DeltaP=',meanPres

  call Timers_start("Gcells_DelP")
   ! FIX FLUXES FOR dDELP/dxi :
#ifdef FLASH_GRID_UG
  ! Don't Fix Fluxes in block Boundaries
  ! Fill Guardcells for DelP: Used in boundary dDelp/dx fluxes:
  gcMask = .FALSE.
  gcMask(DELP_VAR) = .TRUE.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,  &
      maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask, &
      selectBlockType=ACTIVE_BLKS)
 ! ---------------------------------------------------------------------
#else
  ! fix fluxes at block boundaries
  ! fix dp gradient fluxes at block boundaries
  call ins_fluxfix_p(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                     DELP_VAR,blockCount,blockList)
#endif
  call Timers_stop("Gcells_DelP")


  call Timers_start("Corrector")
  ! CORRECTOR STEP:
  ! --------- ---
  alfadt = ins_alfa*dt
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take depth to be 1:
     del(DIR_Z) = 1.
#endif

     dtdydz = alfadt*(del(DIR_Y)*del(DIR_Z))**(-1.)   !(dy*dz)**-1.
     dtdxdz = alfadt*(del(DIR_X)*del(DIR_Z))**(-1.)   !(dx*dz)**-1.
     dtdxdy = alfadt*(del(DIR_X)*del(DIR_Y))**(-1.)   !(dx*dy)**-1.

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(DIR_X)-NGUARD
     ey = dataSize(DIR_Y)-NGUARD*K2D
     ez = dataSize(DIR_Z)-NGUARD*K3D

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM ==3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

#ifdef FLASH_GRID_UG
     ! UNIFORM Grid:
     ! west face
     facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez) =                  &
          facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez) -             &
          alfadt/del(DIR_X)*(solnData(DELP_VAR,sx,sy:ey,sz:ez)- &
                             solnData(DELP_VAR,sx-1,sy:ey,sz:ez))

     ! east face
     facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez) =                  &
          facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez) -             &
          alfadt/del(DIR_X)*(solnData(DELP_VAR,ex+1,sy:ey,sz:ez)- &
                             solnData(DELP_VAR,ex,sy:ey,sz:ez))

     ! south face
     faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez) =                  &
          faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez) -             &
          alfadt/del(DIR_Y)*(solnData(DELP_VAR,sx:ex,sy,sz:ez)- &
                             solnData(DELP_VAR,sx:ex,sy-1,sz:ez))

     ! north face
     faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez) =                  &
          faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez) -             &
          alfadt/del(DIR_Y)*(solnData(DELP_VAR,sx:ex,ey+1,sz:ez)- &
                             solnData(DELP_VAR,sx:ex,ey,sz:ez))


#if NDIM == 3
     ! front face
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz) =                  &
          facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz) -             &
          alfadt/del(DIR_Z)*(solnData(DELP_VAR,sx:ex,sy:ey,sz)- &
                             solnData(DELP_VAR,sx:ex,sy:ey,sz-1))

     ! back face
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1) =                  &
          facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1) -             &
          alfadt/del(DIR_Z)*(solnData(DELP_VAR,sx:ex,sy:ey,ez+1)- &
                             solnData(DELP_VAR,sx:ex,sy:ey,ez))
#endif

#else

     ! AMR GRID:
     ! update block boundary velocities using corrected fluxes
     ! X direction:
     call Grid_getFluxData(blockID, IAXIS, &
                           flxint_u, dataSize)


     ! west face
     facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez) =       &
          facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez) -  &
          dtdydz*flxint_u(VELC_FLUX,sx,sy:ey,sz:ez)

     ! east face
     facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez) =       &
          facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez) -  &
          dtdydz*flxint_u(VELC_FLUX,ex+1,sy:ey,sz:ez)


     ! Y direction:
     call Grid_getFluxData(blockID, JAXIS, &
                           flxint_v, dataSize)

     ! south face
     faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez) =       &
          faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez) -  &
          dtdxdz*flxint_v(VELC_FLUX,sx:ex,sy,sz:ez)

     ! north face
     faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez) =       &
          faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez) -  &
          dtdxdz*flxint_v(VELC_FLUX,sx:ex,ey+1,sz:ez)

#if NDIM == 3

     ! Z direction:
     call Grid_getFluxData(blockID, KAXIS, &
                           flxint_w, dataSize)

     ! front face
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz) =       &
          facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz) -  &
          dtdxdy*flxint_w(VELC_FLUX,sx:ex,sy:ey,sz)

     ! back face
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1) =       &
          facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1) -  &
          dtdxdy*flxint_w(VELC_FLUX,sx:ex,sy:ey,ez+1)

#endif

#endif


     ! update divergence-free velocities (not on block boundary)
     call ins_corrector( facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
                         solnData(DELP_VAR,:,:,:),&
                         sx,ex,sy,ey,sz,ez,&
                         dt,del(DIR_X),del(DIR_Y),del(DIR_Z),ins_alfa)

     ! update pressure
     solnData(PRES_VAR,:,:,:) = ins_prescoeff*solnData(PRES_VAR,:,:,:) + &
                                solnData(DELP_VAR,:,:,:)

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

  enddo ! End of corrector loop
#if NDIM < 3
  deallocate(facezData)
#endif
  call Timers_stop("Corrector")

  call Timers_start("Gcell_FinalVelocsP")
  ! FILL GUARDCELLS FOR FINAL VELOCITIES AND PRESSURE:
  ! ---- ---------- --- ----- ---------- --- --------
  ! The pressure fill is used to compute distributed forces on
  ! immersed bodies.
  gcMask = .FALSE.
  gcMask(PRES_VAR) = .TRUE.                                ! pressure
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! u
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! v
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! w
#endif
  ins_predcorrflg = .false.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
  call Timers_stop("Gcell_FinalVelocsP")

!!$  ! FIX FLUXES FOR UFINAL: (Only for AMR grids, needed when using
!!$  ! boundary restriction of order > 1, and force_consistency_at_srl_interfaces=True)
!!$  ! --- ------ --- -----
!!$#ifdef FLASH_GRID_PARAMESH
!!$  if (force_consistency) then
!!$  ! Fix fluxes at block boundaries
!!$  call ins_fluxfix(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
!!$                   blockCount,blockList)
!!$  endif
!!$#endif

  enddo ! End  of time substeps loop

#ifdef CHK_BND_VLC
    ! Output interpolated velocity
    call Driver_getNStep(NStep)
    if( mod(NStep,CHK_BND_VLC) == 0) then
    write(indNstep, '(I6.6)') NStep/CHK_BND_VLC
    call gr_intepVel('parts-uv-step-'//indNStep//'-final',0)
    endif
#endif


  ! Restore Interpolation values for guardcell-filling:
  call Grid_setInterpValsGcell(.false.)

  ! ------------------------------------------------------------------------------------------
  ! Check min max divergence:
  ! ----- --- --- ----------
  mxdivv = -10.**(10.)
  mndivv =  10.**(10.)
  maxu   = mxdivv; maxv = maxu; maxw = maxu; maxp = maxu;
  minu   = mndivv; minv = minu; minw = minu; minp = minu;
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

  mxdivv = max( mxdivv,maxval( (facexData(VELC_FACE_VAR,NGUARD+2:nxc,NGUARD+1:nyc-1,NGUARD+1:nzc-1) - &
                    facexData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_X) + &
                   (faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+2:nyc,NGUARD+1:nzc-1) - &
                    faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_Y) + &
                   (facezData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+2:nzc) - &
                    facezData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_Z) ))

  mndivv = min( mndivv,minval( (facexData(VELC_FACE_VAR,NGUARD+2:nxc,NGUARD+1:nyc-1,NGUARD+1:nzc-1) - &
                    facexData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_X) + &
                   (faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+2:nyc,NGUARD+1:nzc-1) - &
                    faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_Y) + &
                   (facezData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+2:nzc) - &
                    facezData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_Z) ))


  maxu = max(maxu,maxval(facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))
  minu = min(minu,minval(facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))

  maxv = max(maxv,maxval(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)))
  minv = min(minv,minval(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)))

  maxw = max(maxw,maxval(facezData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1)))
  minw = min(minw,minval(facezData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1)))

  maxp = max(maxp,maxval(solnData(PRES_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))
  minp = min(minp,minval(solnData(PRES_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))


#elif NDIM == 2

  maxu = max(maxu,maxval(facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))
  minu = min(minu,minval(facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))

  maxv = max(maxv,maxval(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)))
  minv = min(minv,minval(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)))

  maxw = 0.0
  minw = 0.0

  maxp = max(maxp,maxval(solnData(PRES_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))
  minp = min(minp,minval(solnData(PRES_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))

  mxdivv = max( mxdivv,maxval( (facexData(VELC_FACE_VAR,NGUARD+2:nxc,NGUARD+1:nyc-1,1) - &
                    facexData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,1))/del(DIR_X) + &
                   (faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+2:nyc,1) - &
                    faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,1))/del(DIR_Y) ))

  mndivv = min( mndivv,minval( (facexData(VELC_FACE_VAR,NGUARD+2:nxc,NGUARD+1:nyc-1,1) - &
                    facexData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,1))/del(DIR_X) + &
                   (faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+2:nyc,1) - &
                    faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,1))/del(DIR_Y) ))


#endif

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo

  vecmaxaux(1) = mxdivv
  vecminaux(1) = mndivv
  vecmaxaux(2) = maxu
  vecminaux(2) = minu
  vecmaxaux(3) = maxv
  vecminaux(3) = minv
  vecmaxaux(4) = maxw
  vecminaux(4) = minw
  vecmaxaux(5) = maxp
  vecminaux(5) = minp

  call MPI_Allreduce(vecmaxaux, vecmax, 5, FLASH_REAL,&
                     MPI_MAX, ins_meshComm, ierr)

  call MPI_Allreduce(vecminaux, vecmin, 5, FLASH_REAL,&
                     MPI_MIN, ins_meshComm, ierr)


  if (ins_meshMe .eq. MASTER_PE) then
     write(*,*) ' '
     write(*,'(A24,2g14.6)') ' Min , Max  U =',vecmin(2),vecmax(2) !minu,maxu
     write(*,'(A24,2g14.6)') ' Min , Max  V =',vecmin(3),vecmax(3) !minv,maxv
#if NDIM == 3
     write(*,'(A24,2g14.6)') ' Min , Max  W =',vecmin(4),vecmax(4) !minw,maxw
#endif
     write(*,'(A24,2g14.6)') ' Min , Max  P =',vecmin(5),vecmax(5) !minp,maxp
     write(*,'(A24,2g14.6)') ' Min , Max  Divergence =',vecmin(1),vecmax(1) !mndivv,mxdivv
  endif

  !----------------------------------------------------------------------------------------------------



  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1))/count_rate
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total AB Step Time =',ET


END SUBROUTINE ins_ab2rk3
