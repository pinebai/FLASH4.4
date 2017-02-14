!!****if* source/physics/RadTrans/RadTransMain/MGD/Unified/RadTrans
!!
!!  NAME 
!!
!!  RadTrans
!!
!!  SYNOPSIS
!!
!!  call RadTrans( integer(IN) :: nblk,
!!                 integer(IN) :: blklst(nblk),
!!                 real(IN)    :: dt, 
!!       optional, integer(IN) :: pass)
!!
!!  DESCRIPTION 
!!      This subroutine performs the radiatiative transfer calculation
!!      for this step using multigroup diffusion.
!!
!!      In this variant, a system of equations in multiple coupled scalar
!!      unknown functions ("variables") is passed to a HYPRE grid solver
!!      to solve as one big unified linear matrix problem.
!!      In this variant, the set of unknown variables consists of one
!!      energy density per radiation energy group, plus one electron
!!      energy density variable.
!!      In this variant, emission and absorption coupling terms
!!      couple electrons to each of the radiation groups; the radiation
!!      energy groups are coupled to each other only indirectly via
!!      the electron energy.
!!
!!      In this variant, electron heat conduction is not implemented,
!!      although in principle it would be possible.
!!
!!      This variant does not represent ione energy as a variable, thus
!!      cannot account for ion-electron thermal coupling or ion head conduction.
!!
!!      This variant cannot be used with radiation energy groups that
!!      are split across multiple instances of the mesh ('mesh replication'),
!!      all energy groups must exist in UNK at the same time.
!!
!! ARGUMENTS
!!
!!   nblk   : The number of blocks in the list
!!   blklst : The list of blocks on which the solution must be updated
!!   dt     : The time step
!!   pass   : reverses solve direction
!!
!!***
#define DEBUG_GRID_GCMASK
!!#define USE_FLUXLIMITER_FCB
#define USE_SOLVESCALAR_FCB
subroutine RadTrans(nblk, blklst, dt, pass)
  use rt_data, ONLY: rt_useMGD, rt_mgdDomainBC, rt_mgdBounds, rt_mgdFlMode, &
       rt_mgdFlCoef, rt_mgdNumGroups, rt_mgdBcVals, rt_mgdthetaImplct, &
       rt_mgdthetaC, rt_mgdthetaD, rt_tightIonCoupling, &
       rt_timeGroups, rt_groupBarrier, rt_gcMaskSize, rt_gcMask
  use RadTrans_data, ONLY: rt_boltz, rt_speedlt, rt_speedlt3, rt_radconst, rt_meshMe, &
       rt_meshCopyCount, rt_acrossMe, rt_globalComm , &
       rt_dbgContext
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_fillGuardCells, &
       Grid_getDeltas, GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       GRID_PDE_BND_DIRICHLET, &
       Grid_dfsvCreateSystem, Grid_dfsvBeginSystem, Grid_dfsvAddToSystem, &
       Grid_dfsvCompleteSystem, &
       Grid_ascAllocMem,     Grid_ascDeallocMem,   Grid_ascStart, &
       Grid_ascGetBlkPtr,    Grid_ascReleaseBlkPtr
  use Diffuse_interface, ONLY: Diffuse_solveScalar, &
       Diffuse_solveCoupledScalar, &
       Diffuse_fluxLimiter, Diffuse_setContextInfo
  use RadTrans_interface, ONLY: RadTrans_planckInt, RadTrans_sumEnergy
  use Eos_interface, ONLY: Eos_wrapped, Eos_guardCells, Eos, Eos_getTempDataF
  use Opacity_interface, ONLY: Opacity
  use Timers_interface, ONLY : Timers_start, Timers_stop 
#ifdef DEBUG_GRID_GCMASK
  use Logfile_interface, ONLY: Logfile_stampVarMask
#endif

  implicit none

#include "Flash_mpi.h"
#include "Flash.h"
#include "Eos.h"
#include "constants.h"

  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)
  real,    intent(in) :: dt
  integer, intent(in), optional :: pass


  ! Local variables:
  integer :: lb, g, gloc, gvar, iv, i, j, k
  integer :: avar, cvar, dvar, d0var
  integer :: blkLimitsGC(LOW:HIGH,MDIM), blkLimits(LOW:HIGH,MDIM)
  integer :: bcTypes(6)
  integer :: vecLen, mode
  integer :: retryCount
  integer,dimension(VARDESC_SIZE) :: ueleVarDesc,groupVarsDesc,diffCoeffDesc, &
                                     absorpCoeffDesc
  integer,dimension(VARDESC_SIZE) :: emissCoeffDesc ,emissTermDesc
  integer :: groupChunkSize, numGroupChunks, numVarsForHypre
  integer :: numGroupsLoc, numGroups, ch
  real    :: absorb_opac, emit_opac, trans_opac
  real    :: mgdTheta
  real,allocatable :: bcValues(:,:,:)
  real,allocatable :: eosData(:,:)
  real    :: cveleV
  logical :: eosMask(EOS_VARS+1:EOS_NUM)
  data eosMask /EOS_DERIVS * .FALSE./ 
  real    :: erad, eradg, urad, rho, emiss, tele, change
  real    :: xg, xgp1, pxg, pxgp1
  real, pointer   :: blkPtr(:,:,:,:)
  real    :: f
#if defined(USE_FLUXLIMITER_FCB)||defined(USE_SOLVESCALAR_FCB)
  real, POINTER, DIMENSION(:,:,:,:) :: facBptrX, facBptrY, facBptrZ
#endif

  real, parameter        :: thetaDstepPar = 0.5 ! Not yet used
  real                   :: thetaC
  real                   :: thetaD
  real                   :: thetaDstep

  ! Physical constants:
  real :: C  ! Speed of light [cm/s]
  real :: A  ! Radiation constant [ergs/cm^3/K^4]
  real :: KB ! Boltzmann constant [ergs/K]

  character(len=MAX_STRING_LENGTH) :: grp_timer
  integer :: ierr

#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,parameter :: gcMaskLogged =.TRUE.
#endif
#if defined(USE_FLUXLIMITER_FCB) && !defined(USE_SOLVESCALAR_FCB)
  logical, parameter :: returnCcB = .TRUE.
#else
  logical, parameter :: returnCcB = .FALSE.
#endif
  integer :: kk1d, kk2d, kk3d
  !=========================================================================
  if(.not. rt_useMGD) return

  call Timers_start("RadTrans") 

  kk1d = 0; kk2d = 0; kk3d = 0
  if (returnCcB) then
     kk1d = 1; kk2d = K2D; kk3d = K3D
  end if

  ! Load physical constants:
  KB = rt_boltz
  C  = rt_speedlt
  A  = rt_radconst

  ! The specific radiation energy (ERAD_VAR) may have changed due to
  ! hydrodynamic work. The group specific internal energies must be
  ! updated accordingly. This involves adding up the energy in each
  ! group and comparing it to ERAD_VAR. I will use MGDC_VAR as a
  ! temporary storage location for the sum of the group energies.

  ! Zero out MGDC_VAR:
  do lb = 1, nblk
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              blkPtr(MGDC_VAR, i, j, k) = 0.0
           enddo
        end do
     end do
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
  end do

  ! Compute the total radiation energy from the group energies:
  do gloc = 1, NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)
     gvar = MGDR_NONREP_LOC2UNK(gloc)

     do lb = 1, nblk
        call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blklst(lb), blkPtr)

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 blkPtr(MGDC_VAR,i,j,k) = &
                      blkPtr(MGDC_VAR,i,j,k) + blkPtr(gvar,i,j,k)
              enddo
           end do
        end do
        call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     end do
  end do

  call RadTrans_sumEnergy(MGDC_VAR, nblk, blklst)

  ! Scale the group energies according to ERAD_VAR:
  do gloc = 1, NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)
     gvar = MGDR_NONREP_LOC2UNK(gloc)

     do lb = 1, nblk
        call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blklst(lb), blkPtr)

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 if (blkPtr(MGDC_VAR,i,j,k) == 0.0) then
                    f = 0.0
                 else
                    f = blkPtr(ERAD_VAR,i,j,k) / blkPtr(MGDC_VAR,i,j,k)
                 end if
                 blkPtr(gvar,i,j,k) = blkPtr(gvar,i,j,k) * f

!!$                 ! While we are here, convert from specific energies
!!$                 ! to energy densities. This is reversed after the
!!$                 ! diffusion call.
!!$                 ! Nope, moved further down now. - KW
!!$                 blkPtr(gvar,i,j,k) = blkPtr(gvar,i,j,k) * blkPtr(DENS_VAR,i,j,k)
              enddo
           end do
        end do
        call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     end do
  end do

  eosMask(EOS_CVELE) = .NOT. rt_tightIonCoupling
  eosMask(EOS_CV )   = .TRUE.
  eosMask(EOS_DET)   = .TRUE.
  ! Zero out the radiation energy change and radiation temperature:
  do lb = 1, nblk
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           allocate(eosData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),EOS_NUM))
           call Eos_getTempDataF(IAXIS,(/blkLimits(LOW,IAXIS),j,k/),&
                (1-blkLimits(LOW,IAXIS)+blkLimits(HIGH,IAXIS)), &
                blkPtr, CENTER, eosData, MODE_DENS_TEMP_GATHER)
           eosData(:,EOS_DENS) = blkPtr(DENS_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), j, k)
           call Eos(MODE_DENS_TEMP_GATHER,size(eosData,1),&
                eosData, &
                massfrac=transpose(blkPtr(SPECIES_BEGIN:SPECIES_END,&
                         blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),j,k)),&
                mask=eosMask)
           if (rt_tightIonCoupling) then ! always use CV, presumably this includes contributions
                                         ! from physical ions and electrons
              blkPtr(EMIS_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), j, k) = eosData(:,EOS_CV)
           else
              blkPtr(EMIS_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), j, k) = eosData(:,EOS_CVELE)
              where (eosData(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),EOS_CVELE) == 0.0) ! fallback
                 blkPtr(EMIS_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), j, k) = eosData(:,EOS_CV)
              end where
           end if
           deallocate(eosData)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              blkPtr(EMIS_VAR, i, j, k) = blkPtr(EMIS_VAR, i, j, k) * blkPtr(DENS_VAR, i, j, k)
              blkPtr(MGDC_VAR, i, j, k) = 0.0
              blkPtr(ERAD_VAR, i, j, k) = 0.0
           enddo
        end do
     end do
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
  end do

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(rt_gcMask, .FALSE., '[RadTrans]', 'gcNeed')
  end if
#endif
  call Grid_fillGuardCells(CENTER,ALLDIR,minLayers=2,selectBlockType=LEAF, &
                           maskSize=rt_gcMaskSize, mask=rt_gcMask, &
                           doLogMask=.NOT.gcMaskLogged)
#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
#endif

  do lb = 1, nblk
     call Eos_guardCells(MODE_DENS_EI_RECAL_GATHER,blklst(lb),corners=returnCcB,&
       layers=(/1+kk1d,1+kk2d,1+kk3d/), skipSrl=.TRUE.)
  end do

  thetaDstep = thetaDstepPar
  thetaC = rt_mgdthetaC
  thetaD = rt_mgdthetaD
  mgdTheta = rt_mgdthetaImplct
  if (mgdtheta==0.0) thetaC = 0.0 !for now...
  if (mgdtheta==0.0) thetaD = 0.0 !for now...

  numVarsForHypre = 1 + rt_mgdnumgroups
  numGroupsLoc = NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)

  ! store eele*rho in MGDC_VAR
  do lb = 1, nblk
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              rho   = blkPtr(DENS_VAR, i, j, k)
              blkPtr(MGDC_VAR, i, j, k) = blkPtr(EELE_VAR, i, j, k) * rho
           enddo
        end do
     end do
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
  end do

  ueleVarDesc = (/MGDC_VAR,1,CENTER,VD_DUR_PERM/)
  call Grid_dfsvCreateSystem(baseVarDesc=ueleVarDesc,ntotVars=numVarsForHypre, &
       dt=dt, &
       theta=mgdTheta,thetaC=thetaC,thetaD=thetaD, maxChunkSize=groupChunkSize)
  if (groupChunkSize < 0) then
     groupChunkSize = numGroupsLoc
  else
     groupChunkSize = min(groupChunkSize,numGroupsLoc)
  end if
  numGroupChunks = (numGroupsLoc-1) / groupChunkSize + 1 !FOR NOW

  call Grid_dfsvBeginSystem(baseVarDesc=ueleVarDesc)

  if(rt_groupBarrier) call MPI_Barrier(rt_globalComm, ierr)

  call Grid_ascStart
#if defined(USE_FLUXLIMITER_FCB)||defined(USE_SOLVESCALAR_FCB)
  call Grid_ascAllocMem(FACEX,&
       1, &
       groupChunkSize+1, &
       0,kk1d,kk1d, &
       blklst(1:nblk) )
#if NDIM > 1
  call Grid_ascAllocMem(FACEY,&
       1, &
       groupChunkSize+1, &
       0,kk2d,kk2d, &
       blklst(1:nblk) )
#endif
#if NDIM > 2
  call Grid_ascAllocMem(FACEZ,&
       1, &
       groupChunkSize+1, &
       0,kk3d,kk3d, &
       blklst(1:nblk) )
#endif
#endif

  call Grid_ascAllocMem(SCRATCH,&
       1, &
       2*groupChunkSize, &
       0,0,0, &
       blklst(1:nblk) )
  ! Loop over energy group chunks:
  do ch = 1, numGroupsLoc, groupChunkSize
     numGroups = min(numGroupsLoc-ch+1,groupChunkSize) !number of groups in this chunk
     allocate(bcValues(2,numGroups,6))

     gvar  = MGDR_NONREP_LOC2UNK(ch)
     cvar  = CMGD_NONREP_LOC2UNK(ch)
     dvar  = DMGD_NONREP_LOC2UNK(ch)
     d0var = BMGD_NONREP_LOC2UNK(ch)
     groupVarsDesc =   (/gvar,    numGroups,CENTER,VD_DUR_PERM/)
     absorpCoeffDesc = (/cvar,    numGroups,CENTER,VD_DUR_PERM/)
     emissCoeffDesc  = (/dvar,    numGroups,CENTER,VD_DUR_PERM/)
     emissTermDesc   = (/d0var,   numGroups,CENTER,VD_DUR_PERM/)

     ! Loop over energy groups:
     do gloc = ch, ch+numGroups-1
        iv = gloc - ch

        gvar = MGDR_NONREP_LOC2UNK(gloc)
        g = NONREP_LOC2GLOB(gloc, rt_acrossMe, rt_meshCopyCount)

        if(rt_timeGroups) then 
           write(grp_timer, '(a,i4)') "Group ",g
           call Timers_start(trim(grp_timer)) 
        end if

        cvar  = CMGD_NONREP_LOC2UNK(gloc)
        dvar  = DMGD_NONREP_LOC2UNK(gloc)
        d0var = BMGD_NONREP_LOC2UNK(gloc)

        ! Setup boundary conditions:
        bcTypes(:) = rt_mgdDomainBC(g,:)
        bcValues(1,1+iv,:) = rt_mgdBcVals(g,:)
        bcValues(2,1+iv,:) = -1.0

        where (bcTypes == PERIODIC)
           bcTypes = GRID_PDE_BND_PERIODIC
        elsewhere (bcTypes == DIRICHLET)
           bcTypes = GRID_PDE_BND_DIRICHLET
        elsewhere (bcTypes == OUTFLOW .or. bcTypes == REFLECTING)
           bcTypes = GRID_PDE_BND_NEUMANN
        elsewhere (bcTypes == VACUUM)
           bcTypes = VACUUM
        end where

        ! Set the various terms that go into the diffusion calculation
        ! and compute the emission term:
        do lb = 1, nblk
           call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
           call Grid_getBlkPtr(blklst(lb), blkPtr)

           do k = blkLimits(LOW,KAXIS)-2*K3D, blkLimits(HIGH,KAXIS)+2*K3D
              do j = blkLimits(LOW,JAXIS)-2*K2D, blkLimits(HIGH,JAXIS)+2*K2D
                 do i = blkLimits(LOW,IAXIS)-2, blkLimits(HIGH,IAXIS)+2
                    ! First convert from specific energies to energy
                    ! densities. Note that this should happen after the
                    ! Grid_fillGuardCells call to ensure that the
                    ! radiation energies are handled consistently with
                    ! other PER_MASS variables. The conversion is
                    ! reversed after the diffusion call.
                    blkPtr(gvar,i,j,k) = blkPtr(gvar,i,j,k) * blkPtr(DENS_VAR,i,j,k)
                 enddo
              enddo
           enddo

           do k = blkLimits(LOW,KAXIS)-K3D-kk3d, blkLimits(HIGH,KAXIS)+K3D+kk3d
              do j = blkLimits(LOW,JAXIS)-K2D-kk2d, blkLimits(HIGH,JAXIS)+K2D+kk2d
                 do i = blkLimits(LOW,IAXIS)-1-kk1d, blkLimits(HIGH,IAXIS)+1+kk1d

                    ! Get the opacities:
                    call Opacity(blkPtr(:,i,j,k), g, absorb_opac, emit_opac, trans_opac)

                    ! Set the diffusion coefficient before application of flux limiter:
                    blkPtr(COND_VAR, i, j, k) = rt_speedlt3/max(TINY(trans_opac)*rt_speedlt3,trans_opac)

                    if (k.GE.blkLimits(LOW,KAXIS) .AND. k.LE.blkLimits(HIGH,KAXIS) .AND. &
                        j.GE.blkLimits(LOW,JAXIS) .AND. j.LE.blkLimits(HIGH,JAXIS) .AND. &
                        i.GE.blkLimits(LOW,IAXIS) .AND. i.LE.blkLimits(HIGH,IAXIS)) then
                       ! Set the leading coefficient ('factor A' in Diffuse_solveScalar) to 1 for
                       ! radiation diffusion:
                       blkPtr(DFCF_VAR, i, j, k) = 1.0

                       ! Set the absorption term:
                       blkPtr(ABSR_VAR, i, j, k) = C * absorb_opac
                       blkPtr(cvar    , i, j, k) = C * absorb_opac

                       ! Compute the emission term:

                       tele  = blkPtr(TELE_VAR, i, j, k)
                       if (KB*tele > 0.0) then
                          xg    = rt_mgdBounds(g)   / (KB*tele)
                          xgp1  = rt_mgdBounds(g+1) / (KB*tele)
                          call RadTrans_planckInt(xg, pxg)
                          call RadTrans_planckInt(xgp1, pxgp1)
                          emiss = emit_opac * C * A*(tele/PI)**4   &
                               * 15.0 * (pxgp1 - pxg)
                       else
                          emiss = 0.0
                       end if
                       cveleV = blkPtr(EMIS_VAR, i, j, k)
!!$                    blkPtr(EMIS_VAR, i, j, k) = emiss
                       blkPtr(d0var,    i, j, k) = emiss
                       blkPtr(dvar,     i, j, k) = 0.0 !!!DEV: need T-deriv of B_g here!! or use: C * absorb_opac
                       if (emiss .NE. 0.0) then
                          blkPtr(dvar,     i, j, k) = 4.0 * emiss/max(tele,1.e-20)/cveleV  !!!DEV: want better approach...
                       else
                          blkPtr(dvar,     i, j, k) = 0.0
                       end if

#if (0)
                       change = blkPtr(MGDC_VAR, i, j, k)
                       change = change - emiss
                       blkPtr(MGDC_VAR, i, j, k) = change
#endif
                    end if

                    ! Set the value of the flux limit:
                    blkPtr(FLLM_VAR,i,j,k) = rt_mgdFlCoef * C * blkPtr(gvar,i,j,k)
                    blkPtr(FLLM_VAR,i,j,k) = max(0.0, blkPtr(FLLM_VAR,i,j,k))
                 enddo
                 if (k.GE.blkLimits(LOW,KAXIS) .AND. k.LE.blkLimits(HIGH,KAXIS) .AND. &
                     j.GE.blkLimits(LOW,JAXIS) .AND. j.LE.blkLimits(HIGH,JAXIS)) then
!!$                    print*,'cvar:',blkPtr(cvar,     5:12, j, k),blklst(lb)
!!$                    print*,'dvar:',blkPtr(dvar,     5:12, j, k)
                 end if
              enddo
           enddo
           call Grid_releaseBlkPtr(blklst(lb), blkPtr)
        end do

     ! Apply the flux limiter: (also returns 3*lambda (flux limiter factor) in FLLM_VAR)
#     ifdef USE_FLUXLIMITER_FCB
        call Diffuse_fluxLimiterFcB(COND_VAR, gvar, FLLM_VAR, rt_mgdFlMode, nblk, blklst, &
             returnCcB)

#     else
        call Diffuse_fluxLimiter(COND_VAR, gvar, FLLM_VAR, rt_mgdFlMode, nblk, blklst)

#      ifdef USE_SOLVESCALAR_FCB
        do lb = 1, nblk
           call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
           call Grid_getBlkPtr(blklst(lb), blkPtr)
           call Grid_ascGetBlkPtr(blklst(lb),facBptrX,FACEX)
           do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                 do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)+1
                    facBptrX(1+iv,i,j,k) = 0.5* &
                         ( blkPtr(COND_VAR, i-1, j, k) + blkPtr(COND_VAR, i, j, k) )
                 end do
              end do
           end do
           call Grid_ascReleaseBlkPtr(blklst(lb),facBptrX,FACEX)
#if NDIM >= 2
           call Grid_ascGetBlkPtr(blklst(lb),facBptrY,FACEY)
           do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)+1
                 do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                    facBptrY(1+iv,i,j,k) = 0.5* &
                         ( blkPtr(COND_VAR, i, j-1, k) + blkPtr(COND_VAR, i, j, k) )
                 end do
              end do
           end do
           call Grid_ascReleaseBlkPtr(blklst(lb),facBptrY,FACEY)
#if NDIM == 3
           call Grid_ascGetBlkPtr(blklst(lb),facBptrZ,FACEZ)
           do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)+1
              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                 do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                    facBptrZ(1+iv,i,j,k) = 0.5* &
                         ( blkPtr(COND_VAR, i, j, k-1) + blkPtr(COND_VAR, i, j, k) )
                 end do
              end do
           end do
           call Grid_ascReleaseBlkPtr(blklst(lb),facBptrZ,FACEZ)
#endif
#endif     
           call Grid_releaseBlkPtr(blklst(lb), blkPtr)
        end do
#      endif
#     endif
     end do

     ! Add matrix elements for this chunk of groups:
#ifdef USE_SOLVESCALAR_FCB
     diffCoeffDesc = (/1,numGroups,FACES,VD_DUR_GASC/)
     call Grid_dfsvAddToSystem(unkVarsDesc=groupVarsDesc,baseVarDesc=ueleVarDesc, &
                               firstHypreVar=ch, &
                               diffCoeffDesc=diffCoeffDesc, &
                               absorpCoeffDesc=absorpCoeffDesc, &
                               emissCoeffDesc=emissCoeffDesc, &
                               emissTermDesc=emissTermDesc, &
                               bcTypes=bcTypes, bcValues=bcValues)
#if (0)
  call Diffuse_solveScalar(   &
          gvar, -1, DFCF_VAR, &
          bcTypes, bcValues, dt, 1.0, 1.0, mgdtheta, &
          pass, nblk,blklst, ABSR_VAR, EMIS_VAR)
#endif

#else
     diffCoeffDesc = (/COND_VAR,numGroups,CENTER,VD_DUR_PERM/)
     call Grid_dfsvAddToSystem(unkVarsDesc=groupVarsDesc,baseVarDesc=ueleVarDesc, &
                               firstHypreVar=ch, &
                               diffCoeffDesc=diffCoeffDesc, &
                               absorpCoeffDesc=absorpCoeffDesc, &
                               emissCoeffDesc=emissCoeffDesc, &
                               emissTermDesc=emissTermDesc, &
                               bcTypes=bcTypes, bcValues=bcValues)
#endif



     deallocate(bcValues)

     if(rt_timeGroups) call Timers_stop(trim(grp_timer))
  end do
  
  gvar  = MGDR_NONREP_LOC2UNK(1)
  cvar =  CMGD_NONREP_LOC2UNK(1)
  dvar  = DMGD_NONREP_LOC2UNK(1)
  d0var = BMGD_NONREP_LOC2UNK(1)
  groupVarsDesc =   (/gvar,    numGroupsLoc,CENTER,VD_DUR_PERM/)
  absorpCoeffDesc = (/cvar,    numGroupsLoc,CENTER,VD_DUR_PERM/)
  emissCoeffDesc  = (/dvar,    numGroupsLoc,CENTER,VD_DUR_PERM/)
  emissTermDesc   = (/d0var,   numGroupsLoc,CENTER,VD_DUR_PERM/)

  call Grid_dfsvCompleteSystem(baseVarDesc=ueleVarDesc)

  if (rt_groupBarrier) then
     call Timers_start("RadTrans load imbalance")
     call MPI_Barrier(rt_globalComm, ierr)
     call Timers_stop("RadTrans load imbalance")     
  end if

  call Diffuse_setContextInfo(group=0, component=3)

  mgdTheta = rt_mgdthetaImplct
  retryCount = 0
  ! Change the following definition to 0 to disable all the retrystuff - KW
#define MAXRETRY 0
111 continue

  rt_dbgContext%willingToRetry = (retryCount < MAXRETRY)

  call Diffuse_solveCoupledScalar(unkVarsDesc=groupVarsDesc,xtraVarDesc=ueleVarDesc,diffCoeffDesc=diffCoeffDesc, &
                               absorpCoeffDesc=absorpCoeffDesc, &
                               bcTypes=bcTypes, dt=dt, &
                               scaleFact=1.0, theta=mgdtheta, &
                               pass=pass, blockCount=nblk,blockList=blklst &
!!$                               ,iFactorD=ABSR_VAR, iFactorA=DFCF_VAR &
                               )

     
900 format(1x,'[RadTrans] ***** ',A,' IN Diffuse_solveScalar !!! ****',1x,6(I14:1x),L1,(I6))
902 format(1x,'[RadTrans] ** ',A,' IN Diffuse_solveScalar **',1x,6(I14:1x),L1,(I6))
  if (rt_dbgContext%flashErrCode.NE.0) then
     if (rt_meshMe==MASTER_PE) then
        select case (rt_dbgContext%flashErrCode)
        case(2)
           if (retryCount>0) &
                print 900,'WARNING', rt_dbgContext, rt_meshMe
        case(3)
           if (retryCount>0) &
                print 902,'INFO', rt_dbgContext, rt_meshMe
        case default
           print 900,'ERROR DETECTED', rt_dbgContext, rt_meshMe
        end select
     end if
     if (rt_dbgContext%retriable .NE. 0 .AND. rt_dbgContext%willingToRetry) then
        if (retryCount < MAXRETRY) then
           retryCount = retryCount + 1
           if (retryCount==MAXRETRY) then
              f = 0.0
           else
              f = 0.1
           end if
           do lb = 1, nblk
              call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
              call Grid_getBlkPtr(blklst(lb), blkPtr)
              do k = blkLimits(LOW,KAXIS)-K3D, blkLimits(HIGH,KAXIS)+K3D
                 do j = blkLimits(LOW,JAXIS)-K2D, blkLimits(HIGH,JAXIS)+K2D
                    do i = blkLimits(LOW,IAXIS)-1, blkLimits(HIGH,IAXIS)+1
                       ! Scale down (or zero out) the diffusion coefficient:
                       blkPtr(COND_VAR, i, j, k) = f * blkPtr(COND_VAR, i, j, k)
                    enddo
                 enddo
              enddo
              call Grid_releaseBlkPtr(blklst(lb), blkPtr)
           end do
           if (rt_meshMe==MASTER_PE) then
              print*,'[RadTrans] ***** RETRY #',retryCount,' with COND_VAR scaled by',f
           end if
!!$              if (retryCount==MAXRETRY) then
!!$                 mgdTheta = 0.0
!!$              else
!!$                 mgdTheta = 0.5 * mgdTheta
!!$              end if
!!$              if (rt_meshMe==MASTER_PE) then
!!$                 print*,'[RadTrans] ***** RETRY #',retryCount,' with theta=',mgdTheta
!!$              end if
           goto 111
        end if
     end if
  end if

  call Grid_ascDeallocMem(SCRATCH)
#if defined(USE_FLUXLIMITER_FCB)||defined(USE_SOLVESCALAR_FCB)
  call Grid_ascDeallocMem(FACEX)
  if (NDIM > 1) call Grid_ascDeallocMem(FACEY)
  if (NDIM > 2) call Grid_ascDeallocMem(FACEZ)
#endif

     ! Compute the energy absorption:
  do gloc = 1, NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)
     gvar = MGDR_NONREP_LOC2UNK(gloc)
     do lb = 1, nblk
        call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blklst(lb), blkPtr)

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 
                 ! Get the group radiation energy (erg/cm^3)
                 urad = blkPtr(gvar, i, j, k)

                 ! Convert radiation energy density back into specific energy:
                 rho   = blkPtr(DENS_VAR, i, j, k)
                 eradg = urad/rho
                 blkPtr(gvar, i, j, k) = eradg

                 ! Compute total radiation specific energy:
                 erad = blkPtr(ERAD_VAR, i, j, k)
                 erad = erad + eradg
                 blkPtr(ERAD_VAR, i, j, k) = erad
              enddo
           enddo
        enddo
        call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     end do
  end do
  ! Sum radiation energy change and specific energy over all meshes:
  call RadTrans_sumEnergy(ERAD_VAR, nblk, blklst)
!!$  call RadTrans_sumEnergy(MGDC_VAR, nblk, blklst)

  ! Next, update the electron temperature to account for
  ! emission/absorption: !DEV: ???? - KW
  do lb = 1, nblk

     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              rho = blkPtr(DENS_VAR, i, j, k)
#if (0)
              change = blkPtr(MGDC_VAR, i, j, k)
              blkPtr(EELE_VAR,i,j,k) = blkPtr(EELE_VAR,i,j,k) + change*dt/rho
#else
              blkPtr(EELE_VAR,i,j,k) = blkPtr(MGDC_VAR,i,j,k) /rho
#endif
           end do
        end do
     end do

     call Grid_releaseBlkPtr(blklst(lb), blkPtr)

     call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,blklst(lb))
  end do


  ! Not needed here - KW 2012-12-10
!!$  call Grid_fillGuardCells(CENTER,ALLDIR)

  call Timers_stop("RadTrans") 

  return

end subroutine RadTrans
