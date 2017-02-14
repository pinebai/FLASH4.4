!!****if* source/physics/Diffuse/DiffuseMain/Unsplit/diff_advanceTherm
!!
!!  NAME 
!!
!!  diff_advanceTherm
!!
!!  SYNOPSIS
!!
!!  call diff_advanceTherm(integer(IN)                  :: blockCount,
!!                         integer(IN)                  :: blockList(blockCount),
!!                         real(IN)                     :: dt,
!!                         integer, OPTIONAL, intent(IN):: pass)
!!
!!  DESCRIPTION 
!!      This routine advances the heat diffusion equation (i.e., heat conduction).
!!      An implicit scheme is used. 
!!
!!      Supported boundary conditions are: 
!!                PERIODIC, OUTFLOW (tested).
!!                DIRICHLET (untested).
!!
!!
!! ARGUMENTS
!!
!!  blockCount   - The number of blocks in the list
!!  blockList(:) - The list of blocks on which the solution must be updated
!!   dt           : The time step
!!   pass         : pass=1 directional order of solution sweep X-Y-Z, 
!!                  pass=2 directional order of solution sweep Z-Y-X.
!!  dt           - The time step
!!
!! SIDE EFFECTS
!!
!!  Updates certain variables in permanent UNK storage to contain the
!!  updated temperature and some auxiliaries.  Invokes a solver (of the Heat
!!  diffusion equation). On return,
!!     TEMP_VAR:  contains updated temperature for the current simulation time.
!!     EINT_VAR, ENER_VAR, PRES_VAR, etc.: updated accordingly by EOS call
!!
!!  May modify certain variables used for intermediate results by the solvers
!!  invoked. The list of variables depends on the Diffuse implementation.
!!  The following information is subject to change.
!!     COND_VAR:  contains conductivity that was passed to Grid_advanceDiffusion
!!  
!!  NOTES:
!!  
!!  The current implementation of Radiation diffusion is a Gray Approximation.
!!  Types: P1, Simple Diffusion, Flux limited Diffusion, P1/3 and Variable Eddington Factor (VEF)
!!  We have implemented Simple Diffusion.
!!
!!
!! NOTES
!!
!!  The interface of this subroutine must be explicitly known to code that
!!  calls it.  The simplest way to make it so is to have something like
!!    use diff_interface,ONLY: diff_advanceTherm
!!  in the calling routine.
!!***

!!REORDER(4): solnVec

#define DEBUG_GRID_GCMASK

subroutine diff_advanceTherm(blockCount,blockList,dt,pass)


  use Diffuse_data, ONLY : diffusion_cutoff_density
  use Diffuse_data, ONLY : useDiffuse, diff_meshMe, diff_meshcomm,&
       diff_useEleCond, &
       diff_mele, diff_boltz, &
       diff_singleSpeciesA, diff_singleSpeciesZ, diff_avo, diff_mele, &
       diff_eleFlMode, diff_eleFlCoef
  use Diffuse_data, ONLY: diff_useIonCond
  use Diffuse_data, ONLY: diff_ionFlMode
  use Diffuse_data, ONLY: diff_ionFlCoef
  
  use diff_saData, ONLY : diff_boundary, &
       updateDiffuse, &
       diff_scaleFactThermSaTempDiff, diff_scaleFactThermSaTime, &
       diff_eleDomainBC, diff_thetaImplct, diff_updEint
  use diff_saData, ONLY: diff_ionThetaImplct
  use diff_saData, ONLY: diff_ionDomainBC
  use Eos_interface, ONLY : Eos_wrapped, Eos, Eos_getAbarZbar, Eos_getTempData
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Conductivity_interface, ONLY : Conductivity
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_advanceDiffusion, Grid_getBlkIndexLimits, Grid_fillGuardCells, &
      Grid_getDeltas, GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
      GRID_PDE_BND_DIRICHLET

#include "Flash.h"

#ifdef DEBUG_GRID_GCMASK
  use Logfile_interface, ONLY : Logfile_stampVarMask
#endif
  use Diffuse_interface, ONLY: Diffuse_solveScalar, &
       Diffuse_fluxLimiter, &
       Diffuse_setContextInfo

  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Eos.h"
#include "Eos_components.h"


  integer,intent(IN)                       :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  real,intent(in)                          :: dt
  integer, OPTIONAL, intent(IN)            :: pass
  
  integer :: Temptodiffuse, TempForEos, cvToUse
  integer :: firstComponent
  integer :: i,j,k,n
  integer :: lb, blockID
  integer :: bcTypes(6)
  integer :: mode, vecLen
  
  real    :: bcValues(2,6) = 0.
  real    :: cond_zone, diff_coeff, xdens, xtemp
  real    :: massfrac(NSPECIES), Ye 
  real    :: eos_arr(EOS_NUM)
  
  logical :: gcmask(NUNK_VARS)  
  logical :: mask(EOS_VARS+1:EOS_NUM)
  
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec  

  real :: mion
  real :: abar

  integer ::EosMode
  logical,save :: gcMaskLogged =.FALSE.
  
  !=========================================================================  

  if (.not. useDiffuse) return  
  
#ifdef TELE_VAR
  
  if (.not. (diff_useEleCond .or. diff_useIonCond) ) return
  firstComponent= EOSCOMP_ELE
  Temptodiffuse = TELE_VAR 
  cvToUse       = EOS_CVELE
  TempForEos    = EOS_TEMPELE 
  mode          = MODE_DENS_TEMP_GATHER
  
  EosMode = MODE_DENS_TEMP_GATHER
  if (diff_updEint) EosMode = MODE_DENS_EI_GATHER
  
#else  
  firstComponent= EOSCOMP_MATTER
  Temptodiffuse = TEMP_VAR
  cvToUse       = EOS_CV
  TempForEos    = EOS_TEMP  
  mode          = MODE_DENS_TEMP

  EosMode = MODE_DENS_TEMP
  if (diff_updEint) EosMode = MODE_DENS_EI
#endif  
  
  call Timers_start("diff_advanceTherm")

  
  if(diff_useEleCond .or. Temptodiffuse == TEMP_VAR) then
     ! This code handles:
     ! 1. Conduction in the 1T case
     ! 2. Electron conduction in the 3T case

     bcTypes(:) = diff_eleDomainBC(:)

     bcValues = 0.

     where (bcTypes == PERIODIC)
        bcTypes = GRID_PDE_BND_PERIODIC
     elsewhere (bcTypes == DIRICHLET)
        bcTypes = GRID_PDE_BND_DIRICHLET
     elsewhere (bcTypes == OUTFLOW)
        bcTypes = GRID_PDE_BND_NEUMANN
     end where

     gcmask(:) = .FALSE.
     gcmask(Temptodiffuse) = .TRUE.


#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(gcmask, .FALSE., '[diff_advanceTherm]', 'gcNeed')
     end if
#endif
     call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS, &
          mask=gcmask,selectBlockType=LEAF,doLogMask=.NOT.gcMaskLogged)   

     do lb = 1, blockCount
        blockID = blockList(lb)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blockID, solnVec)

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 xtemp = solnVec(Temptodiffuse,i,j,k)
                 xdens = solnVec(DENS_VAR,i,j,k)

                 ! if the density is too low, turn off thermal diffusion
                 if (xdens < diffusion_cutoff_density) then
                    cond_zone = 0.0
                 else
                    call Conductivity(solnVec(:,i,j,k), cond_zone, component=EOSCOMP_ELE)
                 end if


                 ! Compute the specific heat:
                 do n = 1, NSPECIES
                    massfrac(n) = solnVec(SPECIES_BEGIN-1+n,i,j,k)
                 enddo

                 vecLen = 1
                 call Eos_getTempData(IAXIS,(/i,j,k/),vecLen, &
                      solnVec, CENTER, eos_arr, mode)
                 eos_arr(EOS_DENS) = xdens
                 mask = .false.
                 mask(cvToUse) = .true.
                 mask(EOS_DET) = .true.
                 call Eos(mode,vecLen,eos_arr,massfrac,mask)

                 solnVec(COND_VAR,i,j,k) = cond_zone
                 solnVec(DFCF_VAR,i,j,k) = xdens*eos_arr(cvToUse)

                 ! Set the flux limiter:
                 call Eos_getAbarZbar(solnVec(:,i,j,k),Ye=Ye)
                 solnVec(FLLM_VAR,i,j,k) = diff_eleFlCoef * &
                      sqrt(diff_boltz*xtemp/diff_mele) * &
                      diff_boltz*xtemp * &
                      (Ye * diff_avo * xdens)
              enddo
           enddo
        enddo

        call Grid_releaseBlkPtr(blockID, solnVec)
     end do

     call Grid_fillGuardCells(CENTER,ALLDIR)

     call Diffuse_setContextInfo(group=0, component=firstComponent)

     call Diffuse_fluxLimiter(COND_VAR, Temptodiffuse, FLLM_VAR, &
          diff_eleFlMode, blockCount, blockList)

     if (diff_updEint) then
        call diff_getFaceFluxes (blockCount, blockList, Temptodiffuse, & 
             COND_VAR, 0, 1.0-diff_thetaImplct) 
     end if

     call Diffuse_solveScalar(Temptodiffuse,COND_VAR, DFCF_VAR, bcTypes,    &
          bcValues, dt, 1.0, 1.0,           &   
          diff_thetaImplct,pass, blockCount,blockList)

     if (diff_updEint) then

        call Timers_start("diff_updEint")

#ifdef DEBUG_GRID_GCMASK
        if (.NOT.gcMaskLogged) then
           call Logfile_stampVarMask(gcmask, .FALSE., '[diff_advanceTherm]', 'gcNeed')
        end if
#endif
        call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS, &
             mask=gcmask,selectBlockType=LEAF,doLogMask=.NOT.gcMaskLogged)

        call diff_getFaceFluxes (blockCount, blockList, Temptodiffuse, & 
             COND_VAR, 1, diff_thetaImplct) 

        call diff_updateEnergy (blockCount, blockList, dt) 

        call Timers_stop("diff_updEint")

     end if
  end if

  ! **************************
  ! *                        *
  ! *     ION CONDUCTION     *
  ! *                        *
  ! **************************
#ifdef FLASH_3T
  if(diff_useIonCond) then
  
     bcTypes(:) = diff_ionDomainBC(:)

     bcValues = 0.

     where (bcTypes == PERIODIC)
        bcTypes = GRID_PDE_BND_PERIODIC
     elsewhere (bcTypes == DIRICHLET)
        bcTypes = GRID_PDE_BND_DIRICHLET
     elsewhere (bcTypes == OUTFLOW)
        bcTypes = GRID_PDE_BND_NEUMANN
     end where

     gcmask(:) = .FALSE.
     gcmask(TION_VAR) = .TRUE.
     call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS, &
          mask=gcmask,selectBlockType=LEAF,doLogMask=.NOT.gcMaskLogged)   

     do lb = 1, blockCount
        blockID = blockList(lb)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blockID, solnVec)

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 ! if the density is too low, turn off thermal diffusion
                 if (solnVec(DENS_VAR,i,j,k) < diffusion_cutoff_density) then
                    cond_zone = 0.0
                 else
                    ! Get the ion conductivity:
                    call Conductivity(solnVec(:,i,j,k), cond_zone, component=EOSCOMP_ION)
                 end if

                 ! Compute the specific heat:
                 do n = 1, NSPECIES
                    massfrac(n) = solnVec(SPECIES_BEGIN-1+n,i,j,k)
                 enddo

                 eos_arr(EOS_DENS) = solnVec(DENS_VAR,i,j,k)
                 eos_arr(EOS_TEMPION) = solnVec(TION_VAR,i,j,k)
                 eos_arr(EOS_TEMPELE) = solnVec(TELE_VAR,i,j,k)
                 eos_arr(EOS_TEMPRAD) = solnVec(TRAD_VAR,i,j,k)
                 mask = .false.
                 mask(EOS_CVION) = .true.
                 mask(EOS_DET) = .true.
                 call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,massfrac,mask)

                 solnVec(COND_VAR,i,j,k) = cond_zone
                 solnVec(DFCF_VAR,i,j,k) = solnVec(DENS_VAR,i,j,k)*eos_arr(EOS_CVION)

                 ! Set the flux limiter:
                 call Eos_getAbarZbar(solnVec(:,i,j,k), abar=abar)
                 mion = abar / diff_avo
                 solnVec(FLLM_VAR,i,j,k) = diff_ionFlCoef * &
                      sqrt(diff_boltz*solnVec(TION_VAR,i,j,k)/mion) * &
                      diff_boltz*solnVec(TION_VAR,i,j,k) * &
                      (solnVec(DENS_VAR,i,j,k)/mion)
              enddo
           enddo
        enddo

        call Grid_releaseBlkPtr(blockID, solnVec)
     end do

     call Grid_fillGuardCells(CENTER,ALLDIR)

     call Diffuse_setContextInfo(group=0, component=EOSCOMP_ION)

     call Diffuse_fluxLimiter(COND_VAR, TION_VAR, FLLM_VAR, &
          diff_ionFlMode, blockCount, blockList)

     call Diffuse_solveScalar(TION_VAR,COND_VAR, DFCF_VAR, bcTypes, &
          bcValues, dt, 1.0, 1.0, &   
          diff_ionThetaImplct,pass, blockCount,blockList)

  end if
#endif

  ! ********************
  ! *                  *  
  ! *     CALL EOS     *
  ! *                  *
  ! ********************
  call Timers_start("Eos_wrapped")  
  do lb =1, blockCount
     blockID = blockList(lb)     
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)                   
     call Eos_wrapped(EosMode, blkLimits, blockList(lb))         
  enddo  
  call Timers_stop("Eos_wrapped")
  
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
  
  call Timers_stop ("diff_advanceTherm")
  
  return
  
end subroutine diff_advanceTherm
