!!****if* source/physics/Diffuse/DiffuseMain/Split/diff_advanceTherm
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

subroutine diff_advanceTherm(blockCount,blockList,dt,pass)


  use Diffuse_data, ONLY : useDiffuse, diff_meshMe, diff_meshcomm,&
       diff_useEleCond, diff_asol, &
       diff_mele, diff_boltz, &
       diff_singleSpeciesA, diff_singleSpeciesZ, diff_avo, diff_mele, &
       diff_eleFlMode, diff_eleFlCoef
  
  use diff_saData, ONLY : diff_boundary, &
       updateDiffuse, &
       diff_scaleFactThermSaTempDiff, diff_scaleFactThermSaTime, &
       diff_eleDomainBC, diff_thetaImplct
  use Eos_interface, ONLY : Eos_wrapped, Eos, Eos_getAbarZbar, Eos_getTempData
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Conductivity_interface, ONLY : Conductivity
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_advanceDiffusion, Grid_getBlkIndexLimits, Grid_fillGuardCells, &
      Grid_getDeltas, GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
      GRID_PDE_BND_DIRICHLET
  use Diffuse_interface, ONLY: Diffuse_solveScalar, &
       Diffuse_fluxLimiter

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"
#include "Eos.h"


  integer,intent(IN)                       :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  real,intent(in)                          :: dt
  integer, OPTIONAL, intent(IN)            :: pass
  
  integer :: Temptodiffuse, TempForEos, cvToUse
  integer :: i,j,k,n
  integer :: lb, blockID
  integer :: bcTypes(6)
  integer :: mode, vecLen
  
  real    :: outputScaleFact
  real    :: bcValues(2,6) = 0.
  real    :: chi
  real    :: cond_zone, diff_coeff, xdens, xtemp
  real    :: massfrac(NSPECIES), Ye 
  real    :: eos_arr(EOS_NUM)
  
  logical :: gcmask(NUNK_VARS)  
  logical :: mask(EOS_VARS+1:EOS_NUM)
  
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec  
  
  integer ::EosMode

  !=========================================================================  
  
#ifdef TELE_VAR
  
  if (.not. diff_useEleCond) return  
  Temptodiffuse = TELE_VAR 
  cvToUse       = EOS_CVELE
  TempForEos    = EOS_TEMPELE 
  mode          = MODE_DENS_TEMP_GATHER
  
  EosMode = MODE_DENS_TEMP_GATHER

#else  
  
  if (.not. useDiffuse) return
  Temptodiffuse = TEMP_VAR
  cvToUse       = EOS_CV
  TempForEos    = EOS_TEMP  
  mode          = MODE_DENS_TEMP

  EosMode = MODE_DENS_TEMP
#endif  
  
  call Timers_start("diff_advanceTherm")
  
  bcTypes(:) = diff_eleDomainBC(:)
  
  bcValues = 0.
  
  where (bcTypes == PERIODIC)
     bcTypes = GRID_PDE_BND_PERIODIC
  elsewhere (bcTypes == DIRICHLET)
     bcTypes = GRID_PDE_BND_DIRICHLET
  elsewhere (bcTypes == OUTFLOW)
     bcTypes = GRID_PDE_BND_NEUMANN
  end where
  
  outputScaleFact = 1.0  

  gcmask(Temptodiffuse) = .TRUE.
  
  call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS, &
       mask=gcmask,selectBlockType=LEAF)   
  
  do lb = 1, blockCount
     
     blockID = blockList(lb)
     
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     
     call Grid_getBlkPtr(blockID, solnVec)
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              xtemp = solnVec(Temptodiffuse,i,j,k)
              xdens = solnVec(DENS_VAR,i,j,k)

              ! load the mass fractions
              do n = 1, NSPECIES
                 massfrac(n) = solnVec(SPECIES_BEGIN-1+n,i,j,k)
              enddo
              
              call Conductivity(solnVec(:,i,j,k), cond_zone, component=2)
              
              ! Compute CV
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
              
              ! Set abar and zbar:
#ifdef FLASH_MULTISPECIES
              call Eos_getAbarZbar(solnVec(:,i,j,k),Ye=Ye,massFrac=massfrac)
#else
#ifdef YE_MSCALAR
              Ye = solnVec(YE_MSCALAR,i,j,k)
#else
              Ye = diff_singleSpeciesZ / diff_singleSpeciesA
#endif
#endif           
              ! Set electron flux limiter:
              solnVec(FLLM_VAR,i,j,k) = diff_eleFlCoef * &
                   sqrt(diff_boltz*xtemp/diff_mele) * &
                   diff_boltz*xtemp * &
                   (Ye * diff_avo * xdens)
           enddo
        enddo
     enddo
     
     call Grid_releaseBlkPtr(blockID, solnVec)
  end do
     
  call Diffuse_fluxLimiter(COND_VAR, Temptodiffuse, FLLM_VAR, &
       diff_eleFlMode, blockCount, blockList)
     
  call Diffuse_solveScalar(Temptodiffuse,COND_VAR, DFCF_VAR, bcTypes,    &
                           bcValues, dt, outputScaleFact, chi,           &   
                           diff_thetaImplct,pass, blockCount,blockList)
  
  
  call Timers_start("Eos_wrapped")
  
  do lb =1, blockCount
     
     blockID = blockList(lb)
     
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)              
     
     call Eos_wrapped(EosMode, blkLimits, blockList(lb))    
     
  enddo
  
  call Timers_stop("Eos_wrapped")
  
  
  call Timers_stop ("diff_advanceTherm")
  
  return
  
end subroutine diff_advanceTherm
