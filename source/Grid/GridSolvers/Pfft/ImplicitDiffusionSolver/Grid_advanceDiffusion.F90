!!****if* source/Grid/GridSolvers/Pfft/ImplicitDiffusionSolver/Grid_advanceDiffusion
!!
!! NAME
!!  Grid_advanceDiffusion
!!
!! SYNOPSIS
!!
!!  call Grid_advanceDiffusion(integer(IN)           :: iSoln,
!!                             integer(IN)           :: iSrc,
!!                             integer(IN)           :: iCond,
!!                             integer(IN)           :: iConstFact,
!!                             integer(IN)           :: bcTypes(6),
!!                             real(IN)              :: bcValues(2,6),
!!                             real(IN)              :: dt,
!!                             real(IN)              :: chi,
!!                             real(IN)              :: scaleFact,
!!                             logical(IN)           :: solnIsDelta,
!!                             integer(IN), OPTIONAL :: pass)
!!
!! DESCRIPTION
!!
!!  This routine advances a diffusive operator (think head conduction) of the form
!!
!!      A * dV/dt = d/dx(B*dV/dx) + d/dy(B*dV/dy) + d/dx(B*dV/dz)
!!      Heat conduction : rho*cv dT/dt = d/dx(kdT/dx) + ..      
!!  
!! ARGUMENTS
!!
!!  iSoln         : The index for the solution variable (V, Temperature in heat conduction)
!!  iSrc          : discretized RHS (d/dx(B*dV/dx) + d/dy(B*dV/dy) + d/dx(B*dV/dz))
!!  iCond         : The index of the diffusion coefficient (if used, otherwise -1)
!!  iConstFact    : The index of coefficient of temporal derivative (think A)
!!                  rho*cv in Heat conduction.
!!  bcTypes       : Left and Right boundary condition type.
!!                  Supports PERIODIC, OUTFLOW (tested).
!!                           DIRICHLET,PNEUMAN (untested).
!!  bcValues      : an unused argument used to keep the interface standard.
!!  dt            : the time step.
!!  chi           : a factor in the heat diffusion equation, not used anymore.
!!  scaleFact     : the scaling factor for the eventual solution.
!!  solnIsDelta   : Is the solution only a delta that the caller has to apply to the
!!                  temperature, rather than temperature itself.
!!  pass          : pass=1 order of directional sweep X-Y-Z, pass=2 order of directional sweep Z-Y-X.
!!
!!
!! NOTES
!!
!!  It is currently assumed in this implementation
!!  DEV: But not checked! - KW
!!  that boundary condition types at all sides of the domain are the same. 
!!  What this means to the user  is that only the first value in bcTypes is
!!  checked here, and may be assumed to give the type of boundary condition
!!  for all 2*NDIM directions.
!!
!! SEE ALSO
!! 
!!  Diffuse_advance1D
!!  
!!
!!***

subroutine Grid_advanceDiffusion (iVar, iSrc, iFactorB, iFactorA, bcTypes, bcValues, dt, chi, scaleFact, &
     theta, solnIsDelta, iFactorC, iFactorD, pass)       
     
  use Grid_data, ONLY        : gr_meshMe, gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax, gr_meshNumProcs
  
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY: Logfile_stamp
  use Diffuse_interface, ONLY: Diffuse_advance1D

  use gr_pfftData,      ONLY : pfft_inLen, pfft_midLen, pfft_outLen, pfft_usableProc, pfft_setupOnce, &
       pfft_work1, pfft_work2, pfft_workSize, pfft_comm, pfft_procGrid, pfft_globalLen, pfft_transformType
  
  use gr_pfftInterface, ONLY : gr_pfftTranspose, gr_pfftSpecifyTransform
  
  

  use Grid_interface, ONLY  : Grid_getGlobalIndexLimits, Grid_pfftMapToInput, Grid_pfftInit, Grid_pfftMapFromOutput, &
       Grid_pfftFinalize
  
  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "Pfft.h" 
         
  integer, intent(IN) :: iVar
  integer, intent(IN) :: iSrc
  integer, intent(IN) :: iFactorB
  integer, intent(IN) :: iFactorA
  real, intent(IN)    :: dt 
  real, intent(IN)    :: chi
  real, intent(IN)    :: scaleFact
  real, intent(IN)    :: theta
  logical, intent(IN) :: solnIsDelta
  integer, dimension(6),  intent(IN) :: bcTypes
  real   , dimension(2,6),intent(IN) :: bcValues
  integer, intent(IN), OPTIONAL :: pass
  integer, intent(IN), OPTIONAL :: iFactorC
  integer, intent(IN), OPTIONAL :: iFactorD  

  
  !=======================================================================
  integer                             :: inSize, lpass
  logical                             :: needMap, lConstFact,lmatC,lmatD
  real, dimension(:), allocatable     :: inArray, outArray
  real, dimension(:), pointer         :: ConstFactArray1, ConstFactArray2
  real, dimension(:), pointer         :: matCArray, matDArray    
  real                                :: DELX, DELY, DELZ
  integer, dimension(MDIM)            :: localSize, globalSize, transformType
  integer, dimension(0:MDIM)          :: baseDatType
  integer                             :: ML,NL,LL 
  real, pointer    , dimension(:)     :: pfft_work3, pfft_work4
  !=======================================================================

  real                                :: LbcValues(6)

  LbcValues(:) = bcValues(1,:)


  if (solnIsDelta) then
     call Logfile_stamp(&
          'This implementation only supports solnIsDelta=.FALSE.',&
          '[Grid_advanceDiffusion]')
     call Driver_abortFlash('[Grid_advanceDiffusion] solnIsDelta must be FALSE!')
  end if
  
  if (present(pass)) then
     lpass = pass
  else
     lpass = 0
  end if

  
  if (present(iFactorC)) then
     lMatC = (iFactorC > 0)
  else
     lMatC = .FALSE.
  end if
  
  if (present(iFactorD)) then
     lMatD = (iFactorD > 0)
  else
     lMatD = .FALSE.
  end if
  
  
  lConstFact = (iFactorA > 0)

  ! Initialize here
  if(.not.pfft_setupOnce) then   !THIS "IF" WILL ONLY BE INVOKED IN PARAMESH SIMULATIONS.

#ifdef FLASH_GRID_PARAMESH
     needMap=(NDIM > 1)
#else
     needMap= .TRUE.
#endif
  
     call Grid_getGlobalIndexLimits(globalSize)

     call gr_pfftSpecifyTransform(transformType, baseDatType, bcTypes)

     call Grid_pfftInit( NDIM,needMap,globalSize,&
                        localSize,transformType, baseDatType)
  else
     globalSize(:) = pfft_globalLen(:)
     transformType(:) = pfft_transformType(:)
     localSize(:) = pfft_inLen(:)
  end if

  ! Initializations performed in the stand-alone Diffuse unit implementation.
  if(.not.pfft_usableProc) return

  if ((NDIM == 1) .and. (globalSize(1) == 1)) then
     call Driver_abortFlash('[Grid_advanceDiffusion] X-GridSize must be set > 1 in 1D')      
  endif
  
  if ((NDIM == 1) .and. (globalSize(2) /= 1 .or. globalSize(3) /= 1)) then
     call Driver_abortFlash('[Grid_advanceDiffusion] Y-GridSize, Z-GridSize must be set to 1 in 1D')      
  endif
  
  if ((NDIM == 2) .and. (globalSize(3) /= 1)) then
     call Driver_abortFlash('[Grid_advanceDiffusion] Z-GridSize must be set to 1 in 1D')
  endif  
  
  inSize=pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS) 

  !! Use pfft_workSize, to allocate inArray.
  allocate( inArray(max(pfft_workSize,pfft_inLen(1))))
  allocate(outArray(max(pfft_workSize,pfft_inLen(1))))
  
  !! Factor A
  if (lConstFact) then
     allocate(ConstFactArray1(max(pfft_workSize,pfft_inLen(1))))
  else
     nullify(ConstFactArray1)
  end if

  !! Factor C
  if (lMatC) then
     allocate(matCArray(max(pfft_workSize,pfft_inLen(1))))
  else
     nullify(matCArray)
  end if
  
  !! Factor D
  if (lMatD) then
     allocate(matDArray(max(pfft_workSize,pfft_inLen(1))))
  else
     nullify(matDArray)
  end if


 

  ! inArray will store the Temperature and outArray Conductivity  
  call Grid_pfftMapToInput(iSrc , inArray)
  call Grid_pfftMapToInput(iFactorB,outArray)

  if(lConstFact)call Grid_pfftMapToInput(iFactorA,ConstFactArray1)
  if(lMatC)     call Grid_pfftMapToInput(iFactorC,  MatCArray)
  if(lMatD)     call Grid_pfftMapToInput(iFactorD,  MatDArray)   
  

  
  globalSize(:) = pfft_globalLen(:) 
  DELX = (gr_imax - gr_imin)/real(globalSize(1))
  DELY = (gr_jmax - gr_jmin)/real(globalSize(2))
  DELZ = (gr_kmax - gr_kmin)/real(globalSize(3))
  
  !PFFT_FORWARD

  if (lpass.LE.1 .OR. NDIM==1) then
     !Use Implicit Crank Nicholson
     NL = pfft_inLen(IAXIS) ! No of Points in X Direction.
     LL = pfft_inLen(JAXIS) ! No of Points in Y Direction.
     ML = pfft_inLen(KAXIS) ! No of Points in Z Direction.     
     
     call Diffuse_advance1D (NL,LL,ML,SWEEP_X,bcTypes(1:2),LbcValues(1:2),DELX,dt,theta,chi, &
          outArray, inArray, lConstFact, ConstFactArray1, lMatC, matCArray, lMatD, matDArray)
  end if

  if (NDIM ==  1) then
     outArray(:) = inArray(1:product(pfft_outLen)) !!reshape(tmp3DArray, (/product(pfft_outLen)/))
  end if
  
  if (NDIM >= 2) then     
     
     call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,inArray,&
          pfft_work1,pfft_inLen,pfft_midLen,&
          pfft_procGrid(JAXIS),pfft_comm(JAXIS))     
     inArray(:) = pfft_work1(:)     
     
     call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,outArray,&
          pfft_work1,pfft_inLen,pfft_midLen,&
          pfft_procGrid(JAXIS),pfft_comm(JAXIS))
     
     outArray(:) = pfft_work1(:)
    
     if(lConstFact) then
        call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,ConstFactArray1,&
             pfft_work1,pfft_inLen,pfft_midLen,&
             pfft_procGrid(JAXIS),pfft_comm(JAXIS))  
        ConstFactArray1(:) = pfft_work1(:)
     end if
     
     !! Transpose MatC I->J
     if(lMatC) then 
        call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,MatCArray,&
             pfft_work1,pfft_inLen,pfft_midLen,&
             pfft_procGrid(JAXIS),pfft_comm(JAXIS))  
        MatCArray(:) = pfft_work1(:)
     end if
     
     !! Transpose MatD I->J
     if(lMatD) then 
        call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,MatDArray,&
             pfft_work1,pfft_inLen,pfft_midLen,&
             pfft_procGrid(JAXIS),pfft_comm(JAXIS))  
        MatDArray(:) = pfft_work1(:)
     end if
     
     !! All arrays transformed.     

     if (lpass.LE.1 .OR. NDIM==2) then
        ! Transformed Domain
        NL = pfft_midLen(IAXIS) ! No of Points in Y Direction.
        LL = pfft_midLen(JAXIS) ! No of Points in Z Direction.
        ML = pfft_midLen(KAXIS) ! No of Points in X Direction.        
        
        ! The direction of sweep is still X as the entire domain has been transformed.


        call Diffuse_advance1D (NL,LL,ML,SWEEP_Y,bcTypes(3:4),LbcValues(3:4),DELY,dt, theta,chi, &
             outArray, inArray, lConstFact, ConstFactArray1, lMatC, matCArray, lMatD, matDArray)        
     end if   
  endif

  if (NDIM == 3) then
     if(lConstFact) then
        call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,ConstFactArray1,&
             pfft_work1,pfft_midLen,pfft_outLen,&
             pfft_procGrid(KAXIS),pfft_comm(KAXIS))
        ConstFactArray1(:) = pfft_work1(:)
     end if
     
     if(lMatC) then
        call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,MatCArray,&
             pfft_work1,pfft_midLen,pfft_outLen,&
             pfft_procGrid(KAXIS),pfft_comm(KAXIS))
        MatCArray(:) = pfft_work1(:)
     end if     

     if(lMatD) then

        call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,MatDArray,&
             pfft_work1,pfft_midLen,pfft_outLen,&
             pfft_procGrid(KAXIS),pfft_comm(KAXIS))
        MatDArray(:) = pfft_work1(:)
     end if

     call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,inArray,&
          pfft_work1,pfft_midLen,pfft_outLen,&
          pfft_procGrid(KAXIS),pfft_comm(KAXIS))
     inArray(:) = pfft_work1(:)

     call gr_pfftTranspose(PFFT_FORWARD,PFFT_PCLDATA_REAL,outArray,&
          pfft_work1,pfft_midLen,pfft_outLen,&
          pfft_procGrid(KAXIS),pfft_comm(KAXIS))

     outArray(:) = pfft_work1(:)
     
     ! Transformed Domain
     NL = pfft_outLen(IAXIS) ! No of Points in Z Direction.
     LL = pfft_outLen(JAXIS) ! No of Points in X Direction.
     ML = pfft_outLen(KAXIS) ! No of Points in Y Direction.


     call Diffuse_advance1D (NL,LL,ML,SWEEP_Z,bcTypes(5:6),LbcValues(5:6),DELZ,dt, theta,chi, &
          outArray, inArray, lConstFact, ConstFactArray1, lMatC, matCArray, lMatD, matDArray)          

  end if


  !PFFT_REVERSE
  if (NDIM == 3) then
     call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,inArray,&
          pfft_work2,pfft_outLen,pfft_midLen,&
          pfft_procGrid(KAXIS),pfft_comm(KAXIS))

     inArray(:) = pfft_work2(:)


     if (lpass.GT.1) then
        NL = pfft_midLen(IAXIS) ! No of Points in Y Direction.
        LL = pfft_midLen(JAXIS) ! No of Points in Z Direction.
        ML = pfft_midLen(KAXIS) ! No of Points in X Direction.

        
        call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,outArray,&
             pfft_work2,pfft_outLen,pfft_midLen,&
             pfft_procGrid(KAXIS),pfft_comm(KAXIS))
        
        outArray = pfft_work2
        
        if(lConstFact) then
           call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,ConstFactArray1,&
                pfft_work2,pfft_outLen,pfft_midLen,&
                pfft_procGrid(KAXIS),pfft_comm(KAXIS))

           ConstFactArray1 = pfft_work2
        endif
        
        if(lMatC) then
           call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,MatCArray,&
                pfft_work2,pfft_outLen,pfft_midLen,&
                pfft_procGrid(KAXIS),pfft_comm(KAXIS))
           
        end if

        if(lMatD) then                   

           call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,MatDArray,&
                pfft_work2,pfft_outLen,pfft_midLen,&
                pfft_procGrid(KAXIS),pfft_comm(KAXIS))
        end if
        

        call Diffuse_advance1D (NL,LL,ML,SWEEP_Y,bcTypes(3:4),LbcValues(3:4),DELY,dt, theta,chi, &
             outArray, inArray, lConstFact, ConstFactArray1, lMatC, matCArray, lMatD, matDArray)     
     end if
  end if
  
  if (NDIM >= 2) then         
     call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,inArray,&
          pfft_work1,pfft_midLen,pfft_inLen,&
          pfft_procGrid(JAXIS),pfft_comm(JAXIS))
     
     inArray = pfft_work1
     if (lpass.GT.1) then
        NL = pfft_inLen(IAXIS) ! No of Points in X Direction.
        LL = pfft_inLen(JAXIS) ! No of Points in Y Direction.
        ML = pfft_inLen(KAXIS) ! No of Points in Z Direction.         
        
        call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,outArray,&
             pfft_work2,pfft_midLen,pfft_inLen,&
             pfft_procGrid(JAXIS),pfft_comm(JAXIS))        
        outArray(:) = pfft_work2(:)       
        
        if(lConstFact) then
           call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,ConstFactArray1,&
                pfft_work2,pfft_midLen,pfft_inLen,&
                pfft_procGrid(JAXIS),pfft_comm(JAXIS))
           
           ConstFactArray1(:) = pfft_work2(:) 
        endif
        
        if(lMatC) then
           call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,MatCArray,&
                pfft_work2,pfft_midLen,pfft_inLen,&
                pfft_procGrid(JAXIS),pfft_comm(JAXIS))

           MatCArray = pfft_work2
           
        end if
        
        if(lMatD) then                   

           call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,MatDArray,&
                pfft_work2,pfft_midLen,pfft_inLen,&
                pfft_procGrid(JAXIS),pfft_comm(JAXIS))
           MatDArray = pfft_work2
           
        end if
        

        call Diffuse_advance1D (NL,LL,ML,SWEEP_X,bcTypes(1:2),LbcValues(1:2),DELX,dt,theta,chi, &
             outArray, inArray, lConstFact, ConstFactArray1, lMatC, matCArray, lMatD, matDArray)
     end if
     
     outArray(1:product(pfft_inLen)) = inArray(1:product(pfft_inLen))
     
  else
     ! 1D -- No Transpose
  end if
  
  outArray(1:inSize) = outArray(1:inSize)*scaleFact

  ! Map back to the non-uniform mesh
  call Grid_pfftMapFromOutput(iVar, outArray)
  
  deallocate(outArray) 
  deallocate(inArray)  

  if(lConstFact)deallocate (ConstFactArray1)
  if(lMatC)     deallocate (MatCArray)
  if(lMatD)     deallocate (MatDArray)
  
  if(.not.pfft_setupOnce) then
     call Grid_pfftFinalize()
  end if


  
  return
end subroutine Grid_advanceDiffusion
