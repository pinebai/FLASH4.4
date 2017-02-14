!!****h* source/physics/Diffuse/Diffuse_interface
!!
!! NAME
!!   Diffuse_interface
!!
!! SYNOPSIS
!!   use Diffuse_interface,ONLY: Diffuse_init, Diffuse_therm, Diffuse_visc, Diffuse_computeDt
!!
!! DESCRIPTION
!! This is the header file for the diffuse module that defines its
!! public interfaces.
!!***
Module Diffuse_interface
#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"

  interface
    subroutine Diffuse_computeDt (blockID, &
                                  xCenter,xLeft,xRight, dx, uxgrid, &
                                  yCenter,yLeft,yRight, dy, uygrid, &
                                  zCenter,zLeft,zRight, dz, uzgrid, &
                                  blkLimits,blkLimitsGC,        &
                                  solnData,   &
                                  dt_check, dt_minloc )
     

      integer, intent(IN) :: blockID
      integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
      real,INTENT(INOUT)    :: dt_check
      integer,INTENT(INOUT)    :: dt_minloc(5)
      real, pointer, dimension(:,:,:,:) :: solnData

      real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: &
          xCenter,xLeft,xRight
      real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: &
           yCenter,yLeft,yRight
      real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) ::&
           zCenter,zLeft,zRight
      real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: &
           dx, uxgrid
      real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: &
           dy, uygrid
      real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: &
           dz, uzgrid
    end subroutine Diffuse_computeDt
  end interface

  interface
    subroutine Diffuse_finalize()
    end subroutine Diffuse_finalize
  end interface

  interface
    subroutine Diffuse_init()
    end subroutine Diffuse_init
  end interface

  interface
    subroutine Diffuse_species(sweepDir, igeom, blockID,numCells,blkLimits,blkLimitsGC,&
                               leftCoords,rightCoords,&
                               temp_flx, temp_fly, temp_flz)

      integer, intent(IN) :: sweepDir, igeom, blockID, numCells
      integer, dimension(2,MDIM), intent (IN) :: blkLimitsGC, blkLimits
  
  
      real, intent(INOUT), DIMENSION(NFLUXES,                   &
                                     blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                                     blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                                     blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                                     temp_flx, temp_fly, temp_flz
  
      real, intent(IN), DIMENSION(numCells) :: leftCoords ,rightCoords
    end subroutine Diffuse_species
  end interface

  interface

    subroutine Diffuse_therm(sweepDir, igeom, blockID,numCells,blkLimits,blkLimitsGC,&
                             leftCoords,rightCoords,&
                             temp_flx, areaLeft)
      implicit none
      integer, intent(IN) :: sweepDir, igeom, blockID, numCells
      integer, dimension(2,MDIM), intent (IN) :: blkLimitsGC, blkLimits
  
#ifdef FIXEDBLOCKSIZE
      real, intent(INOUT), DIMENSION(NFLUXES,                 &
                               GRID_ILO_GC:GRID_IHI_GC,     &
                               GRID_JLO_GC:GRID_JHI_GC,     &
                               GRID_KLO_GC:GRID_KHI_GC) ::  &
                               temp_flx
      real, intent(IN), DIMENSION(                 &
                               GRID_ILO_GC:GRID_IHI_GC,     &
                               GRID_JLO_GC:GRID_JHI_GC,     &
                               GRID_KLO_GC:GRID_KHI_GC) ::  &
                               areaLeft
      real, intent(IN), DIMENSION(MAXCELLS) :: leftCoords ,rightCoords
#else
      real, intent(INOUT), DIMENSION(NFLUXES,                   &
                                     blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                                     blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                                     blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                                     temp_flx
      real, intent(IN), DIMENSION(                   &
                                     blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                                     blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                                     blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                                     areaLeft
      real, intent(IN), DIMENSION(numCells) :: leftCoords ,rightCoords
#endif
    end subroutine Diffuse_therm
  end interface

  interface
     subroutine Diffuse_visc(sweepDir, igeom, blockID,numCells,blkLimits,blkLimitsGC,&
                             leftCoords,rightCoords,&
                             temp_flx, areaLeft, secondCoord, thirdCoord)
       implicit none
       integer, intent(IN) :: sweepDir, igeom, blockID, numCells
       integer, dimension(2,MDIM), intent (IN) :: blkLimitsGC, blkLimits
  
#ifdef FIXEDBLOCKSIZE
       real, intent(INOUT), DIMENSION(NFLUXES,              &
                               GRID_ILO_GC:GRID_IHI_GC,     &
                               GRID_JLO_GC:GRID_JHI_GC,     &
                               GRID_KLO_GC:GRID_KHI_GC) ::  &
                               temp_flx
       real, intent(IN), DIMENSION(                         &
                               GRID_ILO_GC:GRID_IHI_GC,     &
                               GRID_JLO_GC:GRID_JHI_GC,     &
                               GRID_KLO_GC:GRID_KHI_GC) ::  &
                               areaLeft
       real, intent(IN), DIMENSION(MAXCELLS) :: leftCoords ,rightCoords, &
            secondCoord, thirdCoord
#else
       real, intent(INOUT), DIMENSION(NFLUXES,                   &
                               blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                               temp_flx
       real, intent(IN), DIMENSION(                   &
                                   blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
                                   blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
                                   blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
                                   areaLeft
       real, intent(IN), DIMENSION(numCells) :: leftCoords ,rightCoords, &
            secondCoord, thirdCoord
#endif  
    end subroutine Diffuse_visc
  end interface

  interface
     subroutine Diffuse(blockCount,blockList,dt,pass)
       implicit none
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       real,intent(in) :: dt
       integer, OPTIONAL, intent(IN):: pass
     end subroutine Diffuse
  end interface

  interface Diffuse_advance1D
     subroutine Diffuse_advance1D (nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray,ConstFact,MatC, MatD)
       implicit none
       integer,                intent(IN)    :: nx, ny, nz
       integer,                intent(IN)    :: iDirection
       integer, dimension (2), intent(IN)    :: bcType
       real   , dimension (2), intent(IN)    :: bcValues  
       real   ,                intent(IN)    :: DEL
       real   ,                intent(IN)    :: dt
       real   ,                intent(IN)    :: theta
       real   ,                intent(IN)    :: chi
       real   ,                intent(IN)    :: cond(nx,ny,nz)
       real   ,                intent(INOUT) :: soln(nx,ny,nz)
       logical,                intent(IN)    :: useCFArray
       real   ,OPTIONAL,       pointer       :: ConstFact(:,:,:) !ConstFact(nx,ny,nz)
       real   ,OPTIONAL,       pointer       :: MatC(:,:,:)   !ConstFact(nx,ny,nz) 
       real   ,OPTIONAL,       pointer       :: MatD(:,:,:)   !ConstFact(nx,ny,nz)  
     end subroutine Diffuse_advance1D
     subroutine Diffuse_advance1D_1DArr (nx,ny,nz,iDirection,bcType,bcValues,DEL,dt,theta,chi,cond,soln,useCFArray,ConstFact, &
          useMatCArray,MatC, useMatDArray,MatD)
       implicit none
       integer,                intent(IN)    :: nx, ny, nz
       integer,                intent(IN)    :: iDirection
       integer, dimension (2), intent(IN)    :: bcType
       real   , dimension (2), intent(IN)    :: bcValues  
       real   ,                intent(IN)    :: DEL
       real   ,                intent(IN)    :: dt
       real   ,                intent(IN)    :: theta
       real   ,                intent(IN)    :: chi
       real   ,                intent(IN)    :: cond(nx*ny*nz)
       real   ,                intent(INOUT) :: soln(nx*ny*nz)
       logical,                intent(IN)    :: useCFArray
       real   ,                pointer       :: ConstFact(:) !ConstFact(nx*ny*nz)
       logical,                intent(IN)    :: useMatCArray
       real   ,                pointer       :: MatC(:)   !matC(nx,ny,nz) 
       logical,                intent(IN)    :: useMatDArray
       real   ,                pointer       :: MatD(:)   !matD(nx,ny,nz)  
     end subroutine Diffuse_advance1D_1DArr
  end interface

  interface 
     subroutine Diffuse_solveScalar(iVar, iFactorB, iFactorA, &
          bcTypes, bcValues, dt, scaleFact, chi, theta, pass, &
          blockCount, blockList, iFactorC, iFactorD)
       implicit none
       integer, intent(IN)           :: iVar                                              
       integer, intent(IN)           :: iFactorB
       integer, intent(IN)           :: iFactorA
       integer, intent(IN), optional :: iFactorC
       integer, intent(IN), optional :: iFactorD
       integer, intent(IN)           :: bcTypes(6)
       real,    intent(IN)           :: bcValues(2,6)
       real,    intent(IN)           :: dt
       real,    intent(IN)           :: scaleFact
       real,    intent(IN)           :: chi
       real,    intent(IN)           :: theta
       integer, OPTIONAL, intent(IN) :: pass
       integer,           intent(IN) :: blockCount
       integer, dimension(blockCount),intent(IN) :: blockList
       
     end subroutine Diffuse_solveScalar
     subroutine Diffuse_solveCoupledScalar (unkVarsDesc, xtraVarDesc, diffCoeffDesc, absorpCoeffDesc, bcTypes, bcValues, &
          dt, scaleFact, theta, pass, blockCount, blockList, iFactorD, iFactorA)
       implicit none
       integer, dimension(VARDESC_SIZE), intent(IN):: unkVarsDesc
       integer, dimension(VARDESC_SIZE), intent(IN):: xtraVarDesc
       integer, dimension(VARDESC_SIZE), intent(IN):: diffCoeffDesc
       integer, dimension(VARDESC_SIZE), OPTIONAL,intent(IN):: absorpCoeffDesc
       integer, intent(IN):: bcTypes(6)
       real, intent(IN),OPTIONAL:: bcValues(2,6)
       real, intent(IN):: dt
       real, intent(IN):: scaleFact
       real, intent(IN):: theta
       integer, OPTIONAL,intent(IN):: pass
       integer,           intent(IN) :: blockCount
       integer, dimension(blockCount),intent(IN) :: blockList
       integer, dimension(VARDESC_SIZE), OPTIONAL,intent(IN):: iFactorD
       integer, dimension(VARDESC_SIZE), OPTIONAL, intent(IN):: iFactorA
     end subroutine Diffuse_solveCoupledScalar
  end interface
  
  interface
     subroutine Diffuse_fluxLimiter(idcoef, ifunc, ifl, mode, blkcnt, blklst)
       implicit none
       
       ! Arguments:
       integer, intent(in) :: idcoef
       integer, intent(in) :: ifunc
       integer, intent(in) :: ifl
       integer, intent(in) :: mode
       integer, intent(IN) :: blkcnt
       integer, intent(IN) :: blklst(blkcnt)
       
     end subroutine Diffuse_fluxLimiter

     subroutine Diffuse_computeFluxLimiter(idcoef, ifunc, iDenForIfunc, ifl, iflOut, ieddi3,&
          mode, solnData, lbUI,lbUJ,lbUK,&
          blockID, gcLayers)
       implicit none
       integer, intent(in) :: idcoef
       integer, intent(in) :: ifunc
       integer, intent(in) :: iDenForIfunc
       integer, intent(in) :: ifl
       integer, intent(in) :: iflOut
       integer, intent(in) :: ieddi3
       integer, intent(in) :: mode
       integer, VALUE_INTENT(IN) :: lbUI,lbUJ,lbUK
       real,    intent(INOUT) :: solnData(:,lbUI:,lbUJ:,lbUK:)
       integer, intent(IN) :: blockID
       integer, intent(IN),OPTIONAL :: gcLayers
     end subroutine Diffuse_computeFluxLimiter
  end interface
    
  interface
     subroutine Diffuse_setContextInfo(group,component)
       implicit none
       integer,intent(in),OPTIONAL :: component, group
     end subroutine Diffuse_setContextInfo
  end interface

end Module Diffuse_interface
