!!****ih* source/Grid/GridSolvers/HYPRE/gr_hypreLocalInterface
!!
!! This is a module file for the HYPRE solver in FLASH that defines
!! additional interfaces private to the GridSolvers/HYPRE implementations.
!!
!! NOTES
!!
!!  Adding explicit interfaces here is being tried as an alternative to
!!  writiting executable FORTRAN wrappers for each additional HYPRE routine
!!  we want to call from FLASH.
!!
!! SEE ALSO
!!
!!  gr_hypreF90CAdapters.c
!!***

#include "constants.h"

Module gr_hypreLocalInterface

#if 0
  ! Maybe one day... for now, not all compilers support this.
  interface
     subroutine hypre_sstructinnerprod(fx, fy, fresult, ierr) &
          bind(C,name="HYPRE_SStructInnerProd")
       implicit none
       integer*8,intent(IN)  :: fx, fy
       real,     intent(OUT) :: fresult
       integer,  intent(OUT) :: ierr
    end subroutine hypre_sstructinnerprod

     integer function hypre_pcggetconverged(solver, converged) &
          bind(C,name="HYPRE_PCGGetConverged")
       implicit none
       integer*8,VALUE :: solver ! intent(IN)
       integer,intent(OUT) :: converged
     end function hypre_pcggetconverged

     subroutine hypre_describeerror(hypreErr, description) &
          bind(C,name="HYPRE_DescribeError")
       implicit none
       integer,VALUE :: hypreErr ! intent(IN)
       character(len=1),dimension(MAX_STRING_LENGTH),intent(OUT) :: description
     end subroutine hypre_describeerror

  end interface

#else
  ! fallback... wrapper implementations are in gr_hypreF90CAdapters.c

  interface
     subroutine hypre_sstructinnerprod(fx, fy, fresult, ierr)
       implicit none
       integer*8,intent(IN)  :: fx, fy
       real,     intent(OUT) :: fresult
       integer,  intent(OUT) :: ierr
    end subroutine hypre_sstructinnerprod

     subroutine hypre_pcggetconverged(solver, converged, ierr)
       implicit none
       integer*8,intent(IN) :: solver
       integer,intent(OUT) :: converged
       integer,  intent(OUT) :: ierr
     end subroutine hypre_pcggetconverged

     subroutine hypre_describeerror(hypreErr, description)
       implicit none
       integer,intent(IN) :: hypreErr
       character(len=MAX_STRING_LENGTH),intent(OUT) :: description
     end subroutine hypre_describeerror

  end interface

#endif

  interface
     subroutine gr_hypreAddGraph (blockHandle, blockID, blkPartNo, direction, datasize, CornerID, blkStride, &
          firstHypreVar, numVars)
       implicit none
       integer, intent(IN) :: blockHandle
       integer, intent(IN) :: direction
       integer, intent(IN) :: blockID   
       integer, intent(IN) :: blkPartNo
       integer, intent(IN) :: datasize(MDIM)
       integer, intent(IN) :: CornerID(MDIM)
       integer, intent(IN) :: blkStride(MDIM)
       integer, intent(IN),OPTIONAL :: firstHypreVar, numVars
     end subroutine gr_hypreAddGraph
  end interface

#include "FortranLangFeatures.fh"

  interface
     subroutine gr_hypreGetFaceBFcB (direction, blkLimits, blkLimitsGC, facBptr, flux, iVar)
       implicit none
       integer, intent(IN) :: direction
       integer, intent(IN) :: blkLimits (2,MDIM) 
       integer, intent(IN) :: blkLimitsGC (2,MDIM)
       real, POINTER_INTENT_IN :: facBptr(:,:,:)   
       real, intent(INOUT) :: flux(:,:,:,:)
       integer, intent(IN) :: iVar
     end subroutine gr_hypreGetFaceBFcB
  end interface

  interface
     subroutine gr_hypreExchangeFacB (iFactorB, blockCount, blockList)
       implicit none
       integer,intent(IN) :: iFactorB
       integer,intent(IN) :: blockCount
       integer, dimension(blockCount),intent(IN):: blockList
     end subroutine gr_hypreExchangeFacB
  end interface

  interface
     subroutine gr_hypreExchangeFacBFcB (iFactorB, blockCount, blockList)
       implicit none
       integer,intent(IN) :: iFactorB
       integer,intent(IN) :: blockCount
       integer, dimension(blockCount),intent(IN):: blockList
     end subroutine gr_hypreExchangeFacBFcB
  end interface

  interface
     subroutine gr_hypreMultiExchangeFacB (unkVarsDesc, firstHypreVar,diffCoeffDesc, blockCount, blockList)
       implicit none
       integer,intent(IN) :: unkVarsDesc(VARDESC_SIZE)
       integer,intent(IN) :: firstHypreVar
       integer,intent(IN) :: diffCoeffDesc(VARDESC_SIZE)
       integer,intent(IN) :: blockCount
       integer, dimension(blockCount),intent(IN):: blockList
     end subroutine gr_hypreMultiExchangeFacB
  end interface

  interface
     subroutine gr_hypreMultiAddToMatrix(numVars, firstHypreVar, iFactorB, &
          absorpCoeffDesc, &
          emissCoeffDesc,emissTermDesc, &
          unkVarsDesc, &
          iFactorA, bcTypes, bcValues, dt, &
          alpha, blockCount, blockList, JacobiMatrix)
       implicit none
       integer, intent(IN) :: numVars
       integer, intent(IN) :: firstHypreVar
       integer, intent(IN) :: iFactorB
       integer,intent(in),dimension(:) :: absorpCoeffDesc
       integer,intent(in),dimension(:) :: emissCoeffDesc, emissTermDesc
       integer,intent(in),dimension(:) :: unkVarsDesc
       integer, intent(IN) :: iFactorA
       integer, intent(IN) :: bcTypes(6)
       real,    intent(IN) :: bcValues(:,:,:)
       real,    intent(IN) :: dt
       real,    intent(IN) :: alpha
       integer, intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       logical, intent(IN) :: JacobiMatrix
     end subroutine gr_hypreMultiAddToMatrix
  end interface

#include "Flash.h"  
  interface
     subroutine gr_hypreApplyBcToFace(blkLimits,blkLimitsGC,part,var,iFactorB,bcType,direction, &
          bcValue, dt, theta, del, Lower, scalefactor, faceArea, solnVec, blockID, &
          facBptrX,facBptrY,facBptrZ)
       implicit none
       integer, intent(IN) :: blkLimits (2,MDIM) 
       integer, intent(IN) :: blkLimitsGC (2,MDIM)
       integer, intent(IN) :: part
       integer, intent(IN) :: var
       integer, intent(IN) :: iFactorB
       integer, intent(IN) :: bcType
       integer, intent(IN) :: direction
       real,    intent(IN) :: bcValue(2)
       real,    intent(IN) :: dt
       real,    intent(IN) :: theta
       real,    intent(IN) :: del
       integer, intent(IN) :: Lower(MDIM)
       real,    intent(IN) :: scalefactor
       real,    intent(IN) :: faceArea(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
       real,    intent(IN) :: solnVec(NUNK_VARS, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))   
       integer, intent(IN) :: blockID
       real, POINTER_INTENT_IN, OPTIONAL :: facBptrX(:,:,:,:),facBptrY(:,:,:,:),facBptrZ(:,:,:,:)
     end subroutine gr_hypreApplyBcToFace
  end interface

  interface
     subroutine gr_hypreMultiSetIniGuess (varDesc, firstHypreVar, blockCount, blockList)
       implicit none
       integer,intent(IN) :: varDesc(VARDESC_SIZE)
       integer,intent(IN) :: firstHypreVar
       integer,intent(IN) :: blockCount
       integer,intent(IN) :: blockList (blockCount)
     end subroutine gr_hypreMultiSetIniGuess
  end interface

end Module gr_hypreLocalInterface


