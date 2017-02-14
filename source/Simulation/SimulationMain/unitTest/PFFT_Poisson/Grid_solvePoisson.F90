!!****if* source/Simulation/SimulationMain/unitTest/PFFT_Poisson/Grid_solvePoisson
!!
!! NAME
!!
!!  Grid_solvePoisson
!!
!!
!! SYNOPSIS
!!
!!   Grid_solvePoisson(integer(IN) :: iSoln,
!!                     integer(IN) :: iSrc, 
!!                  integer(6)(IN) :: bcTypes,
!!                   real(2,6)(IN) :: bcValues,
!!                     real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This module implements the 
!!   fft based method for periodic and dirichlet problems
!!   on a uniform grid.  
!!   Isolated problems are not supported
!!
!!
!! ARGUMENTS
!!
!!  iSoln -  index to variable containing potential
!!  iSrc - index to variable containing density
!!  bcTypes - boundary types along various faces,
!!          valid values are: (although only some are implemented)
!!          PERIODIC -- supported
!!          DIRICHLET
!!          OUTFLOW (Neuman)
!!          HYDROSTATIC (Dirichlet value)
!!          ISOLATED
!!  bcValues - the values to boundary conditions, currently not used
!!  poisfact -  factor to be used in calculation
!!
!!***

subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)

  use gr_pfftData, ONLY : pfft_inLen,pfft_outLen,pfft_setupOnce,pfft_usableProc

  use Grid_interface, ONLY : Grid_pfftMapToInput,Grid_getGlobalIndexLimits,&
       Grid_getBlkIndexLimits,Grid_pfftInit, Grid_pfft,Grid_pfftMapFromOutput, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_pfftFinalize
  use gr_pfftInterface, ONLY : gr_pfftDerivs, gr_pfftSpecifyTransform
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"  

  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes(2*MDIM)
  real, intent(in)       :: bcValues(2,2*MDIM)
  real, intent(inout)    :: poisfact

  !--------------------------------------------------------------------------
  real, allocatable, dimension(:) :: inArray,tranArray,outArray
  integer, dimension(MDIM) :: localSize,globalSize,transformType
  integer, dimension(0:MDIM) :: baseDatType
  integer :: inSize,tranSize
  logical :: needMap

  if(.not.pfft_setupOnce) then
     needMap=.true.
     call Grid_getGlobalIndexLimits(globalSize)
     call gr_pfftSpecifyTransform(transformType, baseDatType, bcTypes)
     call Grid_pfftInit(NDIM,needMap,globalSize, &
          localSize,transformType, baseDatType)
  end if

  !Important.  Tests that this processor should be doing work
  if(.not.pfft_usableProc) return   

  inSize=pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
  tranSize=2*pfft_outLen(IAXIS)*pfft_outLen(JAXIS)*pfft_outLen(KAXIS)

  allocate(inArray(inSize+2))
  allocate(outArray(inSize+2))
  allocate(tranArray(tranSize))


  !! Here's the real work of the fft
  ! Converts to uniform mesh (on output, inArray contains uniformly mapped density)
  call Grid_pfftMapToInput(iSrc,inArray)  
  ! Forward transform of density 
  call Grid_pfft(PFFT_FORWARD,inArray,tranArray)
  ! Calculates the transform of iSoln = GPOT
  !  which is the transform of the delSquared(u) = rho in Poisson equation
  call gr_pfftDerivs(tranArray)
  ! Inverse transform of GPOT
  call Grid_pfft(PFFT_INVERSE,tranArray,outArray)
  ! Now multiply by the poisson factor
   outArray(1:inSize) = outArray(1:inSize)*poisfact

  ! Map back to the non-uniform mesh
  call Grid_pfftMapFromOutput(iSoln,outArray)


  deallocate(inArray)
  deallocate(tranArray)
  deallocate(outArray)
  if(.not.pfft_setupOnce) then
     call Grid_pfftFinalize()
  end if

  return
end subroutine Grid_solvePoisson
