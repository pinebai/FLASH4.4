!!****if* source/Grid/GridSolvers/Pfft/gr_pfftInit
!!
!! NAME
!!
!!  gr_pfftInit
!!
!! 
!! SYNOPSIS
!!
!!  call gr_pfftInit()
!!
!!
!! DESCRIPTION
!!
!!  Initialize the Pfft Poisson solver.  Read in any of the
!!  runtime parameters for this solver, and other one-time
!!  initializations.
!!
!!  The solver common data for this implementation
!!  is stored in the gr_pfftData FORTRAN module.
!!
!! PARAMETERS
!!
!!  pfft_setupOnce   LOGICAL
!!
!!
!!***

subroutine gr_pfftInit()
#include "constants.h"
#include "Pfft.h"
#include "Flash.h"

  use Grid_data, ONLY : gr_meshMe, gr_geometry, gr_domainBC, gr_oneRefLev
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       GRID_PDE_BND_DIRICHLET, &
       Grid_getGlobalIndexLimits, Grid_pfftInit
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use gr_pfftData, ONLY : pfft_setupOnce, gr_pfftDiffOpDiscretize
  use Logfile_interface, ONLY : Logfile_stampMessage
  use gr_pfftInterface, ONLY : gr_pfftSpecifyTransform
#if defined(FLASH_GRID_UG)
  use Grid_data, ONLY : gr_axisNumProcs
#elif defined(FLASH_GRID_PARAMESH)
  use Grid_data, ONLY : gr_nBlockX, gr_nBlockY, gr_nBlockZ
  use tree, ONLY : lrefine_min, lrefine_max
#endif

  implicit none

  integer, dimension(MDIM):: transformType, globalLen,localLen
  integer, dimension(0:MDIM):: baseDatType
  integer, dimension(2*MDIM) :: bcTypes
  integer :: nprocj,nprock
  logical :: needMap

  call RuntimeParameters_get("pfft_setupOnce",pfft_setupOnce)
  call RuntimeParameters_get("gr_pfftDiffOpDiscretize",gr_pfftDiffOpDiscretize)
  call Grid_getGlobalIndexLimits(globalLen)


  ! If we have a uniform grid then there is no reason for 
  ! pfft_setupOnce to be false.  Force it to be true, but let the user know.
#if defined(FLASH_GRID_UG)
  if (.not. pfft_setupOnce) then
     pfft_setupOnce = .true.
     if (gr_meshMe == MASTER_PE) then
        print *, "[gr_pfftInit]: UG mode... forcing pfft_setupOnce = .true."
        call Logfile_stampMessage( & 
             "[gr_pfftInit]: UG mode... forcing pfft_setupOnce = .true.")
     end if
  end if
#endif
  
#ifdef DEBUG_PARTICLES
  print*,'pfft_setupOnce is', pfft_setupOnce
#endif

  if(pfft_setupOnce)then
     bcTypes =  reshape(gr_domainBC, (/ 2*MDIM /) )
     where (bcTypes == PERIODIC)
        bcTypes = GRID_PDE_BND_PERIODIC
     elsewhere (bcTypes == REFLECTING)
        bcTypes = GRID_PDE_BND_DIRICHLET
     elsewhere (bcTypes == DIRICHLET)
        bcTypes = GRID_PDE_BND_DIRICHLET
     elsewhere (bcTypes == OUTFLOW)
        bcTypes = GRID_PDE_BND_NEUMANN
     elsewhere (bcTypes == SLIP_INS)
        bcTypes = GRID_PDE_BND_NEUMANN
     elsewhere (bcTypes == NOSLIP_INS)
        bcTypes = GRID_PDE_BND_NEUMANN
     elsewhere (bcTypes == INFLOW_INS)
        bcTypes = GRID_PDE_BND_NEUMANN
     elsewhere (bcTypes == MOVLID_INS)
        bcTypes = GRID_PDE_BND_NEUMANN
     elsewhere (bcTypes == NEUMANN_INS)
        bcTypes = GRID_PDE_BND_NEUMANN
     elsewhere (bcTypes == OUTFLOW_INS)
        bcTypes = GRID_PDE_BND_NEUMANN
     end where
     call gr_pfftSpecifyTransform(transformType, baseDatType, bcTypes)

#if defined(FLASH_GRID_UG)

     needMap=(gr_axisNumProcs(IAXIS)>1)
     nprocj=gr_axisNumProcs(JAXIS)
     nprock=gr_axisNumProcs(KAXIS)

#elif defined(FLASH_GRID_PARAMESH)

     !If lrefine_min /= lrefine_max then the PARAMESH grid is likely to 
     !change.  As such we must give PFFT the chance to setup multiple PFFT grids.
     if (lrefine_min /= lrefine_max) then
        pfft_setupOnce = .false.
        if (gr_meshMe == MASTER_PE) then
           print *, "[gr_pfftInit]: lrefine_min /= lrefine_max... " // &
                "forcing pfft_setupOnce = .false."
           call Logfile_stampMessage( "[gr_pfftInit]: " // &
                "lrefine_min /= lrefine_max... forcing pfft_setupOnce = .false.")
        end if
     end if
     needMap = (NDIM > 1)   ! For 1D, PFFT requires to be run on only 1 proc anyway.

#endif

     if(needMap) then
        call Grid_pfftInit(NDIM,needMap,globalLen,&
             localLen,transformType=transformType, baseDatType=baseDatType)
     else
        call Grid_pfftInit(NDIM,needMap,globalLen,&
             localLen,transformType=transformType, baseDatType=baseDatType,&
             jProcs=nprocj,&
             kProcs=nprock)
     end if
  end if

end subroutine gr_pfftInit
