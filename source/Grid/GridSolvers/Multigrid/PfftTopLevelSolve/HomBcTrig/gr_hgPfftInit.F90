!!****if* source/Grid/GridSolvers/Multigrid/PfftTopLevelSolve/HomBcTrig/gr_hgPfftInit
!!
!! NAME
!!
!!  gr_hgPfftInit
!!
!! SYNOPSIS
!!
!!  gr_hgPfftInit()
!!
!! DESCRIPTION
!!
!!  This routine initializes the data necessary for the Multigrid PFFT extensions.
!!
!! ARGUMENTS
!!
!!***

subroutine gr_hgPfftInit()

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe, gr_domainBC
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get, &
       RuntimeParameters_mapStrToInt
  use gr_hgData , ONLY : gr_hgBndTypes
  use gr_hgPfftData, ONLY : gr_hgPfftSolveFlag, gr_hgPfftLastMappedLevel, &
       gr_hgPfftMaxDirectSolveLevel, gr_hgPfftBcTypes, gr_hgbcTypes

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multigrid.h"

  !--------------------------------------------------------------------------
  integer :: solveFlag, eachBoundary
  character(len=MAX_STRING_LENGTH), dimension(2*MDIM) :: bcTypeStr

  !The hybrid solver performs a direct FFT solve on the finest
  !fully refined level by default.  To alter this behavior set the
  !following runtime parameter to a more coarse level.
  call RuntimeParameters_get("maxDirectSolveLevel", &
       gr_hgPfftMaxDirectSolveLevel)

  !Initialise to -1 to help us spot errors.
  gr_hgbcTypes = -1

  ! Retrieve Poisson solution Boundary Conditions for each face: 
  call RuntimeParameters_get("xl_mg_boundary_type", bcTypeStr(1))
  call RuntimeParameters_get("xr_mg_boundary_type", bcTypeStr(2))
  if (NDIM >= 2) then
     call RuntimeParameters_get("yl_mg_boundary_type", bcTypeStr(3))
     call RuntimeParameters_get("yr_mg_boundary_type", bcTypeStr(4))
  endif
  if (NDIM == 3) then
     call RuntimeParameters_get("zl_mg_boundary_type", bcTypeStr(5))
     call RuntimeParameters_get("zr_mg_boundary_type", bcTypeStr(6))
  endif


  do eachBoundary = 1, 2*NDIM

     call RuntimeParameters_mapStrToInt(bcTypeStr(eachBoundary), gr_hgPfftBcTypes(eachBoundary))

     select case(gr_hgPfftBcTypes(eachBoundary))
     case (PERIODIC) 
        gr_hgbcTypes(eachBoundary)  = MG_BND_PERIODIC
     case (OUTFLOW)
        gr_hgbcTypes(eachBoundary)  = MG_BND_NEUMANN
     case (DIRICHLET)
        gr_hgbcTypes(eachBoundary)  = MG_BND_DIRICHLET
     case default
           write(*,*) 'In gr_hgPfftInit: Poisson Problem Boundary Condition not supported'
           write(*,*) 'gr_hgPfftBcTypes:',gr_hgPfftBcTypes
           write(*,*) 'gr_hgbcTypes:',gr_hgbcTypes
           call Driver_abortFlash('BCs unsupported')        
     end select
  enddo


  ! Define Solveflag


  !We will refer to these module level variables in gr_hgPfftInitGrid:
  gr_hgPfftSolveFlag = -1
  gr_hgPfftLastMappedLevel = NONEXISTENT
  
end subroutine gr_hgPfftInit
