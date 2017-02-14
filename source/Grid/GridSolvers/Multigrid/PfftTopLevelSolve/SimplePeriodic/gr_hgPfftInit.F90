!!****if* source/Grid/GridSolvers/Multigrid/PfftTopLevelSolve/SimplePeriodic/gr_hgPfftInit
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
           call Driver_abortFlash('BCs unsupported')        
     end select
  enddo


  ! Define Solveflag
  ! solveflag  1 => periodic in x, Neuman conditions in y, z  (BLKTRI)
  ! solveflag  2 => periodic in x, z and Neuman in y  (TRIDIAG)
  ! solveflag  3 => periodic in x, y and z
  ! Use BC_type to define the solver, if BC arrangement not supported stop.
  if (NDIM .eq. 1) then

     if(gr_hgbcTypes(1) .eq. MG_BND_PERIODIC ) then    ! in x periodic
        
         solveflag = 3
!         gr_hgBndConditions = MG_BND_PERIODIC          ! Defined to do mean substractions in hg Multigrid  

     else

        if (gr_meshMe .eq. 0) then
           write(*,*) 'Error: Poisson Solution Boundary Conditions Not Supported'
           write(*,*) 'gr_hgbcTypes x =',gr_hgbcTypes(1),gr_hgbcTypes(2)
           call Driver_abortFlash('BCs unsupported')
        endif

     endif

  elseif (NDIM .eq. 2) then


     if(gr_hgbcTypes(1) .eq. MG_BND_PERIODIC .and. &   ! in x periodic
        gr_hgbcTypes(3) .eq. MG_BND_PERIODIC ) then    ! in y periodic
        
         solveflag = 3
!         gr_hgBndConditions = MG_BND_PERIODIC          ! Defined to do mean substractions in hg Multigrid  

     else

        if (gr_meshMe .eq. 0) then
           write(*,*) 'Error: Poisson Solution Boundary Conditions Not Supported'
           write(*,*) 'gr_hgbcTypes x =',gr_hgbcTypes(1),gr_hgbcTypes(2)
           write(*,*) 'gr_hgbcTypes y =',gr_hgbcTypes(3),gr_hgbcTypes(4)
           call Driver_abortFlash('BCs unsupported')
        endif

     endif

  else

     if(gr_hgbcTypes(1) .eq. MG_BND_PERIODIC .and. &   ! in x periodic
            gr_hgbcTypes(3) .eq. MG_BND_PERIODIC .and. &   ! in y periodic
            gr_hgbcTypes(5) .eq. MG_BND_PERIODIC ) then    ! in z periodic

         solveflag = 3
!         gr_hgBndConditions = MG_BND_PERIODIC

     else

        if (gr_meshMe .eq. 0) then
           write(*,*) 'Error: Poisson Solution Boundary Conditions Not Supported'
           write(*,*) 'gr_hgbcTypes x =',gr_hgbcTypes(1),gr_hgbcTypes(2)
           write(*,*) 'gr_hgbcTypes y =',gr_hgbcTypes(3),gr_hgbcTypes(4)
           write(*,*) 'gr_hgbcTypes z =',gr_hgbcTypes(5),gr_hgbcTypes(6)
           call Driver_abortFlash('BCs unsupported')
        endif

     endif

  endif


  !We will refer to these module level variables in gr_hgPfftInitGrid:
  gr_hgPfftSolveFlag = solveFlag
  gr_hgPfftLastMappedLevel = NONEXISTENT
  
end subroutine gr_hgPfftInit
