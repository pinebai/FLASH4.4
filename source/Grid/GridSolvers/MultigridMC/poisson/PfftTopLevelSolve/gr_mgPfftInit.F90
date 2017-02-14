!!****if* source/Grid/GridSolvers/MultigridMC/poisson/PfftTopLevelSolve/gr_mgPfftInit
!!
!! NAME
!!
!!  gr_mgPfftInit
!!
!! SYNOPSIS
!!
!!  gr_mgPfftInit()
!!
!! DESCRIPTION
!!
!!  This routine initializes the data necessary for the Multigrid PFFT extensions.
!!  Same as hg solver Multigrid.
!!
!! ARGUMENTS
!!
!!***

subroutine gr_mgPfftInit()

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_meshMe, gr_domainBC
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get, &
       RuntimeParameters_mapStrToInt

  use gr_mgData , ONLY : gr_mgBndTypes

  use gr_mgPfftData, ONLY : gr_mgPfftSolveFlag, gr_mgPfftLastMappedLevel, &
       gr_mgPfftMaxDirectSolveLevel, gr_mgPfftBcTypes, gr_mgbcTypes

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
       gr_mgPfftMaxDirectSolveLevel)
  !Initialise to -1 to help us spot errors.
  gr_mgbcTypes = -1

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

     call RuntimeParameters_mapStrToInt(bcTypeStr(eachBoundary), gr_mgPfftBcTypes(eachBoundary))

     select case(gr_mgPfftBcTypes(eachBoundary))
     case (PERIODIC) 
        gr_mgbcTypes(eachBoundary)  = MG_BND_PERIODIC
     case (OUTFLOW)
        gr_mgbcTypes(eachBoundary)  = MG_BND_NEUMANN
     case (DIRICHLET)
        gr_mgbcTypes(eachBoundary)  = MG_BND_DIRICHLET
     case default
           write(*,*) 'In gr_mgPfftInit: Poisson Problem Boundary Condition not supported'
           call Driver_abortFlash('BCs unsupported')        
     end select
  enddo


  ! Define Solveflag
  ! solveflag  1 => periodic in x, Neuman conditions in y, z  (BLKTRI)
  ! solveflag  2 => periodic in x, z and Neuman in y  (TRIDIAG)
  ! solveflag  3 => periodic in x, y and z
  ! Use BC_type to define the solver, if BC arrangement not supported stop.
  if (NDIM .eq. 1) then

        if (gr_meshMe .eq. 0) then
           write(*,*) 'Error: 1D Poisson Solution Not Supported'
           call Driver_abortFlash('BCs unsupported')
        endif

  elseif (NDIM .eq. 2) then


     if(gr_mgbcTypes(1) .eq. MG_BND_PERIODIC .and. &   ! in x periodic
        gr_mgbcTypes(3) .eq. MG_BND_PERIODIC ) then    ! in y periodic
        
         solveflag = 3
!         gr_mgBndConditions = MG_BND_PERIODIC         ! Defined to do mean substractions in mg Multigrid  

     else

         solveflag = -1

     endif

  else

     if( gr_mgbcTypes(1) .eq. MG_BND_PERIODIC  .and.  &   ! in x periodic
         gr_mgbcTypes(3) .eq. MG_BND_NEUMANN   .and.  &   ! in y Neuman
         gr_mgbcTypes(5) .eq. MG_BND_NEUMANN  ) then      ! in z Neuman

         solveflag = 1   
!         gr_mgBndConditions = MG_BND_NEUMANN


     elseif(gr_mgbcTypes(1) .eq. MG_BND_PERIODIC .and. &  ! in x periodic
            gr_mgbcTypes(3) .eq. MG_BND_PERIODIC .and. &  ! in y periodic
            gr_mgbcTypes(5) .eq. MG_BND_NEUMANN ) then    ! in z Neuman

         solveflag = 2
!         gr_mgBndConditions = MG_BND_NEUMANN

     elseif(gr_mgbcTypes(1) .eq. MG_BND_PERIODIC .and. &   ! in x periodic
            gr_mgbcTypes(3) .eq. MG_BND_PERIODIC .and. &   ! in y periodic
            gr_mgbcTypes(5) .eq. MG_BND_PERIODIC ) then    ! in z periodic

         solveflag = 3
!         gr_mgBndConditions = MG_BND_PERIODIC

     else

         solveflag = -1

     endif

  endif


  !We will refer to these module level variables in gr_mgPfftInitGrid:
  gr_mgPfftSolveFlag = solveFlag
  gr_mgPfftLastMappedLevel = NONEXISTENT
  
end subroutine gr_mgPfftInit
