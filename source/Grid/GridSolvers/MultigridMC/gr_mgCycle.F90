!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgCycle
!!
!! NAME
!!
!!  gr_mgCycle
!!
!! SYNOPSIS
!!
!!  call gr_mgCycle(integer(in) :: level,
!!                  integer(in) :: img_soln,
!!                  integer(in) :: img_src,
!!                  integer(in) :: img_res,
!!                  integer(in) :: img_corr,
!!                  integer(in) :: img_temp,
!!                  integer(in) :: img_temp2,
!!                   :: mg_solve,
!!                   :: mg_residual,
!!                   :: mg_relax)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   level : 
!!
!!   img_soln : 
!!
!!   img_src : 
!!
!!   img_res : 
!!
!!   img_corr : 
!!
!!   img_temp : 
!!
!!   img_temp2 : 
!!
!!   mg_solve : 
!!
!!   mg_residual : 
!!
!!   mg_relax : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_cycle()

!  Description: Perform one multigrid cycle beginning at a specified mesh
!               level.  This version implements a V-cycle for adaptively
!               refined meshes (Martin, D. and Cartwright, K.  "Solving
!               Poisson's Equation using Adaptive Mesh Refinement," 1996).


      subroutine gr_mgCycle (level, img_soln, img_src, & 
     &                       img_res, img_corr, img_temp, img_temp2, & 
     &                       mg_solve, mg_residual, mg_relax)

!===============================================================================

      use Grid_data, ONLY : gr_meshMe

      use RuntimeParameters_interface, ONLY : RuntimeParameters_get

!  use gr_mgData, only : writeflag_mg

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_getListOfBlocks,   &
                             Grid_getLocalNumBlks

  use tree, only : lrefine,nodetype

      implicit none

#include "Flash.h"
#include "constants.h"


  integer, intent(in) :: level, img_soln, img_src, & 
     &           img_res, img_corr, img_temp, img_temp2

  external mg_solve, mg_residual, mg_relax

  integer :: i 

  real, pointer, dimension(:,:,:,:,:) :: unk

  integer, save :: mgrid_npresmooth,mgrid_npossmooth
  logical, save :: first_call = .true.

  integer SolveLevel

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  integer lb,blockcount,ii,jj,kk,blockID

  real bsize(MDIM),coord(MDIM)

  integer blockList(MAXBLOCKS)
  real, pointer, dimension(:,:,:,:) :: solnData

  integer lnblocks2

  !==========================================================================

  if (first_call) then
     call RuntimeParameters_get('mgrid_npresmooth',    mgrid_npresmooth)
     call RuntimeParameters_get('mgrid_npossmooth',    mgrid_npossmooth)
     first_call = .false.
  end if

  call gr_mgZero (level, img_corr, 0)

  !-------------------------------------------------------------------------------

      SolveLevel = 1

!     Fine -> coarse leg of V.
      if (level > SolveLevel) then

        do i = level, SolveLevel+1, -1

           call gr_mgCopy (i, img_soln, img_temp, 1)

           call gr_mgZero (i-1, img_corr, 0)

           call mg_relax (i, img_res, img_corr, mgrid_npresmooth)

           call gr_mgCorrect (i, img_soln, img_corr, 1)

           call mg_residual (i, img_res, img_corr, img_temp2, 0, 0)

           call gr_mgRestrict (i, img_temp2, img_res)

           call mg_residual (i-1, img_src, img_soln, img_res, 1, 1)

        enddo

      endif

!-------------------------------------------------------------------------------

!               Solve on coarsest level.
      call mg_solve (SolveLevel, img_res, img_corr, mg_relax)
      call gr_mgCorrect (SolveLevel, img_soln, img_corr, 1)

!-------------------------------------------------------------------------------

!               Coarse -> fine leg of V.

      if (level > SolveLevel) then

        do i = SolveLevel+1, level

           call gr_mgProlong (i-1, img_corr, img_corr, 1)

           call mg_residual (i, img_res, img_corr, img_res, 0, 0)

           call gr_mgZero (i, img_temp2, 0)

           call gr_mgZero (i-1, img_temp2, 0)

           call mg_relax (i, img_res, img_temp2, mgrid_npossmooth)

           call gr_mgCorrect (i, img_corr, img_temp2, 0)

           call gr_mgCopy (i, img_temp, img_soln, 1)

           call gr_mgCorrect (i, img_soln, img_corr, 1)

        enddo

      endif
!-------------------------------------------------------------------------------


!===============================================================================

      return
      end






