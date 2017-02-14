!!****if* source/Grid/GridMain/paramesh/paramesh4/Paramesh4dev/flash_avoid_orrery/gr_initParameshArrays
!!
!! NAME
!!
!!  gr_initParameshArrays
!!
!! SYNOPSIS
!!
!!  call gr_initParameshArrays(logical(IN) :: restart,
!!                             integer(IN) :: xlboundary,
!!                             integer(IN) :: xrboundary,
!!                             integer(IN) :: ylboundary,
!!                             integer(IN) :: yrboundary,
!!                             integer(IN) :: zlboundary,
!!                             integer(IN) :: zrboundary
!!                             )
!!
!! DESCRIPTION
!!
!!  Perform early initialization of some Grid data structures.
!!
!!  This routine prepares the Grid for being filled with
!!  meaningful data.
!!
!! ARGUMENTS
!!
!!   restart -   Is the grid being prepared for initialization with
!!               data from a checkpoint file?
!!   xlboundary - boundary condition type of outer domain boundary in lower X direction.
!!   xrboundary - boundary condition type of outer domain boundary in upper X direction.
!!   ylboundary - boundary condition type of outer domain boundary in lower Y direction.
!!   yrboundary - boundary condition type of outer domain boundary in upper Y direction.
!!   zlboundary - boundary condition type of outer domain boundary in lower Z direction.
!!   zrboundary - boundary condition type of outer domain boundary in upper Z direction.
!!
!! SEE ALSO
!!
!!  gr_initParameshDomainBboxes
!!
!!***
subroutine gr_initParameshArrays(restart,&
                                     &  xlboundary, xrboundary, &
                                     &  ylboundary, yrboundary, &
                                     &  zlboundary, zrboundary)

   use paramesh_dimensions
   use physicaldata
   use workspace
   use tree
   use paramesh_mpi_interfaces, ONLY : mpi_amr_global_domain_limits,&
                                       mpi_amr_boundary_block_info
   use paramesh_interfaces, ONLY : amr_refine_derefine,amr_guardcell
   use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, &
                         gr_imin, gr_imax, gr_jmin, gr_jmax, &
                         gr_kmin, gr_kmax, &
                         gr_nblockX, gr_nblockY, gr_nblockZ

   implicit none
#include "constants.h"
#include "Flash.h"
   logical,intent(IN) :: restart
   integer,intent(IN) :: xlboundary, xrboundary
   integer,intent(IN) :: ylboundary, yrboundary
   integer,intent(IN) :: zlboundary, zrboundary
   integer :: lnblocks_old 
   integer :: five,six

   call mpi_amr_global_domain_limits()

#if defined(BITTREE)
   call amr_build_bittree()
   if(.not. restart) then
      call amr_reorder_grid()
      call amr_refine_derefine()
   end if
   call gr_initParameshDomainBboxes(xlboundary, xrboundary, &
                                 &  ylboundary, yrboundary, &
                                 &  zlboundary, zrboundary)
   call amr_morton_process()
#else
! Make sure that blocks produced by divide_domain are in strict
! morton order
   if(restart) then
      call gr_initParameshDomainBboxes(xlboundary, xrboundary, &
                                    &  ylboundary, yrboundary, &
                                    &  zlboundary, zrboundary)
      call amr_morton_process()
   else
      call gr_initParameshDomainBboxes(xlboundary, xrboundary, &
                                    &  ylboundary, yrboundary, &
                                    &  zlboundary, zrboundary)
      call amr_refine_derefine()
   end if
#endif


    !CD: Inform PARAMESH that we no longer need orrery.
    surr_blks_valid = .true.


    if(restart) then
       grid_analysed_mpi=1
       call mpi_amr_boundary_block_info(gr_meshMe, gr_meshNumProcs)
    end if
  ! reset for quadratic interpolation
  
  interp_mask_unk(:) = 1
  interp_mask_work(:) = 1
  
  return
end subroutine gr_initParameshArrays

