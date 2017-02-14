!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_flux_conserve
!! NAME
!!
!!   amr_flux_conserve
!!
!! SYNOPSIS
!!
!!   call amr_flux_conserve (mype, nsub)
!!   call amr_flux_conserve (mype, nsub, flux_dir)
!!
!!   call amr_flux_conserve (integer, integer, optional integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype  
!!     The calling processor.
!!
!!   integer, intent(in) :: nsub          
!!     The current time subcycle. If this is 1 then this info is used to 
!!     reset the temporary boundary flux arrays to 0. This argument only has
!!     an effect if variable time steps are being used.
!!
!!   optional, integer, intent(in) :: flux_dir
!!     Option integer which selects which coordinate direction to apply
!!     the flux conservation operation to:
!!     If flux_dir = 1 -> x direction
!!        flux_dir = 2 -> y direction
!!        flux_dir = 3 -> z direction
!!     If this argument is not specified, then the default behaviour is
!!     to operate on all directions.  Using this argument can be useful
!!     for Strang-split schemes to improve performance.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!! 
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!
!! CALLS
!! 
!!   amr_flux_conserve_udt
!!   amr_flux_conserve_vdt
!!    
!! RETURNS
!!
!!   Does not return anything.  Upon exit the fluxes stored in the arrays
!!   flux_x, flux_y, or flux_z are corrected at jumps in refinement.  This
!!   is either an averaging proceedure or a sum as selected by the user
!!   by adjusting the preprocessor variables which control this behaviour
!!   in 'paramesh_preprocessor.fh'.  This can also be controled through the 
!!   equivalent logical variables read at runtime from the file 
!!   'amr_runtime_parameters' if the 'LIBRARY' preprocessor flag is defined.
!!
!! DESCRIPTION
!!
!!   This is a wrapper routine which makes the appropriate call to the
!!   routines which manage flux conservation at the boundaries between
!!   grid blocks of different refinement level.
!! 
!!   These routines get block boundary data from neighbors who are
!!   parents of leaf blocks. This is required in flux conserving schemes
!!   where the coarser block needs to use the same fluxes and mean pressures
!!   as will be used on the finer blocks across their shared boundary.
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson for 
!!   directional guardcell filling.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      subroutine amr_flux_conserve(mype,nsub,flux_dir)

      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_flux_conserve_udt, & 
     &                                amr_flux_conserve_vdt

      implicit none

      integer, intent(in)  ::  mype,nsub
      integer, optional, intent(in) :: flux_dir

      integer :: lb


!------------------------------------


      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

      if(nodetype(lb).eq.1 .or. advance_all_levels) then

! Store fluxes in temporary storage
       tflux_x(:,:,:,:,lb) = flux_x(:,:,:,:,lb)
       if (ndim >= 2) then
          tflux_y(:,:,:,:,lb) = flux_y(:,:,:,:,lb)
       end if
       if (ndim == 3) then
          tflux_z(:,:,:,:,lb) = flux_z(:,:,:,:,lb)
       end if

      endif

      enddo
      endif

      if (var_dt) then
      call amr_flux_conserve_vdt(mype,nsub) ! called if variable dt
      else
      call amr_flux_conserve_udt(mype,flux_dir)      ! called if uniform dt
      endif

      return
      end subroutine amr_flux_conserve
