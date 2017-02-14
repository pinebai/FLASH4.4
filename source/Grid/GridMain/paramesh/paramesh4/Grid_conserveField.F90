!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_conserveField
!!
!! NAME
!!
!!  Grid_conserveField
!!
!! SYNOPSIS
!!
!!  Grid_conserveField ()
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!! 
!!  Corrects electric fields at refinement jump boundaries to make 
!!  sure electric fields at the common boundaries are consistent 
!!  at both refined and coarse meshes.
!!
!!***

!!REORDER(5): bedge_facey_x, bedge_facez_x, bedge_facex_y,  bedge_facez_y
!!REORDER(5): bedge_facex_z,  bedge_facey_z, tbedge_facey_x, tbedge_facez_x
!!REORDER(5): bedge_facex_y, tbedge_facez_y, tbedge_facex_z, tbedge_facey_z
!!REORDER(4): scratchData

subroutine Grid_conserveField ()

  use paramesh_dimensions, ONLY : il_bnd, jl_bnd, kl_bnd, iu_bnd, ju_bnd, ku_bnd,&
                                  nguard, nxb, nyb, nzb, npgs

  use Grid_interface,      ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr

  use physicaldata,        ONLY :  bedge_facey_x,  bedge_facez_x, &
                                   bedge_facex_y,  bedge_facez_y, &
                                   bedge_facex_z,  bedge_facey_z, &
                                  tbedge_facey_x, tbedge_facez_x, &
                                  tbedge_facex_y, tbedge_facez_y, &
                                  tbedge_facex_z, tbedge_facey_z, &
                                  no_permanent_guardcells,        &
                                  advance_all_levels

  use tree,                ONLY : nodetype, lnblocks

  use paramesh_interfaces, only : amr_edge_average_udt, &
                                  amr_edge_diagonal_check
  use Grid_data, ONLY : gr_meshMe

  implicit none

#include "Flash.h"
#include "constants.h"


#if FLASH_NEDGE_VAR > 0
#include "UHD.h"

  integer :: lb,nguard0
  integer :: nsub,gridDataStruct
  logical :: lfullblock
  real, pointer, dimension(:,:,:,:) :: scratchData

  nguard0=nguard*npgs

  if (lnblocks .gt. 0) then
  ! Loop over all the paramesh blocks
  do lb=1,lnblocks

     ! Loop over the leaf blocks only
     if (nodetype(lb) .eq. 1 .or. advance_all_levels) then
        
        ! Get saved electric field data
        call Grid_getBlkPtr(lb,scratchData,SCRATCH)

        !! Ez on a face normal to x axis
        bedge_facex_z(1,1,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd,lb) = &
             scratchData(EZ_SCRATCH_GRID_VAR,1+nguard0,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)

        bedge_facex_z(1,2,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd,lb) = &
             scratchData(EZ_SCRATCH_GRID_VAR,nxb+nguard0+1,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd)

        !! Ez on a face normal to y axis
        bedge_facey_z(1,il_bnd:iu_bnd+1,1,kl_bnd:ku_bnd,lb) = &
             scratchData(EZ_SCRATCH_GRID_VAR,il_bnd:iu_bnd+1,1+nguard0,kl_bnd:ku_bnd)

        bedge_facey_z(1,il_bnd:iu_bnd+1,2,kl_bnd:ku_bnd,lb) = &
             scratchData(EZ_SCRATCH_GRID_VAR,il_bnd:iu_bnd+1,nyb+nguard0+1,kl_bnd:ku_bnd)

#if NDIM == 3
        !! Ex on a face normal to y axis
        bedge_facey_x(1,il_bnd:iu_bnd,1,kl_bnd:ku_bnd+1,lb) = &
             scratchData(EX_SCRATCH_GRID_VAR,il_bnd:iu_bnd,1+nguard0,kl_bnd:ku_bnd+1)

        bedge_facey_x(1,il_bnd:iu_bnd,2,kl_bnd:ku_bnd+1,lb) = &
             scratchData(EX_SCRATCH_GRID_VAR,il_bnd:iu_bnd,nyb+nguard0+1,kl_bnd:ku_bnd+1)

        !! Ex on a face normal to z axis
        bedge_facez_x(1,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,1,lb) = &
             scratchData(EX_SCRATCH_GRID_VAR,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,1+nguard0)

        bedge_facez_x(1,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,2,lb) = &
             scratchData(EX_SCRATCH_GRID_VAR,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,nzb+nguard0+1)

        !! Ey on a face normal to y axis
        bedge_facex_y(1,1,jl_bnd:ju_bnd,kl_bnd:ku_bnd+1,lb) = &
             scratchData(EY_SCRATCH_GRID_VAR,1+nguard0,jl_bnd:ju_bnd,kl_bnd:ku_bnd+1)

        bedge_facex_y(1,2,jl_bnd:ju_bnd,kl_bnd:ku_bnd+1,lb) = &
             scratchData(EY_SCRATCH_GRID_VAR,nxb+nguard0+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd+1)

        !! Ey on a face normal to z axis
        bedge_facez_y(1,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,1,lb) = &
             scratchData(EY_SCRATCH_GRID_VAR,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,1+nguard0)

        bedge_facez_y(1,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,2,lb) = &
             scratchData(EY_SCRATCH_GRID_VAR,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,nzb+nguard0+1)
#endif

        tbedge_facex_y(:,:,:,:,lb) = bedge_facex_y(:,:,:,:,lb)
        tbedge_facex_z(:,:,:,:,lb) = bedge_facex_z(:,:,:,:,lb)
        tbedge_facey_x(:,:,:,:,lb) = bedge_facey_x(:,:,:,:,lb)
        tbedge_facey_z(:,:,:,:,lb) = bedge_facey_z(:,:,:,:,lb)
        !if(ndim.eq.3.or.l2p5d.eq.1) then
        tbedge_facez_x(:,:,:,:,lb) = bedge_facez_x(:,:,:,:,lb)
        tbedge_facez_y(:,:,:,:,lb) = bedge_facez_y(:,:,:,:,lb)
        !endif

        call Grid_releaseBlkPtr(lb,scratchData,SCRATCH)
     endif
  enddo
  endif

  ! No need to consider full edge-centered data, that is, we only need to operate
  ! on the temporary edge-centered data computed on block boundary faces.
  lfullblock = .false.
  nsub       = 0
  
  ! Take averages on electric fields at the boundaries of the refined mesh
  ! and use the values to restrict electric fields at the boundaries of the
  ! coarse mesh.
  call amr_edge_average_udt(gr_meshMe)

! amr_edge_diagonal_check works for either timestepping strategy.
  call amr_edge_diagonal_check(gr_meshMe)


  if (lnblocks .gt. 0) then
  do lb=1,lnblocks
     if (nodetype(lb) .eq. 1 .or. advance_all_levels) then

        !Get a pointer for scratchData
        call Grid_getBlkPtr(lb,scratchData,SCRATCH)

        !! Ez on a face normal to x axis
        scratchData(EZ_SCRATCH_GRID_VAR,1+nguard0,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) = &
             bedge_facex_z(1,1,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd,lb)

        scratchData(EZ_SCRATCH_GRID_VAR,nxb+nguard0+1,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd) = &
             bedge_facex_z(1,2,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd,lb)

        !! Ez on a face normal to y axis
        scratchData(EZ_SCRATCH_GRID_VAR,il_bnd:iu_bnd+1,1+nguard0,kl_bnd:ku_bnd) = &
             bedge_facey_z(1,il_bnd:iu_bnd+1,1,kl_bnd:ku_bnd,lb)

        scratchData(EZ_SCRATCH_GRID_VAR,il_bnd:iu_bnd+1,nyb+nguard0+1,kl_bnd:ku_bnd) = &
             bedge_facey_z(1,il_bnd:iu_bnd+1,2,kl_bnd:ku_bnd,lb)

#if NDIM == 3
        !! Ex on a face normal to y axis
        scratchData(EX_SCRATCH_GRID_VAR,il_bnd:iu_bnd,1+nguard0,kl_bnd:ku_bnd+1) = &
             bedge_facey_x(1,il_bnd:iu_bnd,1,kl_bnd:ku_bnd+1,lb)

        scratchData(EX_SCRATCH_GRID_VAR,il_bnd:iu_bnd,nyb+nguard0+1,kl_bnd:ku_bnd+1) = &
             bedge_facey_x(1,il_bnd:iu_bnd,2,kl_bnd:ku_bnd+1,lb)

        !! Ex on a face normal to z axis
        scratchData(EX_SCRATCH_GRID_VAR,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,1+nguard0) = &
             bedge_facez_x(1,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,1,lb)

        scratchData(EX_SCRATCH_GRID_VAR,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,nzb+nguard0+1) = &
             bedge_facez_x(1,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,2,lb)

        !! Ey on a face normal to y axis
        scratchData(EY_SCRATCH_GRID_VAR,1+nguard0,jl_bnd:ju_bnd,kl_bnd:ku_bnd+1) = &
             bedge_facex_y(1,1,jl_bnd:ju_bnd,kl_bnd:ku_bnd+1,lb)

        scratchData(EY_SCRATCH_GRID_VAR,nxb+nguard0+1,jl_bnd:ju_bnd,kl_bnd:ku_bnd+1) = &
             bedge_facex_y(1,2,jl_bnd:ju_bnd,kl_bnd:ku_bnd+1,lb)

        !! Ey on a face normal to z axis
        scratchData(EY_SCRATCH_GRID_VAR,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,1+nguard0) = &
             bedge_facez_y(1,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,1,lb)

        scratchData(EY_SCRATCH_GRID_VAR,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,nzb+nguard0+1) = &
             bedge_facez_y(1,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,2,lb)

#endif
        call Grid_releaseBlkPtr(lb,scratchData,SCRATCH)
     endif
  enddo
  endif

#ifndef FLASH_GRID_PARAMESH2
  gridDataStruct = CENTER
#if NFACE_VARS > 0
  gridDataStruct = CENTER_FACES
#endif
  if (no_permanent_guardcells) then
     call gr_commSetup(gridDataStruct)
  end if
#endif

#endif !endif FLASH_NEDGE_VAR > 0

end subroutine Grid_conserveField


