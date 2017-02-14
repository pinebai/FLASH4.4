!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_edge_average
!! NAME
!!
!!   amr_edge_average
!!
!! SYNOPSIS
!!
!!   call amr_edge_average(mype, lfullblocks, nsub)
!!
!!   call amr_edge_average(integer, logical, integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype
!!     The calling processors id.
!!
!!   integer, intent(in) :: nsub
!!     Controls something for variable timestep algorithms ???.
!!
!!   logical, intent(in) :: lfullblock
!!     Controls which arrays this routine works on.  If .true. then the unk_e_x,y,z
!!     arrays are used.  Otherwise, the bedge_* arrays are used.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_edge_average_udt
!!   amr_edge_diagonal_check
!!   amr_edge_average_vdt
!!
!! RETURNS
!!
!!   Nothing returned.
!!
!! DESCRIPTION
!!
!!  This is a wrapper routine which makes the appropriate calls to the
!!  routines which manage edge data consistency at the boundaries between
!!  grid blocks of different refinement level.  This routine is called by the user
!!  when they wish to enforce consistency of edge centered variables (either stored
!!  in the bedge_* arrays or in the unk_e_* arrays) at jumps in refinement.  This 
!!  routine is analogous to the subroutine 'amr_flux_conserve'.
!!
!!  If lfullblock is true, then this routine will update the main edge-centered
!!  datastructure arrays, unk_e_x[y][z], otherwise it simply operates on
!!  the temporary edge-centered data computed on block boundary faces only.
!! 
!! AUTHORS
!!
!!  Written by Peter MacNeice (July 1997) and modified February 2001.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine amr_edge_average(mype,lfullblock,nsub)


      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_edge_average_udt, & 
     &                                amr_edge_diagonal_check & 
     &                               ,amr_edge_average_vdt

      implicit none

      integer, intent(in)  ::  mype,nsub
      logical, intent(in)  ::  lfullblock

      integer :: lb

      integer :: nguard0

!------------------------------------

      nguard0 = nguard*npgs

      if (ndim > 1) then

      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

      if(nodetype(lb).eq.1 .or. advance_all_levels) then


! If this is to be applied to the full edge-centered data then copy the block
! faces in to block boundary edge-centered datastructure.
        if(lfullblock) then

          bedge_facex_y(1:nedges,                                      & 
                        1,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb) =        & 
                unk_e_y(1:nedges,                                      & 
                        1+nguard0,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb)
          bedge_facex_y(1:nedges,                                      & 
                        2,jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb) =        & 
                unk_e_y(1:nedges,nxb+1+nguard0,                        & 
                        jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb)
          If (ndim == 3.Or.(ndim == 2 .and. l2p5d == 1)) Then
          bedge_facex_z(1:nedges,1,jl_bnd:ju_bnd+k2d,                  & 
                        kl_bnd:ku_bnd,lb) =                            & 
                unk_e_z(1:nedges,1+nguard0,jl_bnd:ju_bnd+k2d,          & 
                        kl_bnd:ku_bnd,lb)
          bedge_facex_z(1:nedges,2,jl_bnd:ju_bnd+k2d,                  & 
                        kl_bnd:ku_bnd,lb) =                            & 
                unk_e_z(1:nedges,nxb+1+nguard0,                        & 
                        jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,lb)
          End If
          bedge_facey_x(1:nedges,il_bnd:iu_bnd,1,                      & 
                        kl_bnd:ku_bnd+k3d,lb) =                        & 
                unk_e_x(1:nedges,il_bnd:iu_bnd,1+nguard0*k2d,          & 
                        kl_bnd:ku_bnd+k3d,lb)
          bedge_facey_x(1:nedges,il_bnd:iu_bnd,1+k2d,                  & 
                        kl_bnd:ku_bnd+k3d,lb) =                        & 
                unk_e_x(1:nedges,il_bnd:iu_bnd,                        & 
                        nyb+k2d+nguard0*k2d,                           & 
                        kl_bnd:ku_bnd+k3d,lb)
          If (ndim == 3.Or. (ndim == 2 .and. l2p5d == 1)) Then
          bedge_facey_z(1:nedges,il_bnd:iu_bnd+1,1,                    & 
                        kl_bnd:ku_bnd,lb) =                            & 
                unk_e_z(1:nedges,il_bnd:iu_bnd+1,1+nguard0*k2d,        & 
                        kl_bnd:ku_bnd,lb)
          bedge_facey_z(1:nedges,il_bnd:iu_bnd+1,1+k2d,                & 
                        kl_bnd:ku_bnd,lb) =                            & 
                unk_e_z(1:nedges,il_bnd:iu_bnd+1,                      & 
                        nyb+k2d+nguard0*k2d,                           & 
                        kl_bnd:ku_bnd,lb)

          bedge_facez_x(1:nedges,il_bnd:iu_bnd, & 
     &                  jl_bnd:ju_bnd+k2d,1,lb) =  & 
     &          unk_e_x(1:nedges,il_bnd:iu_bnd, & 
     &                  jl_bnd:ju_bnd+k2d,1+nguard0*k3d,lb)
          bedge_facez_y(1:nedges,il_bnd:iu_bnd+1, & 
     &                  jl_bnd:ju_bnd,1,lb) =  & 
     &          unk_e_y(1:nedges,il_bnd:iu_bnd+1, & 
     &                  jl_bnd:ju_bnd,1+nguard0*k3d,lb)

          if(ndim.eq.3) then
          bedge_facez_x(1:nedges,il_bnd:iu_bnd, & 
     &                  jl_bnd:ju_bnd+k2d,1+k3d,lb) =  & 
     &          unk_e_x(1:nedges,il_bnd:iu_bnd, & 
     &                  jl_bnd:ju_bnd+k2d,nzb+(1+nguard0)*k3d,lb)
          bedge_facez_y(1:nedges,il_bnd:iu_bnd+1, & 
     &                  jl_bnd:ju_bnd,1+k3d,lb) =  & 
     &          unk_e_y(1:nedges,il_bnd:iu_bnd+1, & 
     &                  jl_bnd:ju_bnd,nzb+(1+nguard0)*k3d,lb)
          endif
          endif

        endif

        tbedge_facex_y(:,:,:,:,lb) = bedge_facex_y(:,:,:,:,lb)
        tbedge_facex_z(:,:,:,:,lb) = bedge_facex_z(:,:,:,:,lb)
        tbedge_facey_x(:,:,:,:,lb) = bedge_facey_x(:,:,:,:,lb)
        tbedge_facey_z(:,:,:,:,lb) = bedge_facey_z(:,:,:,:,lb)
        If (ndim == 3.Or.(ndim == 2 .and. l2p5d == 1)) Then
        tbedge_facez_x(:,:,:,:,lb) = bedge_facez_x(:,:,:,:,lb)
        tbedge_facez_y(:,:,:,:,lb) = bedge_facez_y(:,:,:,:,lb)
        endif

        endif

      enddo
      endif

!---------------------------
! Operate on block boundary edge-centered data


      if (var_dt) then
      call amr_edge_average_vdt(mype,nsub) ! called if variable dt
      else
      call amr_edge_average_udt(mype)      ! called if uniform dt
      endif

!---------------------------
! amr_edge_diagonal_check works for either timestepping strategy.
      call amr_edge_diagonal_check(mype)


!---------------------------


! If this is to be applied to the full edge-centered data then copy the 
! modified data on block face boundaries back to the edge-centered 
! datastructure.
      if(lfullblock) then

        if(lnblocks.gt.0) then
        do lb = 1,lnblocks

        if(nodetype(lb).eq.1 .or. advance_all_levels) then

          unk_e_y(1:nedges,1+nguard0,jl_bnd:ju_bnd, & 
     &            kl_bnd:ku_bnd+k3d,lb)     =  & 
     &      bedge_facex_y(1:nedges,1,jl_bnd:ju_bnd, & 
     &                    kl_bnd:ku_bnd+k3d,lb)
          unk_e_y(1:nedges,nxb+1+nguard0,jl_bnd:ju_bnd, & 
     &            kl_bnd:ku_bnd+k3d,lb) =  & 
     &      bedge_facex_y(1:nedges,2,jl_bnd:ju_bnd, & 
     &                    kl_bnd:ku_bnd+k3d,lb)
          if(ndim.eq.3.or.(ndim == 2 .AND. l2p5d.eq.1)) then
          unk_e_z(1:nedges,1+nguard0,jl_bnd:ju_bnd+k2d, & 
     &                    kl_bnd:ku_bnd,lb)     =  & 
     &      bedge_facex_z(1:nedges,1,jl_bnd:ju_bnd+k2d, & 
     &                    kl_bnd:ku_bnd,lb)
          unk_e_z(1:nedges,nxb+1+nguard0,jl_bnd:ju_bnd+k2d, & 
     &                    kl_bnd:ku_bnd,lb) =  & 
     &      bedge_facex_z(1:nedges,2,jl_bnd:ju_bnd+k2d, & 
     &                    kl_bnd:ku_bnd,lb)
          endif
          unk_e_x(1:nedges,il_bnd:iu_bnd,1+nguard0*k2d, & 
     &                     kl_bnd:ku_bnd+k3d,lb)     =  & 
     &      bedge_facey_x(1:nedges,il_bnd:iu_bnd,1, & 
     &                    kl_bnd:ku_bnd+k3d,lb)
          unk_e_x(1:nedges,il_bnd:iu_bnd,ju_bnd+k2d, & 
     &                     kl_bnd:ku_bnd+k3d,lb) =  & 
     &      bedge_facey_x(1:nedges,il_bnd:iu_bnd,2, & 
     &                    kl_bnd:ku_bnd+k3d,lb)
          if(ndim.eq.3.or.(ndim == 2 .AND. l2p5d.eq.1)) then
          unk_e_z(1:nedges,il_bnd:iu_bnd+1,1+nguard0*k2d, & 
     &                      kl_bnd:ku_bnd,lb)     =  & 
     &      bedge_facey_z(1:nedges,il_bnd:iu_bnd+1,1, & 
     &                    kl_bnd:ku_bnd,lb)
          unk_e_z(1:nedges,il_bnd:iu_bnd+1, & 
     &                    nyb+k2d+nguard0*k2d, & 
     &                    kl_bnd:ku_bnd,lb) =  & 
     &      bedge_facey_z(1:nedges,il_bnd:iu_bnd+1,1+k2d, & 
     &                    kl_bnd:ku_bnd,lb)

          unk_e_x(1:nedges,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &            1+nguard0*k3d,lb)     =  & 
     &      bedge_facez_x(1:nedges,il_bnd:iu_bnd, & 
     &                    jl_bnd:ju_bnd+k2d,1,lb)
          unk_e_y(1:nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd, & 
     &            1+nguard0*k3d,lb)     =  & 
     &      bedge_facez_y(1:nedges,il_bnd:iu_bnd+1, & 
     &                    jl_bnd:ju_bnd,1,lb)

          if(ndim.eq.3) then
          unk_e_x(1:nedges,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &            nzb+(1+nguard0)*k3d,lb) =  & 
     &      bedge_facez_x(1:nedges,il_bnd:iu_bnd, & 
     &                    jl_bnd:ju_bnd+k2d,1+k3d,lb)
          unk_e_y(1:nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd, & 
     &            nzb+(1+nguard0)*k3d,lb) =  & 
     &      bedge_facez_y(1:nedges,il_bnd:iu_bnd+1, & 
     &                    jl_bnd:ju_bnd,1+k3d,lb)
          endif
          endif

        endif

        enddo
        endif

      endif

      end if ! if (ndim > 1) 

      return
      end subroutine amr_edge_average
