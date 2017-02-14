!!****ih* source/Grid/localAPI/gr_mgInterface
!!
!! NAME
!!  gr_mgInterface
!!
!! SYNOPSIS
!!  use gr_mgInterface
!!
!! DESCRIPTION
!! This is the interface module for the multigrid solver implementation.
!!
!! Interfaces for subroutines are defined here.
!! Most of these are meant to be called only within the current (Paramesh) multigrid
!! solver implementation.
!! 
!!***

Module gr_mgInterface
  implicit none

  interface
     subroutine  gr_mgBndry(level, ivar, nlayers, leaf_only, iopt, call)
       implicit none
       integer, intent(in) :: level, ivar, nlayers, leaf_only, iopt, call
     end subroutine gr_mgBndry
  end interface

  interface
     subroutine gr_mgCopy (level, ifrom, ito, leaf_only)
      implicit none
      integer, intent(in) :: level, ifrom, ito, leaf_only
     end subroutine gr_mgCopy
  end interface

  interface
     subroutine gr_mgCorrect (level, isoln, icorr, leaf_only)
      implicit none
      integer, intent(in) :: level, isoln, icorr, leaf_only
    end subroutine gr_mgCorrect
  end interface

  interface
     subroutine gr_mgCycle (level, img_soln, img_src, & 
    &                       img_res, img_corr, img_temp, img_temp2, & 
    &                       mg_solve, mg_residual, mg_relax)
      implicit none
      integer, intent(in) :: level, img_soln, img_src, & 
    &            img_res, img_corr, img_temp, img_temp2
      external mg_solve, mg_residual, mg_relax 
     end subroutine gr_mgCycle
  end interface

  interface
     subroutine gr_mgGuardcell(mype2, ivar, nlayers, simtime, idiag, idir)
      implicit none
      integer, intent(in) :: mype2, nlayers, idiag, idir
      integer, intent(in) :: ivar
      real,    intent(in) :: simtime
     end subroutine gr_mgGuardcell
  end interface

  interface
     subroutine gr_mgInit()
      implicit none
     end subroutine gr_mgInit
  end interface

  interface
     subroutine gr_mgInitSlv(bndTypes)
      implicit none
      integer, intent(in) :: bndTypes(6)
     end subroutine gr_mgInitSlv
  end interface

  interface
     subroutine gr_mgFinalize()
      implicit none
     end subroutine gr_mgFinalize
  end interface

  interface
     subroutine gr_mgPfftInit()
      implicit none
     end subroutine gr_mgPfftInit
  end interface

  interface
     subroutine gr_mgPfftFinalize()
      implicit none
     end subroutine gr_mgPfftFinalize
  end interface

  interface
     subroutine gr_mgInitSrc (isrc_dens, poisfact, img_src, img_soln)
      implicit none
      integer, intent(in) :: isrc_dens, img_src, img_soln
      real, intent(in)    :: poisfact
     end subroutine gr_mgInitSrc
  end interface

 interface
     subroutine gr_mgMapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct, &
          axis,face,idest)
       integer, intent(out) :: bcTypeToApply
       integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct,axis,face
       integer,intent(in),OPTIONAL:: idest
     end subroutine gr_mgMapBcType
  end interface

  interface
     subroutine gr_mgNorm (level, ivar, norm, leaf_only)
      implicit none
      integer, intent(in) :: level, ivar, leaf_only
      real, intent(inout) :: norm
     end subroutine gr_mgNorm
  end interface

  interface
     subroutine gr_mgProlong (level, ifrom, ito, add)
      implicit none
      integer, intent(in) :: level, ifrom, ito, add
     end subroutine gr_mgProlong
  end interface

  interface
     subroutine gr_mgRestrict (level, ifrom, ito)
      implicit none
      integer, intent(in) :: level, ifrom, ito
     end subroutine gr_mgRestrict
  end interface

  interface
     subroutine gr_mgZero (level, ivar, leaf_only)
      implicit none
      integer, intent(in) :: level, ivar, leaf_only
     end subroutine
  end interface

  interface
     subroutine mg_restore_nodetypes (level)
      implicit none
      integer, intent(in) :: level
     end subroutine mg_restore_nodetypes
  end interface

  interface
     subroutine gr_mgSolve (isrc_dens, img_soln, poisfact, img_src, img_res, & 
                           img_corr, img_temp, img_temp2, bc_types,mg_solve,  &
                           mg_residual, mg_residual_leafs, mg_relax)
      implicit none
      integer,intent(in)   :: isrc_dens, img_soln, img_src, img_res, img_corr, img_temp, &
                              img_temp2
      integer,intent(in)   :: bc_types(:)
      real,intent(in)      :: poisfact
      external mg_solve, mg_residual, mg_relax
      external mg_residual_leafs
     end subroutine gr_mgSolve
  end interface

  interface
     subroutine gr_mgPfftFinalizeGrid()
      implicit none
     end subroutine gr_mgPfftFinalizeGrid
  end interface


! MC - Pfft subroutines



!------------------------------------------------------

end Module gr_mgInterface
