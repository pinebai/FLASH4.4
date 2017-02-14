!!****ih* source/Grid/localAPI/gr_hgInterface
!!
!! NAME
!!  gr_hgInterface
!!
!! SYNOPSIS
!!  use gr_hgInterface
!!
!! DESCRIPTION
!! This is the interface module for the multigrid solver implementation.
!!
!! Interfaces for subroutines are defined here.
!! Most of these are meant to be called only within the current (Paramesh) multigrid
!! solver implementation.
!! 
!!***

Module gr_hgInterface
  implicit none

  interface
     subroutine gr_hgBndry (level, ivar, nlayers, leafOnly, iopt, call, extrap)
       implicit none
       integer, intent(IN) :: level, ivar, nlayers, leafOnly, iopt, call
       logical, intent(IN) :: extrap
     end subroutine gr_hgBndry
  end interface

  interface
     Subroutine gr_hgGuardCell(myPE, nlayers, extrap)
       implicit none
       integer, INTENT(in) :: myPE, nlayers
       logical, intent(in) :: extrap
     end Subroutine gr_hgGuardCell
  end interface

  interface
     subroutine gr_hgInit
       implicit none
     end subroutine gr_hgInit
  end interface

  interface
     subroutine gr_hgInitSource (isrc, isoln)
       implicit none
       integer, intent(in) :: isrc, isoln
     end subroutine gr_hgInitSource
  end interface

  interface
     subroutine gr_hgLevelAdd(level, ivar1, ivar2, LeafFlag)
       implicit none
       integer, intent(in) :: ivar1, ivar2, level, LeafFlag
     end subroutine gr_hgLevelAdd
  end interface

  interface
     subroutine gr_hgLevelAddScalar(level, ivar, scalar, LeafFlag)
       implicit none
       integer, intent(in) :: ivar, level, LeafFlag
       real, intent(in)    :: scalar
     end subroutine gr_hgLevelAddScalar
  end interface

  interface
     subroutine gr_hgLevelMultiplyScalar(level, ivar, scalar, LeafFlag)
       implicit none
       integer, intent(in) :: ivar, level, LeafFlag
       real, intent(in)    :: scalar
     end subroutine gr_hgLevelMultiplyScalar
  end interface

  interface
     subroutine gr_hgLevelZero(level, ivar, LeafFlag)
       implicit none
       integer, intent(in) :: ivar, level, LeafFlag
     end subroutine gr_hgLevelZero
  end interface

  interface
     subroutine gr_hgMapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct, &
          axis,face,idest)
       integer, intent(out) :: bcTypeToApply
       integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct,axis,face
       integer,intent(in),OPTIONAL:: idest
     end subroutine gr_hgMapBcType
  end interface

  interface
     subroutine gr_hgNorm (level, normType, ivar, norm, leafOnly)
       implicit none
       integer, intent(IN)  :: normType, level, ivar, leafOnly
       real, intent(OUT)    :: norm
     end subroutine gr_hgNorm
  end interface

  interface
     subroutine gr_hgProlongBndries(level, ifrom, ito, ichild)
       integer, intent(in)          :: ifrom, ito, level, ichild
     end subroutine gr_hgProlongBndries
  end interface

  interface
     subroutine gr_hgResidual(level, gr_iSource, gr_iSoln, ires,dt,chi,theta)
       implicit none
       integer, intent(in)          :: level, gr_iSource, gr_iSoln, ires
       real,intent(IN),OPTIONAL     :: dt, chi, theta
     end subroutine gr_hgResidual
  end interface

  interface
     subroutine gr_hgRestoreNodetypes ()
       implicit none
     end subroutine gr_hgRestoreNodetypes
  end interface

  interface
     subroutine gr_hgRestrict(level, ito, ifrom)
       implicit none
       integer, intent(in)          :: ifrom, ito, level
     end subroutine gr_hgRestrict
  end interface

  interface
     subroutine gr_hgSetZeroBoundary(level, gr_iSoln)
       implicit none
       integer, intent(in) :: gr_iSoln, level
     end subroutine gr_hgSetZeroBoundary
  end interface

  interface
     subroutine gr_hgSetExtBoundary (idiag, idir, extrap)
       implicit none
       integer, intent(IN) :: idiag, idir
       logical, intent(in) :: extrap
     end subroutine gr_hgSetExtBoundary
  end interface

  interface
     subroutine gr_hgSetMaxLevel (level)
       implicit none
       integer,intent(IN) :: level
     end subroutine gr_hgSetMaxLevel
  end interface

  interface
     subroutine gr_hgSolve(gr_iSource, gr_iSoln, gr_iSls, gr_iCorr, SolveBlock, bndTypes, src_fact, dt, chi)
       implicit none
       integer, intent(in) :: gr_iSource, gr_iSoln, gr_iSls, gr_iCorr
       integer, intent(in) :: bndTypes(6)
       real, intent(IN),OPTIONAL :: src_fact, dt, chi
       external               SolveBlock
     end subroutine gr_hgSolve
  end interface

  interface
     subroutine gr_hgSolveLevel(level, gr_iSource, gr_iSoln, SolveBlock, LeafFlag, dt, chi, theta)
       implicit none
       integer, intent(in) :: gr_iSource, gr_iSoln, level, LeafFlag
       real, intent(IN),OPTIONAL :: dt, chi, theta
       external               SolveBlock
     end subroutine gr_hgSolveLevel
  end interface

  interface
     subroutine gr_hgFinalize()
       implicit none
     end subroutine gr_hgFinalize
  end interface


!The Huang-Greengard Multigrid solver will soon be able 
!to use PFFT.  These subroutines are given the 
!hgPfft sub-unit prefix.
!------------------------------------------------------
  interface
     subroutine gr_hgPfftInit()
       implicit none
     end subroutine gr_hgPfftInit
  end interface

  interface
     subroutine gr_hgPfftInitGrid(refinementLevel, gridChanged, poisfact)
       implicit none
       integer, intent(IN) :: refinementLevel, gridChanged
       real, intent(IN) :: poisfact
     end subroutine gr_hgPfftInitGrid
  end interface

  interface
     subroutine gr_hgPfftSolveLevel(iSrc, iSoln, level)
       implicit none
       integer, intent(in)    :: iSrc, iSoln, level
     end subroutine gr_hgPfftSolveLevel
  end interface

  interface
     subroutine gr_hgPfftFinalize()
       implicit none
     end subroutine gr_hgPfftFinalize
  end interface

  interface
     subroutine gr_hgPfftFinalizeGrid()
       implicit none
     end subroutine gr_hgPfftFinalizeGrid
  end interface

!------------------------------------------------------

end Module gr_hgInterface
