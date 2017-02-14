!!****if* source/Grid/GridSolvers/gr_solversData
!!
!! NAME
!!
!!  gr_solversData
!!
!! SYNOPSIS
!!  use gr_solversData, ONLY: gr_solversDbgContext
!!
!! DESCRIPTION
!!
!!  Defines storage for some data items that are common to
!!  several GridSolvers implementation.
!!
!!  
!!***

Module gr_solversData 
  
  use gr_interfaceTypeDecl, ONLY: gr_solversDbgContext_t
  
  implicit none

  ! Structure that holds context information on the current operation,
  ! for debugging
  type(gr_solversDbgContext_t), save :: gr_solversDbgContext
  
end Module gr_solversData
