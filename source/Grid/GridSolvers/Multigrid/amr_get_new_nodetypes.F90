!!****if* source/Grid/GridSolvers/Multigrid/amr_get_new_nodetypes
!!
!! NAME
!!  amr_get_new_nodetypes
!!
!! SYNOPSIS
!!
!!  amr_get_new_nodetypes(integer(IN)::nprocs,
!!                        integer(IN)::mype,
!!                        integer(IN)::level
!!
!! DESCRIPTION
!!
!!  Previous to calling this routine, the contents of the 
!!  nodetypes array should be set to the desired maximum
!!  level.  This routine calls paramesh-internal routines
!!  for rebuilding communication patterns such that
!!  the grid may operate with the new nodetypes set.
!!
!!
!! ARGUMENTS
!!
!! nprocs - number of processors
!! mype   - local processor number
!! level  - the present maximum refinement level
!!
!!
!! NOTES
!!
!! adapted from PARAMESH - an adaptive mesh library.
!! Copyright (C) 2003
!!
!! Use of the PARAMESH software is governed by the terms of the
!! usage agreement which can be found in the file
!! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!!***!----------------------------------------------------------------------

subroutine amr_get_new_nodetypes (nprocs,mype,level)

  implicit none

  integer, intent(in) :: nprocs, mype,level
  integer :: tag_offset
  logical :: lec, lnc, lfulltree

  call mpi_amr_read_guard_comm_mg(nprocs, level)
  call mpi_amr_read_prol_comm_mg(nprocs, level)
  call mpi_amr_read_flux_comm_mg(nprocs, level)
  call mpi_amr_read_restrict_comm_mg(nprocs, level)

  return
end subroutine amr_get_new_nodetypes
