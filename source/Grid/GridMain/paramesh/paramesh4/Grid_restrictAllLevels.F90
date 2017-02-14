!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_restrictAllLevels
!!
!! NAME
!!  Grid_restrictAllLevels
!!
!! SYNOPSIS
!! 
!!  Grid_restrictAllLevels()
!!  
!! DESCRIPTION 
!!  Restricts the grid data to all refinement levels. Normally FLASH
!!  only evolves on the leaf blocks, calling this routine makes all
!!  levels have valid data.  This is mostly for visualization purposes to
!!  be able to look at different levels of resolution.
!!  
!!  
!!
!!***


subroutine Grid_restrictAllLevels()

#include "Flash.h"
#include "constants.h"

  use Timers_interface, ONLY: Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_restrictByLevels
  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, gr_restrictAllMethod
  use paramesh_dimensions, ONLY: nfacevar, nvaredge, nvarcorn
  use physicaldata, ONLY: advance_all_levels
  use paramesh_mpi_interfaces, ONLY : mpi_amr_1blk_restrict, &
       mpi_morton_bnd_restrict

  implicit none
  interface
     subroutine mpi_amr_restrict_fulltree(mype,iopt,lcc,lfc,lec,lnc)
       implicit none
       integer, intent(in)  :: mype,iopt
       logical, intent(in)  :: lcc,lfc,lec,lnc
    end subroutine mpi_amr_restrict_fulltree
  end interface

  logical :: lfc, lec, lnc, saved_advanceAllLevels, change_advanceAllLevels
  integer :: tag

  call Timers_start("restrictAll")

  select case (gr_restrictAllMethod)
  case(1)
       call gr_restrictTree()
  case(2)
     if (nfacevar > 0) then
        lfc = .TRUE.
     else
        lfc = .FALSE.
     end if
     if (nvaredge > 0) then
        lec = .TRUE.
     else
        lec = .FALSE.
     end if
     if (nvarcorn > 0) then
        lnc = .TRUE.
     else
        lnc = .FALSE.
     end if

     tag = 123
#ifdef FLASH_GRID_PARAMESH4_0

     call mpi_morton_bnd_restrict(gr_meshMe, gr_meshNumProcs, .TRUE., lec, lnc, tag)
     call mpi_amr_restrict_fulltree(gr_meshMe,1,.true.,lfc,lec,lnc)
     call mpi_morton_bnd_restrict(gr_meshMe, gr_meshNumProcs, .FALSE., lec, lnc, tag)

#else

     change_advanceAllLevels = .NOT. advance_all_levels
     if (change_advanceAllLevels) then
        saved_advanceAllLevels = advance_all_levels
        advance_all_levels = .TRUE.
        call mpi_morton_bnd_restrict(gr_meshMe, gr_meshNumProcs, tag)
        advance_all_levels = saved_advanceAllLevels
     end if
     call mpi_amr_restrict_fulltree(gr_meshMe,1,.true.,lfc,lec,lnc)
     !  call mpi_amr_1blk_restrict(gr_meshMe,1,.true.,lfc,lec,lnc,.false.,.false.)
     if (change_advanceAllLevels) then
        call mpi_morton_bnd_restrict(gr_meshMe, gr_meshNumProcs, tag)
     end if

#endif

  case(3)
     call Grid_restrictByLevels( CENTER_FACES, -1, 1, checkFinestLevel=.TRUE.)
  end select

  call Timers_stop("restrictAll")
  
end subroutine Grid_restrictAllLevels
