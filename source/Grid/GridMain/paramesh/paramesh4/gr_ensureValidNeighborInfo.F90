!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_ensureValidNeighborInfo
!!
!! NAME
!!
!!  gr_ensureValidNeighborInfo
!!
!! 
!! SYNOPSIS
!!
!!  call gr_ensureValidNeighborInfo(integer(IN) :: requiredPattern)
!!
!!
!! DESCRIPTION
!!
!!  This routine makes sure that some metadata kept internally by
!!  the grid implementation is in a state expected and relied upon
!!  by some code units (or Grid subunits).
!!
!!  Currently (2008-11-25) only Paramesh4 Grid implementations
!!  have an implementation of this subroutine that does anything.
!!
!!  This implementation checks whether the last call to
!!  mpi_amr_comm_setup was one that used the metadata exchange
!!  pattern used for guard cells.  If not, the last global
!!  data exchange operation under mpi_amr_comm_setup must have
!!  been for prolongation or restriction or flux or edge data
!!  correction.  In that case call mpi_amr_comm_setup as if
!!  for guard cell exchange, but don't actually communicate any
!!  solution variable data, don't copy anythign into any blocks'
!!  guard cells, don't interpolate, etc.  We COULD, however,
!!  call Grid_fillGuardCells here.
!!
!!  The purpose here is to make sure the paramesh metadatastructures
!!  surr_blks, parent, child, etc. are in the state that is "normal"
!!  and expected by some parts of Grid code, in particular code
!!  in GridParticles.  This "normal" state is the one that obtains
!!  after Grid_fillGuardCells has just been called.
!!
!! ARGUMENTS
!!
!!  requiredPattern - either 0 to mean ANY, i.e., any valid pattern is
!!                    acceptable, or a specific recognized pattern ID value
!!                    to indicate that the caller requires valid data
!!                    exchanged according to that pattern.  This is
!!                    specific to a Grid implementation, see NOTES below.
!!                    At the moment, only 0 and 10 are implemented,
!!                    others act like 0.
!!
!! NOTES
!!
!!  Paramesh4 uses/sets the following pattern ID values:
!!   10  is the comm pattern for guard cell filling,
!!   20  is the comm pattern for prolongation to new leaf blocks,
!!   30  is the comm pattern for flux correction (incl. edge if facevars present?),
!!   40  is the comm pattern for restriction (non-fulltree; fulltree seems to use 10).
!!
!! HISTORY
!!   2008-11-19 created in trunk - kw
!!   2008-11-25 added requiredPattern argument and associated logic - kw
!!***

#include "constants.h"

subroutine gr_ensureValidNeighborInfo(requiredPattern)
  use Driver_interface, only: Driver_abortFlash
  use Logfile_interface, only: Logfile_stamp
  use Timers_interface, only: Timers_start, Timers_stop

  use Grid_data, only: gr_meshMe
  use tree, only: grid_analysed_mpi
  use physicaldata, only: mpi_pattern_id

  implicit none
  integer,intent(IN) :: requiredPattern


  if (grid_analysed_mpi .NE. 1) then
     ! This routine should not get called in a state were Paramesh
     ! is not initialized at all.
     if (gr_meshMe==MASTER_PE) print*,'gr_ensureValidNeighborInfo: grid_analysed_mpi is', grid_analysed_mpi
     call Driver_abortFlash("gr_ensureValidNeighborInfo found grid_analysed_mpi different from 1!")
  endif


  if ((requiredPattern>0 .AND. mpi_pattern_id .NE. requiredPattern)) then
     call Logfile_stamp( mpi_pattern_id, "[gr_ensureValidNeighborInfo] found mpi_pattern_id")
     call Timers_start('no-data comm_setup')
     select case (requiredPattern)
     case(10)
        call gr_commSetUp(-1)
     case default
        call Driver_abortFlash("gr_ensureValidNeighborInfo: requiredPattern not implemented.")   
     end select
     call Timers_stop('no-data comm_setup')
  else if ((requiredPattern>0 .AND. mpi_pattern_id .EQ. requiredPattern)) then
     ! do nothing, not even log
  else                          ! i.e., requiredPattern .le. 0
     select case (mpi_pattern_id)
     case(10)
           ! do nothing, not even log
     case(20,30,40)             ! log the pattern in effect, leave unchanged
        call Logfile_stamp( mpi_pattern_id, "[gr_ensureValidNeighborInfo] found mpi_pattern_id, leaving unchanged")
     case default
        call Logfile_stamp( mpi_pattern_id, "[gr_ensureValidNeighborInfo] found mpi_pattern_id")
        call Timers_start('no-data comm_setup')
        call gr_commSetUp(-1)
        call Timers_stop('no-data comm_setup')
     end select
  end if



end subroutine gr_ensureValidNeighborInfo
