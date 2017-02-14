!!****if* source/Grid/localAPI/gr_ensureValidNeighborInfo
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

  implicit none
  integer,intent(IN) :: requiredPattern


end subroutine gr_ensureValidNeighborInfo
