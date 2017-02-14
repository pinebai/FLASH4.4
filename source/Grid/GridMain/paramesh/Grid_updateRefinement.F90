!!****if* source/Grid/GridMain/paramesh/Grid_updateRefinement
!!
!! NAME
!!
!!  Grid_updateRefinement
!!
!!
!! SYNOPSIS
!!  
!!  call Grid_updateRefinement(integer(IN) :: nstep,
!!                             real(IN)    :: time,
!!                    OPTIONAL,logical(OUT):: gridChanged)
!!
!!
!! DESCRIPTION
!!
!!  Apply the user-defined refinment critera to determine which blocks need 
!!  to be refined and derefined.  Once the blocks are marked, call 
!!  gr_updateGridRefinement to carry out the rest of the routine.
!!  The internal routine does the refinements (amr_refine_derefine)
!!  During this
!!  stage, the blocks are redistributed across processors (if needed).  
!!  After the refinement, the newly created child blocks are filled via
!!  prolongation from the coarse parents.  This prolongation step can use
!!  prolongation routines from paramesh, or defined by the user
!!  Once the prolongation is done, the guardcells are filled.  Finally, the
!!  EOS is called on the block interiors to make them thermodynamically
!!  consistent. The internal routine also calls Particles_updateRefinement to
!!  move the particles to the correct block after the grid refines.
!!
!!
!!
!! ARGUMENTS
!!
!!  nstep : current step number
!!  time  : current evolution time
!!  gridChanged : returns TRUE if grid may actually have changed.
!!
!!***


subroutine Grid_updateRefinement( nstep,time, gridChanged)

  use Grid_data, ONLY : gr_nrefs, gr_meshMe
!!$  use Grid_data, ONLY : gr_blkList, gr_eosMode
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_markRefineDerefine
  use gr_interface, ONLY : gr_updateRefinement
  use paramesh_interfaces, ONLY : amr_restrict
  use Logfile_interface, ONLY : Logfile_open,Logfile_close
 
!!$  use tree, ONLY : newchild, lnblocks
!!$  use paramesh_interfaces, ONLY : amr_refine_derefine, &
!!$                                  amr_prolong
  
  implicit none


#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: nstep
  real, intent(in) :: time
  logical,intent(out),OPTIONAL :: gridChanged

!!$  integer :: i

!!$  integer :: count

  integer :: iopt,iempty, logUnit
  logical :: logUnitLocal=.true.
  
  !=============================================================================
  
  
  !Before we do any refinement, store a globalID made up of the block's cornerID,
  !refinement level and nodetype.  This is necessary to refine and derefine the 
  !particles.  If particles are not included, this routine will be a stub
  call gr_ptFillBlkParticleInfo()
  
  
  ! We only consider refinements every nrefs timesteps.
  if (mod(nstep, gr_nrefs) == 0) then
     
     call Timers_start("tree")  !1 of 2
     
     iopt=1; iempty=1
     call amr_restrict(gr_meshMe,iopt,iempty)
     
     call Timers_start("markRefineDerefine")
     call Grid_markRefineDerefine()
     call Timers_stop("markRefineDerefine")
     
     call Timers_stop("tree")  !1 of 2 (We restart in gr_updateRefinement)
     !internal routine that does the actual amr refinement and
     !other housekeeping
     
     call gr_updateRefinement(gridChanged)


  else
     if (present(gridChanged)) gridChanged = .FALSE.
  end if
  
  return
end subroutine Grid_updateRefinement
