!!****if* source/Simulation/SimulationMain/incompFlow/Grid_updateRefinement
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
!!  consistent. The internal routine also calls gr_updateParticleRefinement to move the particles
!!  to the correct block after the grid refines.
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
  use Grid_interface, ONLY : Grid_markRefineDerefine,Grid_fillGuardcells
  use gr_interface, ONLY : gr_updateRefinement
  use paramesh_interfaces, ONLY : amr_restrict
  use Logfile_interface, ONLY : Logfile_open,Logfile_close
 
  use physicaldata, only : interp_mask_facex,interp_mask_facey,interp_mask_facez, &
                           interp_mask_facex_res,interp_mask_facey_res,interp_mask_facez_res  

  use IncompNS_data, ONLY : ins_predcorrflg,ins_outflowgridChanged

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
  
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
  !=============================================================================
  
  
  !Before we do any refinement, store a globalID made up of the block's cornerID,
  !refinement level and nodetype.  This is necessary to refine and derefine the 
  !particles.  If particles are not included, this routine will be a stub
  call gr_ptFillBlkParticleInfo()
  
  ins_predcorrflg = .false.
  ! We only consider refinements every nrefs timesteps.
  if (gr_nrefs .ne. 0) then
  if (mod(nstep, gr_nrefs) == 0) then
     
     call Timers_start("tree")  !1 of 2
     
     ! Make restriction of velocities and rhs linear (conservative)
     interp_mask_facex(VELC_FACE_VAR) = 2
     interp_mask_facex_res(VELC_FACE_VAR) = 1
     interp_mask_facey(VELC_FACE_VAR) = 2
     interp_mask_facey_res(VELC_FACE_VAR) = 1
     interp_mask_facez(VELC_FACE_VAR) = 2
     interp_mask_facez_res(VELC_FACE_VAR) = 1


     interp_mask_facex(RHDS_FACE_VAR) = 2
     interp_mask_facex_res(RHDS_FACE_VAR) = 1
     interp_mask_facey(RHDS_FACE_VAR) = 2
     interp_mask_facey_res(RHDS_FACE_VAR) = 1
     interp_mask_facez(RHDS_FACE_VAR) = 2
     interp_mask_facez_res(RHDS_FACE_VAR) = 1     

     iopt=1; iempty=1
     call amr_restrict(gr_meshMe,iopt,iempty)
     
     call Timers_start("markRefineDerefine")
     call Grid_markRefineDerefine()
     call Timers_stop("markRefineDerefine")
     
     call Timers_stop("tree")  !1 of 2 (We restart in gr_updateRefinement)
     !internal routine that does the actual amr refinement and
     !other housekeeping
     
     !write(*,*) 'Before Mark Refine derefine'  

     call gr_updateRefinement(gridChanged)

     !write(*,*) 'Before Guardcell fill in updateref'

     ! apply BC and fill guardcells for turbulent viscosity
     ins_outflowgridChanged = gridChanged
     gcMask = .TRUE.
     call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)
     ins_outflowgridChanged = .false.

     !write(*,*) 'Done with Guardcell fill in updateref'

  else
     if (present(gridChanged)) gridChanged = .FALSE.
     ins_outflowgridChanged = gridChanged
  end if
  end if 
 

  return
end subroutine Grid_updateRefinement
