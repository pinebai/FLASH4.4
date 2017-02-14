!!****f* source/Grid/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!!  This routine is normally called by the implementation of
!!  Grid_updateRefinement. It may also get called repeatedly
!!  during the initial construction of the Grid from
!!  Grid_initDomain.
!!
!! ARGUMENTS
!!
!!  none
!! 
!! SIDE EFFECTS
!!  
!!  This routine works by modifying the global (processor-local) flag
!!  arrays
!!      newchild, refine, derefine, stay
!!  imported from the AMR implementation.
!!
!! NOTES
!!  
!!  The default implementation of this routine uses the Flash runtime
!!  parameters refine_var_{1,2,3,4}, refine_cutoff_{1,2,3,4},
!!  derefine_cutoff_{1,2,3,4}, and refine_filter_{1,2,3,4} to determine
!!  refinement.
!!  
!!  A non-directional guardcell fill for all CENTER variables that are
!!  referenced by any of those refine_var_# runtime parameters
!!  (including EOS calls for guardcells, if any refinement variables
!!  refine_var_# require them to be current) is performed in the
!!  default implementation of this routine.  Since the checking of
!!  refinement criteria may reference the values of some variables in
!!  inactive blocks, the caller of Grid_markRefineDerefine should
!!  ensure that those inactive blocks have values that are up to date
!!  (at least for the variables to be used for refinement criteria);
!!  this is done by an explicit restriction call in the default APR
!!  implementation of Grid_updateRefinement.
!!  
!!  Users creating their own implementation of this interface or of
!!  Grid_updateRefinement should make sure that the above remains
!!  true; that is, they should include the appropriate call(s) to
!!  Grid_fillGuardCells (and to a restriction routine if necessary) in
!!  their implementation.
!!
!! SEE ALSO
!!
!!  Grid_updateRefinement
!!  Grid_initDomain
!!  Grid_fillGuardCells
!!***



subroutine Grid_markRefineDerefine()

implicit none
   
end subroutine Grid_markRefineDerefine














