!!****f* source/Grid/Grid_fillGuardCells
!!
!! NAME
!!  Grid_fillGuardCells
!!
!! SYNOPSIS
!!
!!  call Grid_fillGuardCells(integer(IN) :: gridDataStruct,
!!                           integer(IN) :: idir,
!!                  optional,integer(IN) :: minLayers,
!!                  optional,integer(IN) :: eosMode,
!!                  optional,logical(IN) :: doEos,
!!                  optional,integer(IN) :: maskSize,
!!                  optional,logical(IN) :: mask(maskSize),
!!                  optional,logical(IN) :: makeMaskConsistent,
!!                  optional,logical(IN) :: doLogMask,
!!                  optional,integer(IN) :: selectBlockType,
!!                  optional,logical(IN) :: unitReadsMeshDataOnly)
!!
!!
!! DESCRIPTION 
!!
!!  This routine fills the guardcells of the specified data structure of 
!!  each block present in the mesh. If the adjacent blocks are at the same
!!  level of refinement, then the values are copied from the appropriate
!!  internal cells of the neighboring block. If the blocks are at different
!!  levels of refinement, then in addition to copying, there will be 
!!  restriction or prolongation. 
!!
!!  The optional arguments related to the Eos calls are provided because if
!!  a solver relies on a call to Eos to establish thermodynamic equilibrium
!!  in various variables, then either the solver can call Eos explicitly itself,
!!  or more conveniently, it can instruct the gcfill process to do it. This 
!!  feature is especially useful when mesh quantities are used to derive 
!!  for example tracer particle attributes. 
!!  
!!  For performance reason, users may choose to use masking in filling the 
!!  guardcells, which is feature supported only with Paramesh 4. They can 
!!  do so by using the optional arguments mask, maskSize and makeMaskConsistent. However, 
!!  users should exercise extreme caution in using masking. Please see the 
!!  NOTES section for potential side effects. If masking is present, and 
!!  both makeMaskConsistent, and doEos are true, the Eos calculation is done only
!!  if one of the variables selected in the mask is affected by Eos calculation.
!!  A local routine gr_makeMaskConsistent handles this processing.
!!
!!  The argument "gridDataStruct" can take on one of many valid 
!!  values to determine a specific grid data structure on which to apply
!!  the guardcell fill operation. The currently available options are listed with
!!  the arguments. Most users will use CENTER as the option,
!!  since applications typically use the cell centered grid data, and they want
!!  guardcells to be filled for all the variables. When using AMR with paramesh,
!!  the function Grid_markRefineDerefine may use the single variable data
!!  structure provided by paramesh represented by the option "WORK".
!!  More specialized applications, such as the unsplit methods, may want to use
!!  other options. 
!!  The user can also choose to fill guard cells either in a single direction,
!!  or all of them. For most of the Flash solvers supplied with the release,
!!  guard cell are filled in all directions.
!!
!!  
!!
!!
!! ARGUMENTS 
!!  
!!
!!  gridDataStruct - integer constant, defined in "constants.h", 
!!                   indicating which grid data structure 
!!                   variable's guardcells to fill.
!!                   Paramesh has 5 data structures for grid variables, the first
!!                   four include all physical variables defined on the mesh. The 
!!                   fifth one includes a single variable.
!!
!!                   unk                cell centered, 
!!                   facex,facey,facez  face centered along i,j,k 
!!                                      direction respectively
!!                   work               cell centered, single variable.
!!                   
!!                   valid values of gridDataStruct are  
!!                   CENTER             unk only
!!                   WORK               work 
!!                   FACEX              facex
!!                   FACEY              facey
!!                   FACEZ              facez
!!                   FACES              facex,facey, and facez
!!                   CENTER_FACES     unk,facex,facey,facez
!!
!!  idir - direction of guardcell fill.  User can specify ALLDIR for all (x,y,z)
!!         directions, or if for example the algorithm only does one directional
!!         sweep at a time then time can be saved by filling only the guardcell
!!         direction that is needed.  A user would pass in the constants defined
!!         in constants.h IAXIS, JAXIS or KAXIS to fill guardcells in only one 
!!         direction.        
!!         All layers of guardcells in the given direction(s) are filled if one
!!         of IAXIS, JAXIS, or KAXIS is specified, or if ALLDIR is specified
!!         and minLayers (see below) is not present.
!!
!!  minLayers - minimum number of guardcell layers requested for all directions.
!!              The caller requests at least this many layers of
!!              guardcells to be filled.  If idir is given as one of IAXIS,
!!              JAXIS, or KAXIS, this applies to any directions NOT selected by
!!              idir. On the other hand, if idir is given as ALLDIR, then
!!              minLayers appliers to all directions equally.
!!              If not specified, the default is 0 if idir is given as IAXIS,
!!              JAXIS, or KAXIS, meaning no guardcells need to be filled in for
!!              the directions perpendicular to what idir selects; the default
!!              is NGUARD, the full number of available guardcells, if idir
!!              is ALLDIR.
!!              Note that the caller can specify, using minLayers, how many
!!              layers it needs filled, but the implementation may do
!!              more and actually fill all or some additional layers.
!!              The caller must therefore not rely on some guardcells
!!              remaining unchanged.
!!              
!!              This argument is meaningful with PARAMESH4 only.
!!
!!   eosMode -  the EOS mode being used by the solver that is calling the 
!!              routine. The valid values are :
!!              MODE_DEFAULT     the default eos mode being used by the grid unit
!!              MODE_DENS_EI     density and energy input, pressure and temp output
!!              MODE_DENS_PRES   density/pressure input, temperature/energy output
!!              MODE_DENS_TEMP   density/temperature input, pressure/energy output
!!
!!  doEos    - a logical variable indicating if the calling routine wants the
!!             gcfill process to also make sure that Eos is applied to achieve
!!             thermodynamically consistent values of all variables.
!! 
!!             The mask related arguments are useful only with PARAMESH4. They have 
!!             no meaning in PARAMESH2 and UG.
!!  maskSize - the size of the mask array. 
!! 
!!  mask -  It is a one-dimensional logical array 
!!          with indices corresponding to variables in the grid data
!!          structures. If a variable should have its guardcells filled,
!!          the corresponding element in "mask" is true, otherwise it is
!!          false.
!!          The mask is always ignored if the runtime parameter
!!          enableMaskedGCFill is set .FALSE.
!!  
!! makeMaskConsistent - If true when mask is applied, it is made sure that for
!!          all the selected variables in the mask, the ones they are dependent
!!          on are true too. It is also determined whether there is a need to 
!!          apply Eos if doEos argument is true.
!!
!! doLogMask - If present and true when mask is applied, stamp some information
!!          about masking (and whether Eos is or would be called) to the Logfile.
!!          Ignored in implementations that do not support variable masking.
!!
!!  selectBlockType - selects the blocks whose guard cells must be filled
!!            by block type.
!!
!!              This argument is ignored in UG Implementations.
!!
!!              For PARAMESH Grid implementations, recognized values are :
!!
!!              ALL_BLKS    all local blocks on a processor.
!!                          This is not valid in all PARAMESH versions.
!!
!!              ACTIVE_BLKS all currently active blocks, in paramesh
!!              context that means parent and leaf blocks
!!              
!!              LEAF        only LEAF blocks, i.e., bocks whose paramesh
!!                          node type is 1
!!            
!!              These constants are defined in constants.h
!!
!!       Note that if advance_all_levels is set (in a PARAMESH version
!!       that implements this global flag), guard cells in all levels of
!!       blocks are filled anyway.
!!
!!       Note that the caller can specify, using selectBlockType, which
!!       blocks it needs filled, but the implementation may do
!!       more and actually fill guard cells in additional blocks.
!!       The caller must therefore not rely on guardcells
!!       in other blocks remaining unchanged.
!!
!! unitReadsMeshDataOnly - specifies that the unit calling Grid_fillGuardCells
!!                         does not update any internal grid data.  This
!!                         allows us to skip the next guard cell fill because
!!                         the guard cells already contain up to date data.
!!
!! EXAMPLE 
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells( CENTER, IAXIS)
!!
!!     This call will fill all guardcells for all cell-centered 
!!     variables in the x direction.
!!     
!! EXAMPLE 2
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells( WORK, ALLDIR)
!!     
!!     This call fills guardcells along all directions. The operation is applied
!!     to the WORK data structure available in paramesh only.
!!
!! SIDE EFFECTS
!!
!!  After this function returns, all parents of leaf blocks will have current and 
!!  valid solution data (at least for the variables determined by the gridDataStruct
!!  and mask dummy arguments). This is because amr_guardcell calls amr_restrict
!!  internally.  If selectBlockType is used, even more parents of blocks of interest
!!  may get updated by restriction.
!!
!! NOTES
!!
!!  In the default mode, this routine fills guardcells of all the variables
!!  in the specified data structure. However, if masking is used
!!  it is possible to fill guardcells of only 
!!  selected variables. The users must exercise a great deal of caution in 
!!  using masking in the guardcell filling, since masking out some variables may
!!  have unintended consequences.  For example if any conserved variable
!!  is set to true in the mask, density must be true, otherwise the interpolated
!!  values of the the conserved variable at fine-coarse boundaries will be 
!!  wrong. It is highly recommended that the argument makeMaskConsistent be set to 
!!  true when using masking, since that provides error checking for
!!  the variable sets included in the released solvers. However,  
!!  these checks may not be enough for new solvers introduced by the users.
!!
!!  This function, or one of the lower-level functions invoked by it, MUST
!!  be called (with updated solution data) before child blocks may be removed
!!  in the course of derefinement! See side effects, above, for the reason.
!!
!! 
!!***



subroutine Grid_fillGuardCells( gridDataStruct, idir,&
     minLayers, eosMode, doEos, maskSize, mask, makeMaskConsistent,&
     doLogMask,selectBlockType,unitReadsMeshDataOnly)
  
  implicit none
  integer, intent(in) :: gridDataStruct
  integer, intent(in) :: idir
  integer, optional,intent(in) :: minLayers
  integer, optional, intent(in) :: eosMode
  logical, optional,intent(IN) :: doEos
  integer, optional, intent(IN) :: maskSize
  logical, optional, dimension(:),intent(IN) :: mask
  logical, optional,intent(IN) :: makeMaskConsistent
  logical, optional,intent(IN) :: doLogMask
  integer, optional,intent(IN) :: selectBlockType
  logical, optional, intent(IN) :: unitReadsMeshDataOnly
end subroutine Grid_fillGuardCells
