!!****if* source/Grid/GridMain/paramesh/Grid_restrictByLevels
!!
!! NAME
!!  Grid_restrictByLevels
!!
!! SYNOPSIS
!!
!!  call Grid_restrictByLevels(integer(IN) :: gridDataStruct, 
!!                      integer(IN) :: fromLevel,
!!                      integer(IN) :: toLevel,
!!             optional,logical(IN) :: checkFinestLevel,
!!             optional,integer(IN) :: maskSize,
!!             optional,logical(IN) :: mask(*))
!!  
!! DESCRIPTION 
!!
!!  Update solution data by a restriction sweep from blocks at refinement
!!  level given by fromLevel to refinement level toLevel.
!!
!!  The argument "gridDataStruct" can take on one of several valid 
!!  values to determine a specific grid data structure on which to apply
!!  the data restriction operation. The currently available options are listed
!!  with the arguments. User code will probably use either CENTER as the option,
!!  with a mask selecting one or a few variables, or use WORK.
!!
!!
!!  
!!
!! ARGUMENTS 
!!  
!!
!!  gridDataStruct - integer constant, defined in "constants.h", 
!!                   indicating which grid data structure 
!!                   variable's cell data to update.
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
!!                   WORK               work  !! valid only for Paramesh
!!                   May also be supported in PM3/PM4:
!!                   FACES              facex,facey and facez 
!!                   FACEX              facex
!!                   FACEY              facey
!!                   FACEZ              facez
!!                   CENTER_FACES       unk as well as facex,facey, and facez
!!                                      in current implementation of PM2 and
!!                                      UG, this is equivalent to CENTER
!!
!!  fromLevel - refinement level where the sweep should start, the finest level.
!!              It can be given as
!!               - a positive integer, meaning an absolute refinement level
!!                 as in the runtime parameters lrefine_min, lrefine_max, etc.,
!!                 i.e., 1 means the refinement level of the root block(s).
!!               - 0, meaning the same as specifying lrefine_max
!!               - a negative integer, meaning a refinement level relative
!!                 to the highest refinement level that is currently realized
!!                 anywhere on the grid; the highest realized level is -1,
!!                 the next coarse level is -2, etc. 
!!  toLevel -   refinement level where the sweep should end, the coarsest level.
!!              Can be given as a positive integer, 0, or a negative integer,
!!              where the meaning is as described for forLevel.
!!              If toLevel indicates a level equal to or higher (i.e., finer)
!!              than fromLevel, no restriction will take place.
!!
!!  checkFinestLevel - whether to check for the finest level that is realized
!!              on the grid by an initial MPI_ALLREDUCE operation.
!!              Default is true.  If set to false, the code may unnecessarily
!!              sweep over nonexistent levels.
!!              Note the this check will be done anyway if a negative
!!              fromLevel or toLevel is specified.
!!
!!
!!  maskSize - the size of the mask array. It is an optional argument.
!! 
!!  mask -  This is an optional argument. If this argument is missing, 
!!          the routine assumes that all variables in the specified 
!!          data structure should get updated.
!!          If it is present, it is a one-dimensional logical array 
!!          with indices correspoding to variables in the grid data
!!          structures. If a variable should have its cells updated,
!!          the corresponding element in "mask" is true, otherwise it is
!!          false.
!! 
!!
!! EXAMPLE 
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      call Grid_restrictByLevels( CENTER, 5, 3)
!!
!!     This call will update UNK data in interior cells of all non-LEAF blocks
!!     at refinement levels 4 and 3, by a restriction sweep based on LEAF
!!     blocks at refinement levels 5 and 4 in regions of the domain where
!!     such blocks exists and based on the data in non-LEAF blocks at level 5
!!     elsewhere.
!!     
!! EXAMPLE 2
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      call Grid_restrictByLevels( WORK, 5, 3)
!!
!!     This call ...
!!     
!! EXAMPLE 3
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      logical :: restrictMask(NUNK_VARS)
!!
!!      restrictMask = .FALSE.
!!      restrictMask(DENS_VAR) = .TRUE.
!!      call Grid_restrictByLevels( CENTER, 5, 3, maskSize=NUNK_VARS, mask=restrictMask)
!!
!!     This call ...
!!
!! SIDE EFFECTS
!!
!!  After this function returns, all parent and ancestor blocks at refinement levels
!!  between fromLevel (exclusive) and toLevel (inclusive) will have updated
!!  solution data for the variables determined by the gridDataStruct
!!  and mask dummy arguments, based on restriction operations beginning with blocks
!!  at either a refinement level of fromLevel or coarser LEAF blocks (in regions where
!!  the LEAF blocks are less refined than fromLevel) and sweeping up the grid tree.
!!
!!  Data in variables not selected by gridDataStruct and mask will be either unchanged 
!!  (in the PARAMESH Grid [as well as UG]) or may also have undergone restriction
!!  (in Grid implementations where masks are ignored).
!!
!! NOTES
!!
!!  Face variables not tested / not expected to work.
!!
!!  Restriction, as we use the term, refers only to changing the solution data on
!!  blocks, while leaving the structure of the Grid (existence of blocks, AMR tree 
!!  structure) unchanged.  In particular, no blocks are removed by "restriction".
!!
!!  Guard cells are not updated.
!!
!!  With a uniform grid, all calls are no-ops.
!!
!!  Constants CENTER, WORK, etc. are defined in constants.h .
!!  
!!***

subroutine Grid_restrictByLevels( gridDataStruct, fromLevel, toLevel, checkFinestLevel,&
  maskSize,mask)

  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

  integer, intent(in) :: gridDataStruct
  integer, intent(in) :: fromLevel, toLevel
  logical, optional,intent(in) :: checkFinestLevel
  integer, optional,intent(in) :: maskSize
  logical,dimension(*),optional,intent(in) :: mask

  call Driver_abortFlash('Grid_restrictByLevel is currently not implemented &
       &for this version of the PARAMESH Grid.')
end subroutine Grid_restrictByLevels
