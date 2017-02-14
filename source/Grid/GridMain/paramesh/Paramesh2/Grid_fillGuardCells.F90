!!****if* source/Grid/GridMain/paramesh/Paramesh2/Grid_fillGuardCells
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
!!  The argument "gridDataStruct" can take on one of many valid 
!!  values to determine a specific grid data structure on which to apply
!!  the guardcell fill operation. The currently available options are listed with
!!  the arguments. Most users will use CENTER as the option,
!!  since applications typically use the cell centered grid data, and they want
!!  guardcells to be filled for all the variables.
!!  More specialized applications, such as the unsplit methods, may want to use
!!  other options. 
!!  The user can also choose to fill guard cells either in a single direction,
!!  or all of them. For most of the Flash solvers supplied with the release,
!!  guard cell are filled in all directions.
!!
!!  The optional arguments related to the Eos calls are provided because if
!!  a solver relies on a call to Eos to establish thermodynamic equilibrium
!!  in various variables, then either the solver can call Eos explicitly itself,
!!  or more conveniently, it can instruct the gcfill process to do it. This 
!!  feature is especially useful when mesh quantities are used to derive 
!!  for example tracer particle attributes. The the user opts to enable doEos,
!!  they will not have to determine whether Eos call is needed, the code will
!!  do it for them
!!  
!!
!!
!! ARGUMENTS 
!!  
!!
!!  gridDataStruct - integer constant, defined in "constants.h", 
!!                   indicating which grid data structure 
!!                   variable's guardcells to fill.
!!                   PARAMESH2 supports 2 data structures for grid variables, 
!!                   the first one "unk "includes all physical variables 
!!                   defined, and "work" which includes a single variable.
!!
!!                   unk                cell centered, 
!!                   facex,facey,facez  not defined for PARAMESH2
!!                   work               cell centered, single variable.
!!
!!                   For PARAMESH2, valid values of gridDataStruct are  
!!
!!                   CENTER             unk only
!!                   WORK               work 
!!                   CENTER_FACES     unk,facex,facey,facez
!!                   In the PARAMESH2 implementation, CENTER_FACES is the same
!!                   as CENTER, and if gridDataStruct is FACES or FACEX/Y/Z the
!!                   code aborts.
!!
!!  idir - direction of guardcell fill.  User can specify ALLDIR for all (x,y,z)
!!         directions, or if for example the algorithm only does one directional
!!         sweep at a time then time can be saved by filling only the guardcell
!!         direction that is needed.  A user would pass in the constants defined
!!         in constants.h IAXIS, JAXIS or KAXIS to fill guardcells in only one 
!!         direction.        
!!         All layers of guardcells in the given direction(s) are filled.
!!
!!  minLayers - number of guardcell layers requested for all directions.
!!              The caller requests at least this many layers of
!!              guardcells to be filled.  This applies to any
!!              directions NOT selected by idir, so it is ignored if
!!              idir is given as ALLDIR.
!!              If not specified, the default is 0 meaning no guardcells
!!              need to be filled in for the directions perpendicular to
!!              what idir selects.
!!              Note that the caller can specify, using minLayers, how many
!!              layers it needs filled, but the implementation may do
!!              more and actually fill all or some additional layers.
!!              The caller must therefore not rely on some guardcells
!!              remaining unchanged.
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
!!            The next three arguments have no meaning in Paramesh 2
!!
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
!!  dummy argument). This is because amr_guardcell calls amr_restrict
!!  internally.  If selectBlockType is used, even more parents of blocks of interest
!!  may get updated by restriction.
!!
!! NOTES
!!
!!  This function, or one of the lower-level functions invoked by it, MUST
!!  be called (with updated solution data) before child blocks may be removed
!!  in the course of derefinement! See side effects, above, for the reason.
!!  
!!  When using AMR with paramesh, the function Grid_markRefineDerefine
!!  used to use the single variable data structure provided by paramesh
!!  represented by the option "WORK". That is not the case any more in
!!  the current FLASH3 default implementation when using Paramesh2.
!!
!!***

#ifndef DEBUG_GRID
#define DEBUG_GRID
#endif

subroutine Grid_fillGuardCells(gridDataStruct, idir,&
     minLayers,eosMode,doEos, maskSize, mask, makeMaskConsistent,&
     doLogMask,selectBlockType,unitReadsMeshDataOnly)

  use Grid_data, ONLY : gr_blkList, gr_justExchangedGC,gr_eosModeNow, gr_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_getListOfBlocks
  use paramesh_interfaces, ONLY : amr_guardcell
  use tree, ONLY : lnblocks
  use Eos_interface, ONLY : Eos_guardCells
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: gridDataStruct
  integer, intent(in) :: idir
  integer, optional,intent(in) :: minLayers
  integer, optional,intent(in) :: eosMode
  logical,optional,intent(in) :: doEos
  integer, optional, intent(IN) :: maskSize
  logical, optional, dimension(:), intent(IN) :: mask
  logical, optional, intent(IN) :: makeMaskConsistent
  logical, optional, intent(IN) :: doLogMask  ! ignored in this implementation
  integer, optional, intent(IN) :: selectBlockType
  logical, optional, intent(IN) :: unitReadsMeshDataOnly

  logical,save :: facesWarningDone=.FALSE.
  
  integer :: iopt,guard, numLeafBlocks,i, loc_idir, maxNodetype_gcWanted
  integer :: ierr,gcEosMode

#ifdef DEBUG_GRID
  logical:: validDataStructure
  validDataStructure = (gridDataStruct==CENTER).or.&
                       (gridDataStruct==WORK).or.&
                       (gridDataStruct==CENTER_FACES)
  if(.not.validDataStructure)then
     call Driver_abortFlash("GCfill: invalid data structure")
  end if
#endif

  !if guardcells were just exchanged, don't do it again!
!!  if (gr_justExchangedGC) return

#if (0)
! This code SHOULD be used, were it not for a failure of Paramesh2
! to work correctly with idir .ne. 0.
  loc_idir = idir
  ! Request guardcells for all directions if idir specifies directional
  ! guardcell fill and any guardcells for additional directions are also
  ! requested. - KW
  if (present(minLayers)) then
     if (minLayers .gt. 0) loc_idir = 0
  end if
  if (loc_idir .eq. NODIR) return
#endif

  call Timers_start("guardcell Barrier")
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call Timers_stop("guardcell Barrier")
  call Timers_start("guardcell internal")

#if (1)
! This is the workaround for Paramesh2's failure
! to work correctly with idir .ne. 0.
  loc_idir = 0
#endif

  guard = NGUARD

  if ((gridDataStruct==CENTER).or.(gridDataStruct==CENTER_FACES)) then
     iopt = 1
     if(gridDataStruct==CENTER_FACES) then
	if (.NOT. facesWarningDone) then
           if (gr_meshMe==MASTER_PE) then
              print*,"Warning : PM2 supports only cell centered data, "
              print*,"          CENTER_FACES reverts to CENTER in GCfill"
           endif
           facesWarningDone = .TRUE.
        endif
     end if
  elseif(gridDataStruct==WORK) then
     iopt = 2
  else
     iopt = 0
  end if

  if (present(selectBlockType)) then
     select case (selectBlockType)
     case(LEAF)
        maxNodetype_gcWanted = 1
     case(ACTIVE_BLKS)
        maxNodetype_gcWanted = 2
     case(ALL_BLKS)
        maxNodetype_gcWanted = 3
#ifdef DEBUG_GRID
     case default
        call Driver_abortFlash('Grid_fillGuardCells: unrecognized value for selectBlockType!')
#endif
     end select
  else
     maxNodetype_gcWanted = 1 !by default, guard cells are wanted for leaf blocks only!?
  end if

  if ((gridDataStruct==CENTER).or.(gridDataStruct==CENTER_FACES)) then
     call Grid_getListOfBlocks(ACTIVE_BLKS, gr_blkList, numLeafBlocks)
     call gr_primitiveToConserve(gr_blkList,numLeafBlocks)
     call amr_guardcell(gr_meshMe, iopt, guard, loc_idir, maxNodetype_gcWanted)
     call gr_conserveToPrimitive(gr_blkList,numLeafBlocks, .TRUE.)

  else
     call amr_guardcell(gr_meshMe, iopt, guard, loc_idir, maxNodetype_gcWanted)
  end if
  if(present(doEos))then
     if (doEos) then
        if(present(eosMode)) then
           gcEosMode=eosMode
        else
           gcEosMode=gr_eosModeNow
        end if

        if (maxNodetype_gcWanted .EQ. 1) then
           call Grid_getListOfBlocks(LEAF, gr_blkList, numLeafBlocks)
        else if (maxNodetype_gcWanted .EQ. 2) then
           if ((gridDataStruct.NE.CENTER).AND.(gridDataStruct.NE.CENTER_FACES)) &
                  call Grid_getListOfBlocks(ACTIVE_BLKS, gr_blkList, numLeafBlocks)
        else if (maxNodetype_gcWanted .GE. 3) then
           call Grid_getListOfBlocks(ALL_BLKS, gr_blkList, numLeafBlocks)
        end if
        do i=1,numLeafBlocks
           call Eos_guardCells(gcEosMode,gr_blkList(i),.true., skipSrl=.TRUE.)
        end do
     end if
  end if

  gr_justExchangedGC = .true.
  call Timers_stop("guardcell internal")


end subroutine Grid_fillGuardCells
