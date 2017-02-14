!!****if* source/Grid/GridMain/Chombo/Grid_fillGuardCells
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
!!  The argument "gridDataStruct" can take on one of many valid 
!!  values to determine a specific grid data structure on which to apply
!!  the guardcell fill operation. The currently available options are listed with
!!  the arguments. Most users will use CENTER as the option, 
!!  since applications typically use the cell centered grid data, and they want
!!  guardcells to be filled for all the variables.
!!  More specialized applications, such as the unsplit methods, may want to use
!!  other options. 
!!  The user can also choose to fill guard cells either in a single direction,
!!  or all of them. 
!!
!!
!! ARGUMENTS 
!!  
!!
!!  gridDataStruct - integer constant, defined in "constants.h", 
!!                   indicating which grid data structure 
!!                   variable's guardcells to fill.
!!                   UG has 4 data structures for grid variables that
!!                   can have their guardcells filled. 
!!
!!                   unk                all cell centered variables in the grid
!!                   facex,facey,facez  all face centered variables along i,j,k 
!!                                      direction respectively
!!                   
!!                   valid values of gridDataStruct are  
!!                   CENTER             unk only
!!                   WORK               has no meaning in UG
!!                   The face variables are not yet implemented in UG
!!                   FACES              facex,facey, and facez
!!                   FACEX              facex
!!                   FACEY              facey
!!                   FACEZ              facez
!!                   CENTER_FACES     unk,facex,facey,facez
!!
!!  idir - direction of guardcell fill.  User can specify ALLDIR for all (x,y,z)
!!         directions, or if for example the algorithm only does one directional
!!         sweep at a time then time can be saved by filling only the guardcell
!!         direction that is needed.  A user would pass in the constants defined
!!         in constants.h IAXIS, JAXIS or KAXIS to fill guardcells in only one 
!!         direction.        
!!         All layers of guardcells in the given direction(s) are filled.
!!         In the current UG implementation, idir is ignored and a full
!!         guardcell fill in all directions is always performed.
!!
!!         THE REMAINING ARGUMENTS HAVE NO MEANING IN UG
!!
!!  minLayers - number of guardcell layers requested for all directions.
!!
!!   eosMode  - The mode in which eos is to be applied
!!   doEos    - the UG implementation does not act upon this argument
!!   maskSize - the size of the mask array. 
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
!! selectBlockType - IGNORED
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
!!      call Grid_fillGuardCells(CENTER, IAXIS)
!!
!!     This call will fill all guardcells for all cell-centered 
!!     variables in the x direction.
!!     
!! EXAMPLE 2
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells(CENTER_FACES, ALLDIR)
!!     
!!     This call fills guardcells along all directions in both
!!     cell centered and face centered data structures.
!!
!! NOTES
!!
!!   The masking functionality is not yet included in UG
!!  
!!***


#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine Grid_fillGuardCells(gridDataStruct,idir,minLayers,eosMode,doEos&
     &, maskSize, mask, makeMaskConsistent, doLogMask,selectBlockType, &
     unitReadsMeshDataOnly)
  use Grid_data, ONLY : gr_bndOrder, gr_allPeriodic, gr_justExchangedGC, &
       gr_eosModeNow, gr_blkList
#ifndef FLASH_GRID_UG
  use Grid_data, ONLY : gr_convertToConsvdInMeshInterp
#endif
  use gr_bcInterface, ONLY : gr_bcApplyToAllBlks
  use chombo_f_c_interface, ONLY : ch_fill_guardcells
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_guardCells
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_getListOfBlocks

  implicit none
  integer, intent(in) :: gridDataStruct
  integer, intent(in) :: idir
  integer, optional,intent(in) :: minLayers
  integer, optional,intent(in) :: eosMode
  logical, optional,intent(in) :: doEos
  integer, optional, intent(IN) :: maskSize
  logical, optional, dimension(:),intent(IN) :: mask
  logical, optional,intent(IN) :: makeMaskConsistent
  logical, optional,intent(IN) :: doLogMask  ! ignored in this implementation
  integer,optional,intent(in) :: selectBlockType
  logical, optional, intent(IN) :: unitReadsMeshDataOnly

  integer :: n, axis
  logical :: isWork=.false.
  integer :: listBlockType = ALL_BLKS
  integer,dimension(MDIM) :: layersArray
  integer :: gcEosMode, i
  integer :: numBlocks
  logical :: needEos

  if (gridDataStruct /= CENTER .and. gridDataStruct /= CENTER_FACES) then
     !DEV CD.  I am accepting CENTER_FACES for the time being because it
     !is passed by Grid_markRefineDerefine.  I do not support FACE variables
     !yet so CENTER_FACES is just CENTER for now.
     call Driver_abortFlash("[Grid_fillGuardCells]: Non-center not yet coded")
  end if

  if(present(eosMode)) then
     gcEosMode=eosMode
  else
     gcEosMode=gr_eosModeNow
  end if

  needEos=.true.
  layersArray = NGUARD


#ifndef FLASH_GRID_UG
  if((gridDataStruct==CENTER_FACES).or.(gridDataStruct==CENTER)) then
    call Grid_getListOfBlocks(listBlockType, gr_blkList, numBlocks)
    call gr_primitiveToConserve(gr_blkList,numBlocks)
  end if
#endif


  call ch_fill_guardcells()


  if(.not.gr_allPeriodic) then
     do n = 0,NDIM-1
        axis = gr_bndOrder(NDIM-n)
        call gr_bcApplyToAllBlks(axis,isWork)
     end do
  end if


#ifndef FLASH_GRID_UG
  if ((gridDataStruct==CENTER_FACES).or.(gridDataStruct==CENTER)) then
    call gr_conserveToPrimitive(gr_blkList,numBlocks, .TRUE.)
    if (gr_convertToConsvdInMeshInterp) then
      call gr_sanitizeDataAfterInterp(gr_blkList, numBlocks, 'after gc filling', layersArray)
     end if
  end if
#endif

  gr_justExchangedGC = .true.


  if(present(doEos)) then
     if(doEos.and.needEos) then
        call Timers_start("eos gc")
        call Grid_getListOfBlocks(listBlockType, gr_blkList, numBlocks)
        do i = 1,numBlocks
           call Eos_guardCells(gcEosMode, gr_blkList(i), corners=.true., &
                layers=layersArray)
        end do
        call Timers_stop("eos gc")
     end if
  end if

end subroutine Grid_fillGuardCells
