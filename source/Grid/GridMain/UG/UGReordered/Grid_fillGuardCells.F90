!!****if* source/Grid/GridMain/UG/UGReordered/Grid_fillGuardCells
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
!!  Fills guard cells for all blocks.  
!! 
!!  This routine has a flexible interface and offers
!!  the user a number of options depending on the arguments passed in the interface.  
!!  The user can choose to fill guardcells for all
!!  variables or some of the  variables, fill guardcells in all directions
!!  or just a single direction and fill guardcells on a given level (to be implemented).
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
!!               unk                cell centered, 
!!               facex,facey,facez  face centered along i,j,k direction respectively
!!               work               cell centered, single variable.
!!                   Uniform grid supports the first four data structure, there is 
!!                   no "work" data structure in UG
!!                   With those datastructures, the ones included for specific 
!!                   values of gridDataStruct are  
!!
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
!!         In this current UGReordered implementation, idir is ignored and a
!!         full guardcell fill in all directions is always performed, except
!!         for boundary conditions at the edges of the physical domain.
!!
!!         THE REMAINING ARGUMENTS HAVE NO MEANING IN UG
!!
!!  minLayers - number of guardcell layers requested for all directions.
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
!!              In the current UG implementation, minLayers is ignored
!!              and a full guardcell fill in all directions is always
!!              performed, except for bondary conditions at the edges of
!!              the physical domain, where the code currently always behaves
!!              as if minLayers=0 had been specified.
!! 
!!  eosMode -   the EOS mode being used by the solver that is calling the 
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
!!     This call fills guardcells along all directions and the
!!     operation is applied to the WORK data structure available in
!!     paramesh only.  Thus this call is not supported in this
!!     implementation.
!!
!! NOTES
!!
!!   The masking functionality is not yet included in UG
!!  
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif


subroutine Grid_fillGuardCells( gridDataStruct,idir,minLayers,eosMode,doEos&
     &, maskSize, mask, makeMaskConsistent, doLogMask,selectBlockType,&
     unitReadsMeshDataOnly)

  use Grid_interface, ONLY : Grid_getBlkIndexLimits
  use gr_bcInterface, ONLY : gr_bcApplyToAllBlks
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_axisComm, gr_exch, &
        gr_numDataStruct, gr_gridDataStruct, &
        gr_offset, gr_justExchangedGC,gr_bndOrder,gr_allPeriodic, gr_meshMe

 implicit none
#include "constants.h"
#include "Flash.h"

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

  integer,dimension(MDIM) :: guard
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: blockID = 1
  integer :: i,j,k,beginDataType,endDataType,n
  logical :: isWork=.false.
  integer :: axis, localGridDataStruct

  !These are arrays to hold the sending start point indicies and receiving
  !start point indicies for the MPI data types in the shift data subroutine
  !They are used to specify the starting incidies in the unk array

  integer, dimension(MDIM, MDIM+1) :: sendRight, sendLeft, recvRight, recvLeft

  !if we have just exchanged guardcells then we don't need to do in again.
  !FOR FUTURE: not yet implemented gr_justExchanged
!!$  if (gr_justExchanged) then
!!$     return
!!$  end if


#ifdef DEBUG_GRID
  if((gr_domainBC(LOW,IAXIS) == PERIODIC).and.(gr_domainBC(HIGH,IAXIS) /= PERIODIC)) then
     call Driver_abortFlash("Gaurd Cell fill : one edge is periodic, one is not")
  end if

  if((gr_domainBC(LOW,JAXIS) == PERIODIC).and.(gr_domainBC(HIGH,JAXIS) /= PERIODIC)) then
     call Driver_abortFlash("Gaurd Cell fill : one edge is periodic, one is not")
  end if

  if((gr_domainBC(LOW,KAXIS) == PERIODIC).and.(gr_domainBC(HIGH,KAXIS) /= PERIODIC)) then
     call Driver_abortFlash("Gaurd Cell fill : one edge is periodic, one is not")
  end if

  if((gr_domainBC(LOW,IAXIS) /= PERIODIC).and.(gr_domainBC(HIGH,IAXIS) == PERIODIC)) then
     call Driver_abortFlash("Gaurd Cell fill : one edge is periodic, one is not")
  end if

  if((gr_domainBC(LOW,JAXIS) /= PERIODIC).and.(gr_domainBC(HIGH,JAXIS) == PERIODIC)) then
     call Driver_abortFlash("Gaurd Cell fill : one edge is periodic, one is not")
  end if


  if((gr_domainBC(LOW,KAXIS) /= PERIODIC).and.(gr_domainBC(HIGH,KAXIS) == PERIODIC)) then
     call Driver_abortFlash("Gaurd Cell fill : one edge is periodic, one is not")
  end if

#endif

  
  localGridDataStruct=gridDataStruct
  if((gridDataStruct==CENTER_FACES).and.(NFACE_VARS==0))&
       localGridDataStruct=CENTER
  
  select case(localGridDataStruct)
  case (CENTER)
     beginDataType=CENTER_DATATYPE
     endDataTYPE=CENTER_DATATYPE
  case (CENTER_FACES)
     beginDataType=CENTER_DATATYPE
     if(NDIM==1)endDataType=FACEX_DATATYPE
     if(NDIM==2)endDataType=FACEY_DATATYPE
     if(NDIM==3)endDataType=FACEZ_DATATYPE
  case(FACES)
     beginDataType=FACEX_DATATYPE
     if(NDIM==1)endDataType=FACEX_DATATYPE
     if(NDIM==2)endDataType=FACEY_DATATYPE
     if(NDIM==3)endDataType=FACEZ_DATATYPE
  case(FACEX)
     beginDataType=FACEX_DATATYPE
     endDataType=FACEX_DATATYPE
  case(FACEY)
     if(NDIM>1) then
        beginDataType=FACEY_DATATYPE
        endDataType=FACEY_DATATYPE
     else
        if(gr_meshMe==MASTER_PE)print*,'warning: trying to fill face along Y for 1D problem'
     end if
  case(FACEZ)
     If(NDIM>2) then
        beginDataType=FACEZ_DATATYPE
        endDataType=FACEZ_DATATYPE
     else
        if(gr_meshMe==MASTER_PE)print*,'warning: trying to fill face along Z for 2D problem'
     end if
  end select
  

  do i = beginDataType,endDataType
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,gr_gridDataStruct(i))
     recvLeft(:,:)=1
     sendLeft(:,:)=1
     recvRight(:,:)=1
     sendRight(:,:)=1 !! do a default initialization of all starting points

     do j = 1,NDIM
        sendRight(j, j) = blkLimits(HIGH,j)-blkLimits(LOW,j)+2-gr_offset(i,j)
        
        !index of first interior cell to be sent to be GC on block to left
        sendLeft(j, j) = blkLimits(LOW,j)+gr_offset(i,j)
        
        recvLeft(j, j) = blkLimits(HIGH,j)+1 !recv index of GC 
        call gr_shiftData(gr_axisComm(j), gr_exch(i,j), &
             sendRight(j,:), sendLeft(j,:), &
             recvRight(j,:),recvLeft(j,:),gr_gridDataStruct(i))
     end do
  end do
  if(.not.gr_allPeriodic) then
     do n = 0,NDIM-1
        axis = gr_bndOrder(NDIM-n)
        call gr_bcApplyToAllBlks(axis,isWork)
     end do
  end if
  gr_justExchangedGC = .true.
end subroutine Grid_fillGuardCells

