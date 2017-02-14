!!****if* source/Grid/GridMain/Grid_putPointData
!!
!! NAME
!!  Grid_putPointData
!!
!! SYNOPSIS
!!
!!  Grid_putPointData(integer(IN) :: blockID,
!!                    integer(IN) :: gridDataStruct,
!!                    integer(IN) :: structIndex,
!!                    integer(IN) :: beginCount, 
!!                    integer(IN) :: position(MDIM),
!!                    real(IN)   :: datablock)
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Puts a point of simulation data for a single variable into
!!  specified data structure
!!  
!!  The user is
!!  allowed to specify if index counting should begin at the exterior edge
!!  of the block, (that is including guardcells)
!!  or the interior edge of the block 
!!
!!  
!! ARGUMENTS 
!!
!!  blockID : the local blockid
!!
!!  gridDataStruct : integer value specifying data structure. 
!!                   The options are defined in constants.h and they are :
!!                   CENTER cell centered variables
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along KAXIS
!!                   WORK   single, cell centered variable, valid only
!!                          for paramesh
!!                   SCRATCH scratch space that can fit cell and face centered variables
!!                   SCRATCH_CTR scratch space for cell centered variables
!!                   SCRATCH_FACEX scratch space facex variables
!!                   SCRATCH_FACEY scratch space facey variables
!!                   SCRATCH_FACEZ scratch space facez variables
!!
!!  structIndex : integer value that specifies which variable to put into storage.
!!             for example: DENS_VAR, PRES_VAR as defined in Flash.h 
!!
!!  
!!
!!  beginCount : tells the routine where to start indice counting.  beginCount can
!!               be set to INTERIOR or EXTERIOR.  If INTERIOR is specified
!!               guardcell indices are not included and index 1 is the first interior cell. 
!!               If EXTERIOR is specified
!!               the first index, 1, is the left most guardcell.  See example
!!               below for more explanation.  (For most of the FLASH architecture code,
!!               we use EXTERIOR.  Some physics routines, however, find it helpful 
!!               only to work on the internal parts of the blocks (without
!!               guardcells) and wish to keep loop indicies  
!!               going from 1 to NXB without having to worry about finding 
!!               the correct offset for the number of guardcells.) 
!!               (INTERIOR and EXTERIOR are defined in constants.h)
!!
!!
!!  position(MDIM):
!!           specifies the point to return
!!   
!!           position(1) = i
!!           position(2) = j
!!           position(3) = k
!!
!!           If a problem is only 2d, position(3) is ignored.  For 1d problems
!!           position(2) and position(3) are ignored.
!!
!!
!!  datablock : a real value containing the data
!!
!!
!!
!!
!! EXAMPLE  
!!  
!!    Here is a 3d block example putting a point of data, 
!!    for each block on a local processor.  We will put the point at 
!!    position i=5, j=6 and k=7.  beginCount is set to EXTERIOR
!!    meaning that the first cell of the block including guardcells is 
!!    index = 1.  If in this example, NGUARD = 4 then position 5 is the
!!    first interior cell in a dimension.
!!
!! 
!!    (very hard to draw this, especially the j axis. 
!!     picture is just meant to help show where 
!!     counting begins when "beginCount" is set to EXTERIOR)
!!    
!!     j
!!    
!!     1         - - - - - - - -  
!!     2         - - - - - - - - 
!!     3         - - - - - - - - 
!!     4         - - - - - - - - 
!!     5 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     6 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     7 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     8 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     9 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    10 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    11 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    12 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    13         - - - - - - - - 
!!    14         - - - - - - - - 
!!    15         - - - - - - - - 
!!    16         - - - - - - - - 
!! i     1 2 3 4 5 6 7 8 9 10111213141516 
!!
!!
!!    #include "Flash.h"
!!    #include "constants.h"
!!
!!    ...
!!       
!!      integer ::    position(MDIM)
!!      integer ::    blockID
!!      real    ::    dataBlock
!!
!!          position(1) = 5    
!!          position(2) = 6
!!          position(3) = 7
!!
!!          
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_putPointData(blockID, CENTER, DENS_VAR, EXTERIOR, &
!!                               position, dataBlock)
!!  
!!          end do
!!
!!
!!  
!!    In this 2d block example we will put a point of data for the pressure variable
!!    for each block on a local processor.
!!    beginCount is set to INTERIOR, meaning that all the position indices
!!    will start where index 1 is the first interior cell of the block.
!!    In this example we will put a point where i=4, j=5
!!
!!    (hard to draw, but this is the idea, stars (*) are the cells to return
!!     notice the where indice counting starts when beginCount is set to INTERIOR)
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!     j      - - - - - - - - 
!!     1 ----|-|-|-|-|-|-|-|-|----
!!     2 ----|-|-|-|-|-|-|-|-|----
!!     3 ----|-|-|-|-|-|-|-|-|----
!!     4 ----|-|-|-|-|-|-|-|-|----
!!     5 ----|-|-|-|*|-|-|-|-|----
!!     6 ----|-|-|-|-|-|-|-|-|----
!!     7 ----|-|-|-|-|-|-|-|-|----
!!     8 ----|-|-|-|-|-|-|-|-|----
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!         i  1-2-3-4 5-6-7-8 
!!
!!
!! 
!!    #include "Flash.h"
!!    #include "constants.h"
!!
!!    ...
!!       
!!      integer ::    position(MDIM)
!!      integer ::    blockID
!!      real    ::    dataBlock
!!       
!!          position(1) = 4    
!!          position(2) = 5
!!          position(3) = 1 !ignored since only 2d problem
!!
!!
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_putPointData(blockID, CENTER, PRES_VAR, INTERIOR, &
!!                               position, dataBlock)
!!  
!!          end do
!!
!!
!!
!!
!!
!!
!! 
!!
!!***

!!REORDER(5): unk, facevar[xyz]
!!REORDER(4): solnData

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_putPointData(blockid, gridDataStruct, structIndex, beginCount, &
     position, datablock)

  use Grid_data, ONLY : gr_iguard, gr_jguard, gr_kguard
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr
  use gr_interface, ONLY : gr_getInteriorBlkPtr,gr_releaseInteriorBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"
#ifdef Grid_releaseBlkPtr
! disabling per-block drift logging because this routine gets called too much
! and clogs the drift logs with too much noise.  see: drift
#undef Grid_releaseBlkPtr
#endif


  integer, intent(IN) :: blockid, structIndex, beginCount, gridDataStruct
  integer, dimension(MDIM), intent(IN) :: position
  real, intent(IN) :: datablock
  real, pointer, dimension(:,:,:,:) :: solnData



  integer ::  var, i, j, k, ii
  integer,dimension(MDIM) :: begOffset,dataLen
  integer :: imax, jmax, kmax

  logical :: isget
  logical :: getIntPtr

#ifdef DEBUG_GRID
  isget=.true.
  call gr_checkDataType(blockID,gridDataStruct,imax,jmax,kmax,isget)


  !verify we have a valid blockid
  if((blockid<1).or.(blockid>MAXBLOCKS)) then
     print*,' Put Data : invalide blockid '
     call Driver_abortFlash("Put Data : invalid blockid ")
  end if


  
  !verify beginCount is set to a valid value
  if((beginCount /= INTERIOR) .and. (beginCount /= EXTERIOR)) then
     print *, "Grid_putPointData: beginCount set to improper value"
     print *, "beginCount must = INTERIOR or EXTERIOR (defined in constants.h)"
     call Driver_abortFlash("beginCount must = INTERIOR or EXTERIOR (defined in constants.h)")
  end if



  !verify that indicies aren't too big or too small for the block
  if(beginCount == EXTERIOR) then
    
     if (position(1) > imax) then
        call Driver_abortFlash("Grid_putPointData position(1) index larger than block")
     end if

     if ((NDIM > 1) .and. (position(2) > jmax)) then
        call Driver_abortFlash("Grid_putPointData position(2) index larger than block")
     end if
    
     if ((NDIM > 2) .and. (position(3) > kmax)) then
        call Driver_abortFlash("Grid_putPointData position(3) index larger than block")
     end if
    
     if (position(1) < 1) then
        print*,'position is ',position(1)
        call Driver_abortFlash("Grid_putPointData position(1) index smaller than 1")
     end if

     if ((NDIM > 1) .and. (position(2) < 1)) then
        call Driver_abortFlash("Grid_putPointData position(2) index smaller than 1")
     end if
    
     if ((NDIM > 2) .and. (position(3) < 1)) then
        call Driver_abortFlash("Grid_putPointData position(3) index smaller than 1")
     end if
        
  else !beginCount == INTERIOR

     if ((position(1) + gr_iguard -1) > imax) then
        call Driver_abortFlash("Grid_putPointData position(1) index larger than block")
     end if

     if ((NDIM > 1) .and. ((position(2) + gr_jguard -1) > jmax)) then
        call Driver_abortFlash("Grid_putPointData position(2) index larger than block")
     end if
    
     if ((NDIM > 2) .and. ((position(3) + gr_kguard -1) > kmax)) then
        call Driver_abortFlash("Grid_putPointData position(3) index larger than block")
     end if
    
     if (position(1) < 1) then
        call Driver_abortFlash("Grid_putPointData position(1) index smaller than 1")
     end if

     if ((NDIM > 1) .and. (position(2) < 1)) then
        call Driver_abortFlash("Grid_putPointData position(2) index smaller than 1")
     end if
    
     if ((NDIM > 2) .and. (position(3) < 1)) then
        call Driver_abortFlash("Grid_putPointData position(3) index smaller than 1")
     end if

  end if
  
#endif

  dataLen=0
  call gr_getDataOffsets(blockID,gridDataStruct,position,dataLen,beginCount,begOffset,getIntPtr)


  k = 1
  if(NDIM > 2)k = position(KAXIS) + begOffset(KAXIS) 
  
  j = 1
  if(NDIM > 1)j = position(JAXIS) + begOffset(JAXIS)
  
  i = position(IAXIS) + begOffset(IAXIS)
  
  if(getIntPtr) then
     call gr_getInteriorBlkPtr(blockID,solnData,gridDataStruct)
     solnData(structIndex,i,j,k) = datablock
     call gr_releaseInteriorBlkPtr(blockID,solnData,gridDataStruct)
  else
     call Grid_getBlkPtr(blockID,solnData,gridDataStruct)
!!$     if(gridDataStruct==SCRATCH) then
!!$        solnData(structIndex,i,j,k) = datablock
!!$     else
!!$     end if
     solnData(structIndex,i,j,k) = datablock
     call Grid_releaseBlkPtr(blockID,solnData,gridDataStruct)
  end if
  return
end subroutine Grid_putPointData
