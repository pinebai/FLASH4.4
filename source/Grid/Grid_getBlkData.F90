!!****f* source/Grid/Grid_getBlkData
!!
!! NAME
!!  Grid_getBlkData
!!
!! SYNOPSIS
!!
!!  Grid_getBlkData(integer(IN) :: blockID,
!!                  integer(IN) :: gridDataStruct,
!!                  integer(IN) :: structIndex,
!!                  integer(IN) :: beginCount, 
!!                  integer(IN) :: startingPos(MDIM),
!!                  real(OUT)   :: datablock(dataSize(1),dataSize(2), dataSize(3)),
!!                  integer(IN) :: dataSize(3))
!!  
!! DESCRIPTION 
!!  
!!  This routine allows the user to get an entire block or any contiguous
!!  part of a block of data depending on the arguments passed.
!!  The data could be either a single variable in one of the grid data 
!!  structures such as cell centered/face centered/scratch, or it could
!!  be quantities derived from the size of the cell, such as cell volume/
!!  face area for all the cells in the specified portion of the block
!!
!!  The user is also allowed to specify if index counting 
!!  should begin at the exterior edge of the block (that is including 
!!  guardcells), or the interior edge of the block.
!!
!!  For 3 dimensional problems, a 3d block is returned.
!!  For 2 dimensional problems, a 2d block/plane is returned.
!!  For 1 dimensional problems, a 1d block/row is returned.
!!  
!! ARGUMENTS 
!!
!!  blockID : the local blockid
!!
!!
!!  gridDataStruct : integer value specifying the type of data desired.
!!             Valid options are either one of the Grid data structures,
!!             or one of the derived quantities such as the cell volume
!!             or the cell area.
!!
!!             The options are defined in constants.h and they are :
!!                   CENTER cell centered variables
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!                   WORK   single, cell centered variable, valid only
!!                          for paramesh
!!                   SCRATCH scratch space that can fit cell and face centered variables
!!                   SCRATCH_CTR scratch space for cell centered variables
!!                   SCRATCH_FACEX scratch space facex variables
!!                   SCRATCH_FACEY scratch space facey variables
!!                   SCRATCH_FACEZ scratch space facez variables
!!                   CELL_VOLUME volumes of specified cells 
!!                   CELL_FACEAREA face area for the specified cells on 
!!                                 specified face 
!!                   
!!  structIndex : integer value that provides index into the specified data
!!                structure. When gridDataStruct is one of the grid data structures
!!                structIndex translates to a specific variable in that 
!!                data structure, for example: DENS_VAR or PRES_VAR (define in
!!                Flash.h) if gridDataStruct = CENTER.
!!
!!                For gridDataStruct=CELL_VOLUME it has no meaning
!!
!!                For gridDataStruct = CELL_FACEAREA, it can take one of 
!!                (ILO_FACE,IHI_FACE,JLO_FACE,JHI_FACE,KLO_FACE,KHI_FACE)
!!                defined in constants.h
!!
!!
!!  beginCount : tells the routine where to start index counting.  beginCount can
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
!!  startingPos(MDIM):
!!           specifies the starting position of the block being fetched
!!   
!!           startingPos(1) = i
!!           startingPos(2) = j
!!           startingPos(3) = k
!!
!!           If a problem is 2 dimensions startingPos(3) is irrelevant and
!!           ignored.  For 1d problems startingPos(2) and startingPos(3) 
!!           are ignored.
!!
!!
!!  datablock : a real 3 dimensional array containing the data returned
!!              The dimensions for datablock are 
!!              datablock(irange, jrange, krange) 
!!              Various compilers require the dimensions of
!!              datablock to be specified explicitly.  They are defined by 
!!              the next argument "dataSize".  
!!
!!
!!  dataSize : an integer array specifying the dimensions for datablock
!!          
!!          dataSize(1) holds the number of cells returned in the i direction
!!
!!          dataSize(2) holds the number of cells returned in the j direction
!!                      if 1 dim problem set dataSize(2) = 1
!!
!!          dataSize(3) holds the number of cells returned in the k direction
!!                      if 1 or 2 dim problem set dataSize(3) = 1
!!
!!
!!
!! EXAMPLE  
!!
!!   
!!    EXAMPLE 1:  
!!    Here is a 3d block example getting the entire block of cell centered 
!!    data, including guardcells,
!!    for each block on a local processor.  beginCount is set to EXTERIOR
!!    meaning that the first cell of the block including guardcells is 
!!    index = 1.  If in this example, NGUARD = 4 then position 5 is the
!!    first interior cell in a dimension.
!!
!!    (Can't really draw 3d for this example, but this
!!     picture is meant to help show where 
!!     counting begins when "beginCount" is set to EXTERIOR)
!!    
!!     j
!!    
!!    16         - - - - - - - -  
!!    15         - - - - - - - - 
!!    14         - - - - - - - - 
!!    13         - - - - - - - - 
!!    12 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    11 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    10 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     9 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     8 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     7 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     6 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     5 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     4         - - - - - - - - 
!!     3         - - - - - - - - 
!!     2         - - - - - - - - 
!!     1         - - - - - - - - 
!! i     1 2 3 4 5 6 7 8 9 10111213141516 
!! 
!! 
!!    #include "Flash.h"
!!    #include "constants.h"
!!
!!    ...
!!       
!!      integer ::    startingPos(MDIM)
!!      integer ::    dataSize(3)
!!      integer ::    blockID
!!      real    ::    dataBlock(:,:,:)
!!       
!!          startingPos(1) = 1 !getting the entire extent of dim   
!!          startingPos(2) = 1 !getting the entire extent of dim
!!          startingPos(3) = 1 !getting the entire extent of dim
!!
!!          dataSize(1) = blkLimitsGC(IAXIS) !This is equivalent to NXB + 2*NGUARD 
!!                                           !in FIXEDBLOCKSIZE mode
!!          dataSize(2) = blkLimitsGC(JAXIS) !This is equivalent to NYB + 2*NGUARD
!!                                           !in FIXEDBLOCKSIZE mode
!!          dataSize(3) = blkLimitsGC(KAXIS) !This is equivalent to NZB + 2*NGUARD
!!                                           !in FIXEDBLOCKSIZE mode
!!
!!          allocate(datablock(dataSize(1), dataSize(2), dataSize(3)))
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_getBlkData(blockID, CENTER, DENS_VAR, EXTERIOR, &
!!                               startingPos, dataBlock, dataSize)
!!  
!!          end do
!!
!!
!!
!!    EXAMPLE 2:
!!    In this 2d block example we will get cell volume of only a portion 
!!    of the interior part of a block for each block on the local processor.
!!    beginCount is set to INTERIOR, meaning that all the startPos indices
!!    will start where index 1 is the first interior cell of the block.
!!    In this example we will get only the last 2 columns of interior cells.
!!    i start position = 7, j start position = 1.  (Because this is only a 2d
!!    problem Grid_getPlaneData could also be used to get these cells.)
!!
!!    (Hard to draw, but this is the idea, stars (*) are the cells to return.
!!     Notice where index counting starts when beginCount is set to INTERIOR)
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!     j      - - - - - - - - 
!!     8 ----|-|-|-|-|-|-|*|*|----
!!     7 ----|-|-|-|-|-|-|*|*|----
!!     6 ----|-|-|-|-|-|-|*|*|----
!!     5 ----|-|-|-|-|-|-|*|*|----
!!     4 ----|-|-|-|-|-|-|*|*|----
!!     3 ----|-|-|-|-|-|-|*|*|----
!!     2 ----|-|-|-|-|-|-|*|*|----
!!     1 ----|-|-|-|-|-|-|*|*|----
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!         i  1-2-3-4 5-6-7-8 
!!
!! 
!!    #include "Flash.h"
!!    #include "constants.h"
!!
!!    ...
!!       
!!      integer ::    startingPos(MDIM)
!!      integer ::    dataSize(3)
!!      integer ::    blockID
!!      real    ::    dataBlock(:,:,:)
!!       
!!          startingPos(1) = 7    
!!          startingPos(2) = 1
!!          startingPos(3) = 1 !value ignored since only 2d
!!
!!          dataSize(1) = 2 !just getting 2 cells in the i dir
!!          dataSize(2) = 8 !just getting 8 cells in the j dir
!!          dataSize(3) = 1 !needs to be set to 1, since only 2d problem
!!
!!          allocate(datablock(dataSize(1), dataSize(2), dataSize(3)))
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_getBlkData(blockID, CELL_VOLUME, 0, INTERIOR, &
!!                               startingPos, dataBlock, dataSize)
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


#define DEBUG_GRID


subroutine Grid_getBlkData(blockID, gridDataStruct, structIndex, beginCount, &
     startingPos, datablock, dataSize)

  implicit none
#include "constants.h"
  integer, intent(in) :: blockID, structIndex, beginCount, gridDataStruct
  integer, dimension(MDIM), intent(in) :: startingPos
  integer, dimension(3), intent(in) :: dataSize
  real, dimension(datasize(1), dataSize(2), dataSize(3)),intent(out) :: datablock
  datablock = 0.0
  return

end subroutine Grid_getBlkData

